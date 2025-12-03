use rust_htslib::bam::HeaderView;
use rust_htslib::bam::{Read, Reader, header, Record, Format::Bam, record::Aux, Writer};
use std::{collections::HashMap, path::Path};
use std::sync::mpsc;
use anyhow::{bail, Result as anyResult};
use seq::{SequenceTable, SequenceLookup};
use crate::metrics::Metrics;
use crate::seq;
use crate::util;

/// Represents a thread for reading BAM records
/// Designed to handle BAM file reading in a separate thread
/// It provides a mechanism to receive batches of BAM records through a channel
///
/// # Fields:
/// - thread: The handle to the thread that performs the BAM file reading
/// - rx: A receiver channel used to retrieve batches of BAM records
/// 
/// # Arguments:
/// - bam_reader: An instance of Reader from the rust_htslib library, which is used to read BAM files
/// - read_processing_batch_size: Number of records to process in a single batch
/// - buffer_size: Size of the channel buffer, which determines how many batches can be queued before blocking
pub struct BamReaderThread {
    thread: std::thread::JoinHandle<()>,
    pub rx: mpsc::Receiver<Vec<Record>>,
}

impl BamReaderThread {
    pub fn new(mut bam_reader: Reader, read_processing_batch_size: usize, buffer_size: usize) -> Self {
        let (tx, rx) = mpsc::sync_channel(buffer_size);
        let mut batch = Vec::with_capacity(read_processing_batch_size);
        let thread = std::thread::spawn(move || {
            let mut record = Record::new(); // Reuse same bam record
            // Using records iterator yields memory corruption issue https://github.com/rust-bio/rust-htslib/issues/479
            while let Some(r) = bam_reader.read(&mut record) {
                r.expect("Failed to parse record");
                batch.push(record.clone());
                // take yields the same memory corruption issue
                //batch.push(std::mem::take(&mut record));
                if batch.len() == read_processing_batch_size {
                    // take ensures that batch is cleared after sending, thus making it reusable
                    if tx.send(std::mem::take(&mut batch)).is_err() {
                        // Error in the receiver thread; shutting down
                        return;
                    }
                }
            }
            if !batch.is_empty() {
                // Can ignore error, since thread is done
                let _ = tx.send(batch);
            }
        });
        Self {
            thread,
            rx,
        }
    }

    pub fn finish(self) {
        drop(self.rx);
        self.thread.join().expect("Error closing BamReaderThread");
    }
}

/// Interface for reading BAM files in batches
/// Initialize the BAM file reader and start a background thread to process records in batches
///
/// # Fields:
/// - bam_reader_thread: Instance of BamReaderThread that handles reading BAM records in a separate thread
/// - header: Header of input BAM file
/// 
/// # Arguments:
/// - bam: Path to input BAM file
/// - read_processing_batch_size: Number of records to process in a single batch
/// - buffer_size: Size of the channel buffer, which determines how many batches can be queued before blocking
pub struct BamReader {
    pub bam_reader_thread: BamReaderThread,
    header: HeaderView,
}

impl BamReader {
    pub fn new(bam: &Path, read_processing_batch_size: usize, buffer_size: usize) -> Self {
        let bam_reader = Reader::from_path(bam).expect("Failed to open BAM file");
        let header = bam_reader.header().clone();
        let bam_reader_thread = BamReaderThread::new(bam_reader, read_processing_batch_size, buffer_size);
        
        Self {
            bam_reader_thread,
            header,
        }
    }
    pub fn get_next_record_batch(&self) -> Option<Vec<Record>> {
        // Returns None when EOF
        self.bam_reader_thread.rx.recv().ok()
    }
}

/// Interface for writing BAM files in batches
/// Initialize the BAM file writer and start a background thread to write records in batches
///
/// # Fields:
/// - bam_writer_thread: Instance of BamWriterThread that handles writing BAM records in a separate thread
/// 
/// # Arguments:
/// - bam: Path to input BAM file
/// - out_dir: Output directory where the tagged BAM file will be written
/// - bam_reader: Reference to an instance of BamReader, which provides the header for the BAM file
/// - buffer_size: Size of the channel buffer, which determines how many batches can be queued before blocking
pub struct BamWriter {
    pub bam_writer_thread: BamWriterThread,
}

impl BamWriter {
    pub fn new(bam: &Path, out_dir: &Path, bam_reader: &BamReader, buffer_size: usize) -> Self {
        let mut bam_writer = rust_htslib::bam::Writer::from_path(
            out_dir.join(bam.file_name().unwrap()).with_extension("tagged.bam"),
            &header::Header::from_template(&bam_reader.header),
            Bam,
        ).expect("Failed to create BAM writer");
        bam_writer.set_threads(4).unwrap();
        let bam_writer_thread = BamWriterThread::new(bam_writer, buffer_size);
        Self {
            bam_writer_thread,
        }
    }
}

/// Represents a thread for writing BAM records
/// Designed to handle BAM file writing in a separate thread
/// It provides a mechanism to send batches of BAM records through a channel
///
/// # Fields:
/// - thread: The handle to the thread that performs the BAM file writing
/// - tx: A sender channel used to send batches of BAM records
/// 
/// # Arguments:
/// - bam_writer: An instance of Writer from the rust_htslib library, which is used to write BAM files
/// - buffer_size: Size of the channel buffer, which determines how many batches can be queued before blocking
pub struct BamWriterThread {
    thread: std::thread::JoinHandle<()>,
    tx: mpsc::SyncSender<Vec<Record>>,
}

impl BamWriterThread {
    pub fn new(mut bam_writer: Writer, buffer_size: usize) -> Self {
        let (tx, rx) = mpsc::sync_channel::<Vec<Record>>(buffer_size);
        let thread = std::thread::spawn(move || {
            while let Ok(batch_of_records) = rx.recv() {
                for record in batch_of_records {
                    bam_writer.write(&record).expect("Failed to write BAM record");
                }
            }
        });
        Self { thread, tx }
    }

    pub fn write(&self, bam_batch: Vec<Record>) {
        self.tx.send(bam_batch).expect("Unable to send batch of BAM records");
    }

    pub fn finish(self) {
        drop(self.tx);
        self.thread.join().expect("Error closing BamWriterThread");
    }
}

/// Interface for creating a tagged BAM file
/// This struct provides members that enable reading from an input BAM file,
///  processing records to add tags based on a target sequence map, and write those records to a new BAM file
/// 
/// # Fields:
/// - bam_reader: Instance of BamReader that reads records from the input BAM file
/// - bam_writer: Instance of BamWriter that writes tagged records to the output BAM file
/// - target_map: Used for looking up target sequences
/// - out_tag: The tag to be added to the BAM records
/// 
/// # Arguments:
/// - bam: Path to input BAM file
/// - target_map: Used for looking up target sequences
/// - out_tag: The tag to be added to the BAM records
/// - out_dir: Output directory where the tagged BAM file will be written
pub struct CreateTaggedBam<'a> {
    pub bam_reader: BamReader,
    pub bam_writer: BamWriter,
    pub target_map: SequenceTable,
    pub out_tag: &'a [u8],
}
impl<'a> CreateTaggedBam<'a> {
    pub fn new(
        bam: &Path,
        target_map: SequenceTable,
        out_tag: &'a str,
        out_dir: &'a Path,
        read_processing_batch_size: usize,
        buffer_size: usize
    ) -> Self {
        let bam_reader = BamReader::new(bam, read_processing_batch_size, buffer_size);
        let bam_writer = BamWriter::new(bam, out_dir, &bam_reader, buffer_size);
        Self {
            bam_reader,
            bam_writer,
            target_map,
            out_tag: out_tag.as_bytes(),
        }
    }
    // Search for target in bam record based on most frequent start position
    pub fn compute_tag_to_add_to_bam_record(
        &self,
        record_to_write: &mut Record,
        most_freq_start_pos: usize,
        seq: &mut Vec<u8>,
        metrics: &mut Metrics,
    ) {
        seq.extend(record_to_write.seq().as_bytes());
        match self.target_map
            .lookup(&seq[most_freq_start_pos..most_freq_start_pos + self.target_map.min_length])
        {
            SequenceLookup::Exact(name) => {
                self.push_tag(name, record_to_write, self.out_tag).expect("Failed to add tag to BAM record");
                metrics.exact_count += 1;
            }
            SequenceLookup::ErrorOf(name) => {
                self.push_tag(name, record_to_write, self.out_tag).expect("Failed to add tag to BAM record");
                metrics.mismatch_count += 1;
            }
            _ => {}
        }
        seq.clear();
    }
    
    fn push_tag(&self, name: &str, record_to_write: &mut Record, out_tag: &[u8]) -> anyResult<()> {
        let tag = Aux::String(name);
        if let Err(e) = record_to_write.push_aux(out_tag, tag) {
            bail!("Failed to add tag to BAM record: {}", e);
        }
        Ok(())
    }
}

/// Interface for making a hashmap of target start positions to their frequencies in the BAM file
/// This struct provides members that enable reading from an input BAM file,
///  processing records to find the starting position of a target in the record, and counting its frequency
/// 
/// # Fields:
/// - bam_reader: Instance of BamReader that reads records from the input BAM file
/// - target_map: Used for looking up target sequences
/// - num_reads_to_find_start_pos: Number of reads to process to build the frequency hashmap
/// 
/// # Arguments:
/// - bam: Path to input BAM file
/// - target_map: Used for looking up target sequences
/// - read_processing_batch_size: Number of records to process in a single batch
/// - buffer_size: Size of the channel buffer, which determines how many batches can be queued before blocking
/// - num_reads_to_find_start_pos: Number of reads to process to build the frequency hashmap
pub struct CreateFrequencyHashmap {
    pub bam_reader: BamReader,
    pub target_map: SequenceTable,
    pub num_reads_to_find_start_pos: usize,
}

impl CreateFrequencyHashmap {
    pub fn new(
        bam: &Path,
        target_map: SequenceTable,
        read_processing_batch_size: usize,
        buffer_size: usize,
        num_reads_to_find_start_pos: usize
    ) -> Self {
        let bam_reader = BamReader::new(bam, read_processing_batch_size, buffer_size);
        Self {
            bam_reader,
            target_map,
            num_reads_to_find_start_pos,
        }
    }
    // Make hashmap of target start positions to their frequencies in the BAM file
    pub fn construct_target_start_pos_to_frequency_hashmap(&self) -> HashMap<usize, usize> {
        let mut target_position_frequency: HashMap<usize, usize> = HashMap::new();
        while let Some(bam_record_batch) = self.bam_reader.get_next_record_batch() {
            let mut read_count = 0; // Counter to track the number of input reads processed
            for record in bam_record_batch {
                let bam_record_seq = record.seq().as_bytes();
                let target_len = self.target_map.min_length;
                let record_len = bam_record_seq.len();
                if record_len > target_len { // Prevent out of bounds error
                    for i in 0..=record_len - target_len {
                        let subslice = &bam_record_seq[i..i + target_len];
                        match self.target_map.lookup(subslice) {
                            SequenceLookup::Exact(_) => {
                                // Assign score of 3 to exact matches to prioritize them
                                util::increment_frequency_of_target_start_pos(&mut target_position_frequency, i, 3);
                            }
                            SequenceLookup::ErrorOf(_) => {
                                // Assign score of 1 to mismatches
                                util::increment_frequency_of_target_start_pos(&mut target_position_frequency, i, 1);
                            }
                            _ => {}
                        }
                    }
                }
                read_count += 1;
                if read_count == self.num_reads_to_find_start_pos {
                    return target_position_frequency;
                }
            }
        }
        // This return is only triggered when input bam has less than @num_reads_to_find_start_pos reads
        target_position_frequency
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Sequence;
    use rust_htslib::bam::header::{HeaderRecord, Header};
    use tempfile::NamedTempFile;
    
    fn create_test_bam_writer() -> (Writer, NamedTempFile) {
        // Create a temp bam file just for this test
        let tmpfile = NamedTempFile::new().unwrap();
        let hd_record = HeaderRecord::new("CO\ttest".as_bytes());
        let mut header = Header::new();
        header.push_record(&hd_record);
        let path = tmpfile.path();
        let writer = Writer::from_path(&path, &header, Bam).unwrap();
        // Return tmpfile object so that it does not go out of scope
        // If it goes out of scope, the file is deleted
        (writer, tmpfile)
    }
    fn create_test_record(read_name: &str, seq: &str) -> Record {
        let mut record = Record::new();
        let quality_scores: Vec<u8> = vec![b'I'; seq.len()];
        record.set(read_name.as_bytes(), None, seq.as_bytes(), &quality_scores);
        record.set_tid(-1);
        record.set_mtid(-1);
        record.set_pos(-1);
        record.set_mpos(-1);
        record.set_unmapped();
        record
    }
    #[test]
    fn test_construct_target_start_pos_to_frequency_hashmap() {
        let mut seq_table = SequenceTable::default();
        let seq1 = Sequence::new("ACGT".to_string()).unwrap();
        seq_table.add_seq(&seq1, "target1");
        seq_table.min_length = seq1.len();
        
        let (mut bam_writer, tmpfile) = create_test_bam_writer();
        
        bam_writer.write(&create_test_record("read1", "ACGTACGT")).unwrap();
        bam_writer.write(&create_test_record("read2", "ACGTGCGT")).unwrap();
        // Close the writer to flush the records to the file
        // Need to manually do this since bam_writer does not go out of scope in this function
        //  and we need to read the file later in this function
        drop(bam_writer);
        
        let path = tmpfile.path();
        let create_frequency_hashmap = CreateFrequencyHashmap::new(
            &path,
            seq_table,
            1,
            1,
            2,
        );

        let frequency_map = create_frequency_hashmap.construct_target_start_pos_to_frequency_hashmap();
        // Two exact matches at the 0th position mean that the frequency is 6 (3 for each exact match)
        assert_eq!(frequency_map.get(&0), Some(6).as_ref());

    }
}
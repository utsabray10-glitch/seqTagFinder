#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]

use clap::{App, Arg, ArgMatches};
use anyhow::{Context, Result as anyResult};
use std::path::{Path, PathBuf};
use std::fs;

mod target;
mod seq;
mod bam;
mod util;
mod metrics;

fn main() -> anyResult<()> {
    let config = Config::from_args()?;
    let target_map = target::TargetProcessor::process(&config.whitelist)?;
    run(
        config.bams,
        config.num_reads,
        &config.out_dir,
        &target_map.target_map,
        &config.out_tag,
        config.read_processing_batch_size,
        config.buffer_size
    );
    Ok(())
}

struct Config { 
    bams: Vec<PathBuf>,
    num_reads: usize,
    out_dir: PathBuf,
    whitelist: PathBuf,
    out_tag: String,
    read_processing_batch_size: usize,
    buffer_size: usize,
}

impl Config {
    fn from_args() -> anyResult<Self> { 
        let args = Self::accept_args();
        Self::parse_args(&args)
    }
    fn accept_args() -> ArgMatches<'static> {
        App::new("seqTagFinder")
            .version(clap::crate_version!())
            .author("Utsab Ray <utsab.ray@scale.bio>")
            .about("Detect tags in BAM files")
            .arg(Arg::from_usage("[bams] --bams [FILE1.bam, FILE2.bam]....")
                .help("BAM files to search for sequences in. Whitespace separated list of BAM files.")
                .required(true))
            .arg(Arg::from_usage("--num_reads <NUM> 'Number of reads to look at in each BAM file while determining position of target sequence in read'")
                .default_value("100000"))
            .arg(Arg::from_usage("--out_dir <OUTPUT_DIR> 'Output directory name'")
                .default_value("taggedBams"))
            .arg(Arg::from_usage("--whitelist <TARGET.txt> 'Whitelist file containing sequences to search for in BAM files'")
                .required(true))
            .arg(Arg::from_usage("--tag_in_output_bam <STRING> 'Tag which will have detected target sequences in output BAM files'")
                .default_value("SP"))
            .arg(Arg::from_usage("--read_processing_batch_size <NUM> 'Number of reads to collect in a single batch for processing'")
                .default_value("100"))
            .arg(Arg::from_usage("--buffer_size <NUM> 'Number of batches of reads a thread will collect before sending over the queue'")
                .default_value("10"))
            .get_matches()
    }
    
    fn parse_args(args: &ArgMatches) -> anyResult<Self> {
        let bams: Vec<PathBuf> = args
            .values_of("bams")
            .unwrap()
            .map(PathBuf::from)
            .collect();
        let num_reads: usize = args
            .value_of("num_reads")
            .unwrap()
            .parse::<usize>()
            .context("Invalid number provided for num_reads")?;
        let out_dir: PathBuf = args
            .value_of("out_dir")
            .unwrap()
            .parse::<PathBuf>()
            .context("Invalid output directory provided")?;
        fs::create_dir_all(&out_dir).context("Failed to create output directory")?;
        let whitelist: PathBuf = args
            .value_of("whitelist")
            .unwrap()
            .parse::<PathBuf>()
            .context("Invalid whitelist file provided")?;
        let out_tag: String = args
            .value_of("tag_in_output_bam")
            .unwrap().to_string();
        let read_processing_batch_size = args
            .value_of("read_processing_batch_size")
            .unwrap()
            .parse::<usize>()
            .context("Invalid number provided for read_processing_batch_size")?;
        let buffer_size = args
            .value_of("buffer_size")
            .unwrap()
            .parse::<usize>()
            .context("Invalid number provided for buffer_size")?;
        Ok(Self { bams, num_reads, out_dir, whitelist, out_tag, read_processing_batch_size, buffer_size })
    }
}

pub fn run(
    bams: Vec<PathBuf>,
    num_reads: usize,
    out_dir: &Path,
    target_map: &seq::SequenceTable,
    out_tag: &str,
    read_processing_batch_size: usize,
    buffer_size: usize
) {
    let mut all_metrics: Vec<metrics::Metrics> = Vec::new();
    for bam in bams {
        let most_freq_start_pos_obj = bam::CreateFrequencyHashmap::new(
            &bam,
            target_map.clone(),
            read_processing_batch_size,
            buffer_size,
            num_reads,
        );
        let create_tagged_bam_obj = bam::CreateTaggedBam::new(
            &bam,
            target_map.clone(),
            out_tag,
            out_dir,
            read_processing_batch_size,
            buffer_size,
        );
        let target_position_frequency = most_freq_start_pos_obj.construct_target_start_pos_to_frequency_hashmap();
        let mut metrics = metrics::Metrics::new(target_position_frequency.clone(), bam.clone());
        if let Some(most_freq_start_pos) = util::get_most_frequently_occuring_key(&target_position_frequency) {
            let mut seq = Vec::new();
            while let Some(mut bam_record_batch) = create_tagged_bam_obj.bam_reader.get_next_record_batch() {
                for record in &mut bam_record_batch {
                    metrics.read_count += 1;
                    create_tagged_bam_obj.compute_tag_to_add_to_bam_record(
                        record,
                        most_freq_start_pos,
                        &mut seq,
                        &mut metrics,
                    );
                }
                // take ensures that batch is cleared after sending, thus making it reusable
                let replacement_batch = std::mem::take(&mut bam_record_batch);
                create_tagged_bam_obj.bam_writer.bam_writer_thread.write(replacement_batch);
            }
        } else {
            // If targets are not found in the BAM file, copy the original BAM to the output directory without modification
            fs::copy(&bam, out_dir.join(&bam.file_name().unwrap())).unwrap();
            // Delete the empty tagged BAM file that gets created when BamWriter::new is called
            fs::remove_file(out_dir.join(bam.file_name().unwrap()).with_extension("tagged.bam")).unwrap();
        }
        most_freq_start_pos_obj.bam_reader.bam_reader_thread.finish();
        create_tagged_bam_obj.bam_reader.bam_reader_thread.finish();
        create_tagged_bam_obj.bam_writer.bam_writer_thread.finish();
        all_metrics.push(metrics);
    }
    metrics::write(all_metrics, out_dir).unwrap();
}

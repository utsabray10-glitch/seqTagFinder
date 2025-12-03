# seqTagFinder

Find tags in BAM files

## Setup

### Install Rust
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Confirm cargo is installed
```bash
cargo --version
```

## Execution

The dependencies are mentioned in Cargo.toml. Execute ```bash
cargo build``` to install dependencies. Alternatively execute ```bash
cargo run -- <command line args>``` to execute the tool

## Command Line Arguments

1. `--bams`: BAM files to search for sequences in. Whitespace separated list of BAM files.
2. `--num_reads`: Number of reads to look at in each BAM file while determining position of target sequence in read
3. `--out_dir`: Output directory path that will contain the output BAM files
4. `--whitelist`: Whitelist file containing sequences to search for in BAM files
5. `--tag_in_output_bam`: Tag which will have detected target sequences in output BAM files
6. `--read_processing_batch_size`: Number of reads to collect in a single batch for processing
7. `--buffer_size`: Number of batches of reads a thread will collect before sending over the queue

## Methodology

For each BAM file, and for each target in the whitelist file, an average start position of the target in the read is computed. For this computation only the top @num_reads are looked at. It is necessary to compute the average start position, otherwise we will have to use other algorithms when figuring out whether the target is present in the read or not. For our current data size, this 2 pass approach is efficient enough.

Once the average start position of each target is obtained, one thread opens the BAM file, reads records and sends batches of records over the queue. Then the main thread sends each record from the batch to another thread that finds the target present in the record, based on the average start position. Once the target is found, the record, along with the target is sent to the writer thread, to write to the output BAM file.

## Output

BAM files in @out_dir with each read annotated with @tag_in_output_bam that designates the target found in that read. Also produces a metrics.json file with total read counts and the number of exact vs mismatches found. We allow for 1 bp mismatch.

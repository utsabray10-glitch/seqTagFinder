use std::io::{BufRead, BufReader};
use anyhow::{anyhow, Context, Result as anyResult};
use std::path::Path;
use std::fs;

use crate::seq;
use seq::{Sequence, SequenceTable};

pub struct TargetProcessor {
    pub target_map: SequenceTable,
}
impl TargetProcessor {
    pub fn process(targets: &Path) -> anyResult<Self> {
        let target_map = Self::read_target_whitelist(targets)?;
        Self::trim_seqs_by_len_in_target_map(target_map) 
    }
    fn read_target_whitelist(target_whitelist: &Path) -> anyResult<SequenceTable> {
        let mut target_lookup = SequenceTable::default();
        let file = fs::File::open(target_whitelist)
            .context(anyhow!("Failed to open target whitelist file {:?}", target_whitelist))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let words: Vec<String> = line?.split_whitespace().map(std::string::ToString::to_string).collect();
            if words.len() > 1 {
                for word in &words[1..] {
                    let target_seq = Sequence::new(word.to_string())?;
                    target_lookup.add_seq(&target_seq, &words[0]);
                }
            }
        }
        Ok(target_lookup)
    }
    fn trim_seqs_by_len_in_target_map(untrimmed_target_map: SequenceTable) -> anyResult<Self> {
        let mut target_map: SequenceTable = SequenceTable::default();
        // Trim all sequences to the minimum length of the sequences in the whitelist
        let min_length = untrimmed_target_map.all_whitelist_combinations
            .keys()
            .map(super::seq::Sequence::len)
            .min()
            .ok_or_else(|| anyhow::Error::msg("Whitelist map is empty".to_string()))?;
        for (target_seq, alias) in untrimmed_target_map.all_whitelist_combinations {
            let trimmed_seq = target_seq.seq.into_iter().take(min_length).collect::<Vec<u8>>();
            target_map.all_whitelist_combinations.insert(Sequence { seq: trimmed_seq.clone() }, alias);
        }
        target_map.min_length = min_length;
        Ok(Self { target_map } )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_seqs_by_len_in_target_map() {
        let mut target_map = SequenceTable::default();
        target_map.add_seq(&Sequence::new("ACGT".to_string()).unwrap(), "target1");
        target_map.add_seq(&Sequence::new("AAGTG".to_string()).unwrap(), "target2");
        target_map.add_seq(&Sequence::new("ACG".to_string()).unwrap(), "target3");

        let trimmed_processor = TargetProcessor::trim_seqs_by_len_in_target_map(target_map).unwrap();
        assert_eq!(trimmed_processor.target_map.min_length, 3);
        assert!(trimmed_processor.target_map.all_whitelist_combinations.contains_key(&Sequence::new("AAG".to_string()).unwrap()));
    }
}
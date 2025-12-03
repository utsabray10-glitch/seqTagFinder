use std::{borrow::Borrow, collections::hash_map::Entry};
use std::collections::HashMap;
use anyhow::{bail, Result as anyResult};


// Data structure for matching expected target sequences up to 1 mismatch
/// Used for error correction of target sequences from reads
#[derive(Clone, Default)]
pub struct SequenceTable {
    pub all_whitelist_combinations: HashMap<Sequence, SequenceLookup>,
    pub min_length: usize, // Minimum length of sequences in the hashmap
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SequenceLookup {
    Exact(String), // Exact match to target sequence
    ErrorOf(String), // 1 bp mismatch to target sequence
    Ambiguous, // Multiple possible mismatches to target sequence
    NoMatch, // No match to target sequence
}
impl SequenceTable {
    pub fn add_seq(&mut self, refseq: &Sequence, name: &str) {
        let refseq_seq = &refseq.seq;
        self.all_whitelist_combinations
            .insert(refseq.clone(), SequenceLookup::Exact(name.to_string()));
        // Enter all sequence neighbors (1 mismatch) into hash
        for i in 0..refseq_seq.len() {
            for single_base in b"ACGTN" {
                if *single_base != refseq_seq[i] {
                    let mut one_mismatch = refseq_seq.clone();
                    one_mismatch[i] = *single_base;
                    let one_mismatch_seq = Sequence { seq: one_mismatch };
                    match self.all_whitelist_combinations.entry(one_mismatch_seq) {
                        // If entry does not exist, enter in hashmap
                        Entry::Vacant(e) => {
                            e.insert(SequenceLookup::ErrorOf(name.to_string()));
                        }
                        // If entry exists, and the entry is Exact, keep it
                        // If entry exists, and the entry is ErrorOf, set it to Ambiguous,
                        // as there are multiple possible mismatches
                        Entry::Occupied(mut e) => 
                            if let SequenceLookup::ErrorOf(_) = e.get() {
                                *(e.get_mut()) = SequenceLookup::Ambiguous;
                        }
                    }
                }
            }
        }
    }

    pub fn lookup(&self, seq: &[u8]) -> &SequenceLookup {
        self.all_whitelist_combinations
            .get(seq)
            .unwrap_or(&SequenceLookup::NoMatch)
    }

}
#[derive(Clone, Debug, Default, Hash, PartialEq, Eq)]
pub struct Sequence {
    pub seq: Vec<u8>
}
impl Borrow<[u8]> for Sequence {
    fn borrow(&self) -> &[u8] {
        &self.seq
    }
}
impl Sequence {
    pub fn new(mut seq: String) -> anyResult<Self> {
        if !seq.is_ascii() {
            bail!("Invalid characters in sequence '{}'", seq);
        }
        seq.make_ascii_uppercase();
        if !seq
            .trim_start_matches(&['A', 'C', 'T', 'G', 'N'][..])
            .is_empty()
        {
            bail!("Unknown base in sequence '{}'", seq);
        }
        Ok(Self{seq: seq.into_bytes()})
    }

    pub const fn len(&self) -> usize {
        self.seq.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_seq() {
        let mut seq_table = SequenceTable::default();
        let seq1 = Sequence::new("ACGT".to_string()).unwrap();
        seq_table.add_seq(&seq1, "target1");
        // Input sequence is of length 4. Each base can be replaced with 4 other bases (other than itself)
        // Total combinations = 4*4 = 16 + 1 (the input sequence itself) = 17
        assert_eq!(seq_table.all_whitelist_combinations.len(), 17);
    }

    #[test]
    fn test_lookup() {
        let mut seq_table = SequenceTable::default();
        let seq1 = Sequence::new("ACGT".to_string()).unwrap();
        seq_table.add_seq(&seq1, "target1");
        assert_eq!(seq_table.lookup(b"ACGT"), &SequenceLookup::Exact("target1".to_string()));
        assert_eq!(seq_table.lookup(b"ACGA"), &SequenceLookup::ErrorOf("target1".to_string()));
        assert_eq!(seq_table.lookup(b"ACGTN"), &SequenceLookup::NoMatch);
        assert_eq!(seq_table.lookup(b"AAAA"), &SequenceLookup::NoMatch);
    }
}
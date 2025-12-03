use std::collections::HashMap;

pub fn get_most_frequently_occuring_key(input_hashmap: &HashMap<usize, usize>) -> Option<usize> {
    input_hashmap
        .iter()
        .max_by_key(|&(_, count)| count)
        .map(|(&key, _)| key)
}

pub fn increment_frequency_of_target_start_pos(target_position_frequency: &mut HashMap<usize, usize>, pos: usize, score: usize) {
    target_position_frequency
        .entry(pos)
        .and_modify(|count| *count += score)
        .or_insert(score);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_most_frequently_occuring_key() {
        let mut map = HashMap::new();
        map.insert(1, 2);
        map.insert(2, 3);
        map.insert(3, 1);
        assert_eq!(get_most_frequently_occuring_key(&map), Some(2));
        // When there are multiple keys with the same maximum frequency, it returns the first one encountered (copilot's words, https://doc.rust-lang.org/beta/std/cmp/fn.max_by_key.html is uninformative)
        // That does not seem to be deterministic since the iteration order of a HashMap is not guaranteed to be in insertion order
        // Uncommenting the code below results in a failed test because the result seems to change everytime. Sometimes the test passes, sometimes it fails
        //map.insert(4, 3);
        //assert_eq!(get_most_frequently_occuring_key(&map), Some(4));
    }

    #[test]
    fn test_increment_frequency_of_target_start_pos() {
        let mut map = HashMap::new();
        increment_frequency_of_target_start_pos(&mut map, 1, 5);
        assert_eq!(map.get(&1), Some(&5));
        increment_frequency_of_target_start_pos(&mut map, 1, 3);
        assert_eq!(map.get(&1), Some(&8));
    }
}
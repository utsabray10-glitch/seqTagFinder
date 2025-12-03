use std::path::{Path, PathBuf};
use std::collections::HashMap;

pub struct Metrics {
    pub input_bam: PathBuf,
    pub read_count: u64,
    pub target_position_frequency: HashMap<usize, usize>,
    pub exact_count: u64,
    pub mismatch_count: u64,
}
impl Metrics {
    pub fn new(target_position_frequency: HashMap<usize, usize>, input_bam: PathBuf) -> Self {
        Self {
            input_bam,
            read_count: 0,
            target_position_frequency,
            exact_count: 0,
            mismatch_count: 0,
        }
    }
}
pub fn write(metrics: Vec<Metrics>, out_dir: &Path) -> anyhow::Result<()> {
    let mut all_metrics = json::JsonValue::new_array();
    for metric in metrics {
        let mut metric_json = json::JsonValue::new_object();
        let fname = metric.input_bam.to_str().unwrap();
        metric_json[fname]["read"] = metric.read_count.into();
        metric_json[fname]["exact"] = metric.exact_count.into();
        metric_json[fname]["mismatch"] = metric.mismatch_count.into();
        
        // Convert HashMap to JsonValue
        let mut target_position_frequency_json = json::JsonValue::new_object();
        for (key, value) in &metric.target_position_frequency {
            target_position_frequency_json[key.to_string()] = json::JsonValue::from(*value);
        }
        metric_json[fname]["target_position_frequency"] = target_position_frequency_json;

        all_metrics.push(metric_json)?;
    }
    
    let mut out = std::fs::File::create(out_dir.join("metrics.json"))?;
    all_metrics.write_pretty(&mut out, 4)?;
    Ok(())
}
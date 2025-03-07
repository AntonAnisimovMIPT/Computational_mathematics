use std::fs::{File, create_dir_all};
use std::io::{Write, BufWriter};
use std::path::Path;

pub fn save_residuals_to_csv(residuals: &[f64], filename: &str) {
    
    let path = Path::new(filename);
    if let Some(parent) = path.parent() {
        create_dir_all(parent).expect("Failed to create directories");
    }

    let file = File::create(filename).expect("Unable to create file");
    let mut writer = BufWriter::new(file);

    writeln!(writer, "iteration,residual").expect("Unable to write header");
    for (iteration, residual) in residuals.iter().enumerate() {
        writeln!(writer, "{},{}", iteration, residual).expect("Unable to write data");
    }
}
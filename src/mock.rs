
use std::f64::consts::PI;

use num_complex::Complex;

fn _sine(frequency: f64, phase: f64, t: f64) -> f64 {
    (2.0 * PI * frequency * t + phase).sin()
}

pub fn mock_sine(frequency_list: Vec<f64>, phase_list: Vec<f64>, duration: usize, sample_rate: f64) -> Vec<f64> {
    let sample_count = sample_rate as usize * duration; // 样本数
    let mut wave = Vec::with_capacity(sample_count); // 创建一个向量来存储正弦信号

    for n in 0..sample_count {
        let t = n as f64 / sample_rate; // 计算当前样本的时间
        let mut value = 0.0;
        for (i, frequency) in frequency_list.iter().enumerate() {
            let phase = phase_list[i];
            value += _sine(*frequency, phase, t);
        }
        wave.push(value); // 将计算的正弦值推入向量

    }

    wave
}

fn _cosine(frequency: f64, phase: f64, t: f64) -> f64 {
    (2.0 * PI * frequency * t + phase).cos()
}

pub fn mock_cosine(frequency_list: Vec<f64>, phase_list: Vec<f64>, duration: usize, sample_rate: f64) -> Vec<f64> {
    let sample_count = sample_rate as usize * duration; // 样本数
    let mut wave = Vec::with_capacity(sample_count); // 创建一个向量来存储正弦信号

    for n in 0..sample_count {
        let t = n as f64 / sample_rate; // 计算当前样本的时间
        let mut value = 0.0;
        for (i, frequency) in frequency_list.iter().enumerate() {
            let phase = phase_list[i];
            value += _cosine(*frequency, phase, t);
        }
        wave.push(value); // 将计算的正弦值推入向量

    }

    wave
}

pub fn find_frequency_in_spectrum(spectrum: Vec<(f64, Complex<f64>)>, threshold: Option<f64>) -> Vec<(f64, Complex<f64>)> {
    let threshold = threshold.unwrap_or(100.0);
    spectrum.into_iter()
        .filter(|v| v.1.norm() >= threshold)
        .collect()
}

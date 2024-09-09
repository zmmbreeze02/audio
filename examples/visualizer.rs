use spectrum_analyzer::{samples_fft_to_spectrum, FrequencyLimit};
use spectrum_analyzer::windows::hann_window;
use spectrum_analyzer::scaling::divide_by_N_sqrt;

/// Minimal example.
fn main() {
    // Produce a sinusoid of maximum amplitude.
    // 正弦波的公式是：y(t) = A * sin(2 * π * f * t + ϕ)
    // - y(t) 是时间 t 时的波形值
    // - A 是振幅，表示波的最大值，这里是 1
    // - f 是频率，表示单位时间内波形重复的次数，这里是 440hz，在西方音乐上440hz的声音音调为标准音高
    // - 2πf 是角频率，表示单位时间内波形的角位移
    // - ϕ 是相位角，表示波形相对于时间轴的偏移，这里是 0
    // - t 是时间，在这里是 sample_clock / sample_rate
    let sample_rate = 48000f32;
    let sample_length = 2i32.pow(10);
    let mut samples: Vec<f32> = Vec::new();
    for n in 0..sample_length {
        let sample_clock = (n + 1) as f32 % sample_rate;
        samples.push((sample_clock * 200.0 * 2.0 * std::f32::consts::PI / sample_rate).sin());
    }
    // println!("{:?}, len={}", samples, samples.len());


    // YOU need to implement the samples source; get microphone input for example
    // let samples: &[f32] = &[0.0, 3.14, 2.718, -1.0, -2.0, -4.0, 7.0, 6.0];
    // apply hann window for smoothing; length must be a power of 2 for the FFT
    // 2048 is a good starting point with 44100 kHz
    let hann_window = hann_window(&samples);
    // calc spectrum
    let spectrum_hann_window = samples_fft_to_spectrum(
        // (windowed) samples
        // &samples,
        &samples,
        // sampling rate
        sample_rate as u32,
        // optional frequency limit: e.g. only interested in frequencies 50 <= f <= 150?
        FrequencyLimit::Range(0.0, 1000.0),
        // optional scale
        Some(&divide_by_N_sqrt),
    );
    match spectrum_hann_window {
        Ok(t) => {
            for (fr, fr_val) in t.data().iter() {
                println!("{}Hz => {}", fr, fr_val)
            }
        },
        Err(e) => {
            println!("{:?}", e);
        },
    }

}
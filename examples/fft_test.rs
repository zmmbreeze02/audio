
use std::ops::Range;
use csv::Writer;
use charming::{
    component::Legend,
    component::Axis,
    element::AxisType,
    series::Line,
    Chart, HtmlRenderer
};
use audio::fft::calc_spectrum_by_fft;
use audio::mock::mock_sine;
use num_complex::Complex;
use audio::window::hanning;

fn write_csv(file_name: &str, frequency_axis: &[String], data: &[f64]) -> anyhow::Result<()> {
    let mut writer = Writer::from_path(file_name)?;
    writer.write_record(&["frequency", "db"])?;
    for (i, data) in data.iter().enumerate() {
        writer.write_record([&frequency_axis[i], &data.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

fn draw_line_chart(file_name: &str, axis: Vec<String>, data: Vec<f64>) -> anyhow::Result<()> {
    let line_chart = Chart::new()
    .legend(Legend::new().top("bottom"))
    .x_axis(
        Axis::new()
            .type_(AxisType::Category)
            .data(axis),
    )
    .y_axis(Axis::new().type_(AxisType::Value))
    .series(Line::new().data(data));
    let mut renderer = HtmlRenderer::new("chart", 1920, 750);
    renderer.save(&line_chart, file_name)?;

    anyhow::Ok(())
}

fn draw_spectrum_chart(file_name: &str, data: Vec<(f64, Complex<f64>)>) -> anyhow::Result<()> {
    let sample_count = data.len() / 2;
    let data = data.split_at(sample_count).0.to_vec();
    let data_with_magnitude: Vec<f64> = data.iter().map(|r| r.1.norm()).collect();
    let max_magnitude = data_with_magnitude.iter().fold(0.0, |max, &value| {
        if value > max {
            value
        } else {
            max
        }
    });
    let data_with_log: Vec<f64> = data_with_magnitude.iter().map(|mag|
        20.0 * (mag / max_magnitude).log10()
    ).collect();
    let frequency_axis: Vec<String> = data.iter().map(|(k,_)| {
        format!("{}hz", k)
    }).collect();
    write_csv(&file_name.replace(".html", ".csv"), &frequency_axis, &data_with_log)?;
    draw_line_chart(file_name, frequency_axis, data_with_log)
}

fn draw_time_chart(file_name: &str, sample_rate: f64, data: Vec<f64>) -> anyhow::Result<()> {
    let sample_count = data.len();
    let time_axis: Vec<String> = (0..=sample_count).map(|i| {
        format!("{:.1}s", i as f64 / sample_rate)
    }).collect();
    draw_line_chart(file_name, time_axis, data)
}

fn test_sin_wave(frequency_list: Vec<f64>, sample_rate: f64, duration: usize, silence_rage: Option<Range<usize>>) -> anyhow::Result<()> {
    let mut wave = mock_sine(frequency_list.clone(), vec![0.0], duration, sample_rate);
    if let Some(range) = silence_rage.clone() {
        for i in range {
            wave[i] = 0.0;
        }
    }
    let hanning_window = hanning(wave.len());
    for i in 0..wave.len() {
        wave[i] *= hanning_window[i];
    }

    let spectrum = calc_spectrum_by_fft(&wave, sample_rate)?;
    // let r = find_frequency_in_spectrum(spectrum, None);
    let name = format!(
        "{}hz_{}s_{}hz{}.html", sample_rate, duration,
        frequency_list.iter().map(|n| n.to_string()).collect::<Vec<String>>().join("-"),
        silence_rage.map(|r| format!("_({}-{})", r.start, r.end)).unwrap_or("".to_string())
    );
    draw_time_chart(&format!("sine_wave_{}", name), sample_rate, wave)?;
    draw_spectrum_chart(&format!("sine_spectrum_{}", name), spectrum)?;
    anyhow::Ok(())
}

fn main() -> anyhow::Result<()> {
    test_sin_wave(vec![100.0], 1024.0, 1, None)?;
    test_sin_wave(vec![100.0], 1024.0, 1, Some(100..200))?;
    anyhow::Ok(())
}

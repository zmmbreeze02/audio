
use std::f64::consts::PI;
use num_complex::{Complex, Complex64};
use charming::{
    component::Legend,
    component::Axis,
    element::AxisType,
    series::Line,
    Chart, HtmlRenderer
};
use audio::dft::dft;

fn sine(frequency: f64, phase: f64, t: f64) -> f64 {
    (2.0 * PI * frequency * t + phase).sin()
}

fn dot_product(v1: (f64, f64), v2: (f64, f64)) -> f64 {
    v1.0 * v2.0 + v1.1 * v2.1
}

fn magnitude(v: (f64, f64)) -> f64 {
    dot_product(v, v).sqrt()
}

fn phase(v: (f64, f64)) -> f64 {
    v.1.atan2(v.0)
}

fn draw_line_chart(file_name: String, axis: Vec<String>, data: Vec<f64>) -> anyhow::Result<()> {
    let line_chart = Chart::new()
    .legend(Legend::new().top("bottom"))
    .x_axis(
        Axis::new()
            .type_(AxisType::Category)
            .data(axis),
    )
    .y_axis(Axis::new().type_(AxisType::Value))
    .series(Line::new().data(data));
    let mut renderer = HtmlRenderer::new("chart", 1000, 800);
    renderer.save(&line_chart, file_name)?;

    anyhow::Ok(())
}


const SAMPLE_RATE: f64 = 1000.0; // 采样率
const DURATION: f64 = 10.0; // 信号时长，单位秒
const FREQUENCY: f64 = 5.0; // 信号频率，单位Hz

fn main() -> anyhow::Result<()> {
    let sample_count = (SAMPLE_RATE * DURATION) as usize; // 样本数
    let mut sine_wave = Vec::with_capacity(sample_count); // 创建一个向量来存储正弦信号

    for n in 0..sample_count {
        let t = n as f64 / SAMPLE_RATE; // 计算当前样本的时间
        let value = (2.0 * std::f64::consts::PI * FREQUENCY * t).cos();
        sine_wave.push(value); // 将计算的正弦值推入向量
    }
    let sine_r_map: Vec<(f64, f64)> = dft(&sine_wave)?.into_iter().map(|r| (r.re, r.im)).collect();
    let sine_magnitude_map: Vec<f64> = sine_r_map.to_vec().into_iter().map(|r| magnitude(r)).collect();
    let sine_phase_map: Vec<f64> = sine_r_map.to_vec().into_iter().map(|r: (f64, f64)| phase(r)).collect();

    println!("");
    let (max_index, max_value) = sine_magnitude_map[1..(sample_count/2)]
        .iter()
        .enumerate()
        .fold((0, sine_magnitude_map[0]), |(index, max_value), (new_index, new_value)| {
            let new_value = *new_value;
            if new_value > max_value {
                (new_index, new_value)
            } else {
                (index, max_value)
            }
        });
    println!("Max Sine Frequency {}hz, Value {}", (max_index + 1) as f64 * SAMPLE_RATE / sample_count as f64, max_value);


    let time_axis: Vec<String> = (0..=sample_count).map(|i| format!("{:.1}s", i as f64/SAMPLE_RATE as f64)).collect();
    draw_line_chart("./sine_time.html".to_string(), time_axis.to_vec(), sine_wave)?;
    let frequency_axis: Vec<String> = (0..=sample_count).map(|i| format!("{}hz", (i + 1) as f64 * SAMPLE_RATE / sample_count as f64)).collect();
    draw_line_chart("./sine_frequency.html".to_string(), frequency_axis.to_vec(), sine_magnitude_map)?;
    draw_line_chart("./sine_phase.html".to_string(), frequency_axis.to_vec(), sine_phase_map)?;
    draw_line_chart("./sine_real.html".to_string(), frequency_axis.to_vec(), sine_r_map.to_vec().into_iter().map(|(r, _)| r).collect())?;
    draw_line_chart("./sine_imag.html".to_string(), frequency_axis.to_vec(), sine_r_map.into_iter().map(|(_, i)| i).collect())?;

    anyhow::Ok(())
}


use std::f64::consts::PI;
use rand::Rng;

fn euler(frequency: f64, t: f64) -> (f64, f64) {
    let theta = -2.0 * PI * frequency * t;
    (theta.cos(), theta.sin())
}

fn cosine(frequency: f64, phase: f64, t: f64) -> f64 {
    (2.0 * PI * frequency * t + phase).cos()
}

fn sine(frequency: f64, phase: f64, t: f64) -> f64 {
    (2.0 * PI * frequency * t + phase).sin()
}

fn tangent(frequency: f64, phase: f64, t: f64) -> f64 {
    (PI * frequency * t + phase).tan()
}

fn product(v1: (f64, f64), v2: (f64, f64)) -> (f64, f64) {
    (v1.0 * v2.0 - v1.1 * v2.1, v1.0 * v2.1 + v1.1 * v2.0)
}

fn dot_product(v1: (f64, f64), v2: (f64, f64)) -> f64 {
    v1.0 * v2.0 + v1.1 * v2.1
}

fn magnitude(v: (f64, f64)) -> f64 {
    dot_product(v, v).sqrt()
}

fn add(v1: (f64, f64), v2: (f64, f64)) -> (f64, f64) {
    (v1.0 + v2.0, v1.1 + v2.1)
}

fn main() -> anyhow::Result<()> {
    let sample_count = 8000; // 10s
    let src_sine_frequency = 100.0;
    let src_sine_phase = rand::thread_rng().gen_range(1..=10) as f64;
    let src_cosine_frequency = 1500.0;
    let src_cosine_phase = rand::thread_rng().gen_range(1..=10) as f64;
    // let src_tan_frequency = 1.0;
    // let src_tan_phase = rand::thread_rng().gen_range(1..=10) as f64;

    let mut sine_r_arr = Vec::<(f64, f64)>::new(); 
    let mut cosine_r_arr = Vec::<(f64, f64)>::new(); 
    let mut tan_r_arr = Vec::<(f64, f64)>::new();
    for f in 0..=sample_count {
        let mut sin_r = (0.0, 0.0);
        let mut cos_r = (0.0, 0.0);
        let mut tan_r = (0.0, 0.0);
        for n in 0..=sample_count {
            let t = n as f64 / sample_count as f64;
            let euler_vector = euler(f as f64, t);
            let sin_value = sine(src_sine_frequency, 0.0, t);
            // println!("E = ({}, {}), S = {}", euler_vector.0, euler_vector.1, sin_value);
            // let rr = product((sin_value, 0.0), euler_vector);
            let rr = (sin_value * euler_vector.0, sin_value * euler_vector.1);
            // println!("R = {}, {}", rr.0, rr.1);
            sin_r = add(sin_r, rr);

            // let cos_value = cosine(src_cosine_frequency, 0.0, t);
            // let rr = product((cos_value, 0.0), euler_vector);
            // cos_r = add(cos_r, rr);

        //     let tan_value = tangent(src_tan_frequency, src_tan_phase, t);
        //     let rr = product((tan_value, 0.0), euler_vector);
        //     tan_r = add(tan_r, rr);
        }

        sine_r_arr.push(sin_r);
        cosine_r_arr.push(cos_r);
        tan_r_arr.push(tan_r);
        if f == 1900 || f == 100 {
            println!("F={}: Sine ({},{}), {}", f, sin_r.0, sin_r.1, magnitude(sin_r));
        }
    }

    println!("");
    let (max_index, max_value) = sine_r_arr
        .iter()
        .enumerate()
        .fold((0, magnitude(sine_r_arr[0])), |(index, max_value), (new_index, new_value)| {
            let new_value = magnitude(*new_value);
            if new_value > max_value {
                (new_index, new_value)
            } else {
                (index, max_value)
            }
        });
    println!("Max Sine Frequency {}hz, Value {}", max_index, max_value);

    // let (max_index, max_value) = cosine_r_arr
    //     .iter()
    //     .enumerate()
    //     .fold((0, magnitude(cosine_r_arr[0])), |(index, max_value), (new_index, new_value)| {
    //         let new_value = magnitude(*new_value);
    //         if new_value > max_value {
    //             (new_index, new_value)
    //         } else {
    //             (index, max_value)
    //         }
    //     });
    // println!("Max Cosine Frequency {}hz, Value {}", max_index, max_value);

    // let (max_index, max_value) = tan_r_arr
    //     .iter()
    //     .enumerate()
    //     .fold((0, magnitude(tan_r_arr[0])), |(index, max_value), (new_index, new_value)| {
    //         let new_value = magnitude(*new_value);
    //         if new_value > max_value {
    //             (new_index, new_value)
    //         } else {
    //             (index, max_value)
    //         }
    //     });
    // println!("Max Tangent Frequency {}hz, Value {}", max_index, max_value);

    anyhow::Ok(())
}

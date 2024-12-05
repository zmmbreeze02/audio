use num_complex::Complex;
use thiserror::Error;
use std::f64::consts::PI;
use std::collections::HashMap;



#[derive(Error, Debug)]
pub enum FFTError {
    #[error("Size of samples is not enough {0}")]
    NotEnoughSamples(usize),
    #[error("Size of samples is not power of two {0}")]
    NotPowerOfTwo(usize),
    #[error("unknown data store error")]
    Unknown,
}

pub fn fft(samples: &[f64]) -> Result<Vec<Complex<f64>>, FFTError> {
    let len = samples.len();
    if len <= 1 {
        return Err(FFTError::NotEnoughSamples(len));
    }
    if !len.is_power_of_two() {
        return Err(FFTError::NotPowerOfTwo(len));
    }

    let mut result = vec![Complex::ZERO; len];
    let mut cache = HashMap::new();
    for k in 0..len {
        result[k] = _fft_recursion_k(len, k, samples.to_vec(), &mut cache);
    }

    Ok(result)
}

fn _bit_reverse() {

}





// result = W_N^k = e^{-j * 2\pi * k / N}
fn _gen_w_n(k: usize, len: usize) -> Complex<f64> {
    if k == 0 || k == len {
        return Complex::new(1.0, 0.0);
    }
    if len == 2 * k {
        return Complex::new(-1.0, 0.0);
    }
    if len == 4 * k {
        return Complex::new(0.0, -1.0);
    }
    if len * 3 == 4 * k {
        return Complex::new(0.0, 1.0);
    }
    let len_f64 = len as f64;
    let k_f64 = k as f64;
    let theta = -2.0 * PI * k_f64 / len_f64;
    Complex::new(theta.cos(), theta.sin())
}

fn _fft_recursion_k(position: usize, k: usize, samples: Vec<f64>, cache: &mut HashMap<(usize, usize, usize), (Complex<f64>, Complex<f64>)>) -> Complex<f64> {
    let len = samples.len();
    if len == 1 {
        // X[k] = x[0] * W_N^{0*k} = x[0]
        return Complex::new(samples[0], 0.0);
    }

    if k >= len/2 {
        // X[k] = X_{even}[k] - W_N^k * X_{odd}[k]
        // use cached X_{odd}[k] and W_N^k * X_{even}[k]
        if let Some((l, r)) = cache.get(&(len, k - len/2, position)) {
            return l - r;
        }
        unreachable!();
    }

    // println!("k:{}, len:{}, pos:{}, samples:{:?}", k, len, position, samples);

    // X[k] = X_{even}[k] + W_N^k * X_{odd}[k]
    // split samples into odd and even part
    let (even_indices, odd_indices): (Vec<_>, Vec<_>) = samples
        .iter()
        .enumerate()
        .partition(|&(index, _)| index % 2 == 0);
    let odd_samples: Vec<f64> = odd_indices.into_iter().map(|(_, value)| *value).collect();
    let even_samples: Vec<f64> = even_indices.into_iter().map(|(_, value)| *value).collect();
    // combine odd and even part
    let l = _fft_recursion_k(position, k, even_samples, cache);
    let r = _gen_w_n(k, len) * _fft_recursion_k(position | (len/2), k, odd_samples, cache);

    cache.insert((len, k, position), (l, r));
    l + r
}

// calc spectrum as `Vec<(frequency, complex_result)>`
pub fn calc_spectrum_by_fft(samples: &[f64], sample_rate: f64) -> Result<Vec<(f64, Complex<f64>)>, FFTError> {
    let sample_count = samples.len() as f64;
    let spectrum = fft(samples)?
        .into_iter()
        .enumerate()
        .map(|(k, c)| (k as f64 * sample_rate / sample_count, c))
        .collect();
    Ok(spectrum)
}


#[cfg(test)]
mod tests {
    use super::{calc_spectrum_by_fft, FFTError};
    use super::super::mock::{mock_sine, mock_cosine, find_frequency_in_spectrum};

    #[test]
    fn test_calc_spectrum_by_fft() -> Result<(), FFTError> {
        // let sine = mock_sine(vec![1.0], vec![0.0], 1, 8.0);
        // println!("original: {:?}", sine);
        // let spectrum = calc_spectrum_by_fft(&sine, 8.0)?;
        // let sp: Vec<(f64,f64)> = spectrum.clone().into_iter().map(|(k, v)| (k, v.norm())).collect();
        // println!("result: {:?}", sp);
        // let r = find_frequency_in_spectrum(spectrum, None);
        // // println!("\n{:?}", r);
        // assert_eq!(r.len(), 1);
        // assert_eq!(r[0].0, 1.0);

        let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0], vec![0.0], 2, 1024.0), 1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        // println!("\n{:?}", r);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].0, 5.0);

        let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0, 10.0], vec![0.0, 10.0], 2, 1024.0), 1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        // println!("{:?}", r);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);

        let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0, 10.0, 7000.0], vec![0.0, 10.0, 10000.0], 2, 16.0*1024.0), 16.0*1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        // println!("{:?}", r);
        assert_eq!(r.len(), 3);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);
        assert_eq!(r[2].0, 7000.0);

        let spectrum = calc_spectrum_by_fft(&mock_cosine(vec![5.0, 10.0, 7000.0], vec![0.0, 10.0, 10000.0], 2, 16.0*1024.0), 16.0*1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        // println!("{:?}", r);
        assert_eq!(r.len(), 3);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);
        assert_eq!(r[2].0, 7000.0);

        Ok(())
    }

}


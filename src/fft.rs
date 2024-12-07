use num_complex::Complex;
use rand::seq::index;
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

/**
 * Use Gold Rader bit reversal algorithm to swap the elements of samples.
 * 
 * This algorithm can do `input += 1` in the order of bit_reverse.
 * 1. find the bit position of the first 0 as `n`.
 * 2. Turn all the higher bits after the `n`th bit to 0.
 * 3. Set the `n`th bit into 1.
 * 
 * For example: input is 6 = 0110
 * 1. n = 1, bit = 1000
 * 2. output = 1110
 * 3. output = 1110
 * 
 * For example: input is 13 = 1101
 * 1. n = 3, bit = 0010
 * 2. output = 0001
 * 3. output = 0011
 * 
 * For example: input is 13 = 1000
 * 1. n = 2, bit = 0100
 * 2. output = 0000
 * 3. output = 0100
 */
fn _bit_reverse(samples: &mut [f64]) {
    let len = samples.len();
    let max_index = len - 1;
    let mut reserved_index = 0;

    for index in 1..len {
        // use the previous `reserved_index`, and "grow it up" as the order of bit_reverse
        let previous_reserved_index = reserved_index;
    
        // `bit` must be 2^n
        // `n` represent of the bit position of the first 0
        // (of course in the backward order)
        let mut bit = len;

        loop {
            // move backward
            bit >>= 1;

            // find the bit position of the first 0, 
            // `bit` must <= `max_index - previous_reserved_index`, for example:
            // 1111 - 1101 = 0010, bit should be 0010, and 0010 = 0010
            // 1111 - 1010 = 0101, bit should be 0100, and 0100 < 0101
            // 1111 - 0000 = 1111, bit should be 1000, and 1000 < 1111
            if bit <= max_index - previous_reserved_index {
                break;
            }
        }

        // mask of the `bit`, for example:
        // 1000 => 0111
        // 0100 => 0011
        let mask = bit - 1;
        // Set the `n`th bit of the `previous_reserved_index` to 1
        // and turn all the higher bits after the `n`th bit to 0
        reserved_index = (previous_reserved_index & mask) + bit;

        if reserved_index > index {
            samples.swap(index, reserved_index);
        }
    }
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
    use super::{calc_spectrum_by_fft, FFTError, _bit_reverse};
    use super::super::mock::{mock_sine, mock_cosine, find_frequency_in_spectrum};


    #[test]
    fn test_bit_reverse() {
        let mut samples: Vec<f64> = (0..4).map(|v| v as f64).collect();
        _bit_reverse(&mut samples);
        // samples.iter().enumerate().for_each(|(i, v)| println!("{} = {} = {:02b}", i, v.clone() as usize,v.clone() as usize));
        assert_eq!(samples, vec![0.0, 2.0, 1.0, 3.0]);

        let mut samples: Vec<f64> = (0..8).map(|v| v as f64).collect();
        _bit_reverse(&mut samples);
        // samples.iter().enumerate().for_each(|(i, v)| println!("{} = {} = {:03b}", i, v.clone() as usize,v.clone() as usize));
        assert_eq!(samples, vec![0.0, 4.0, 2.0, 6.0, 1.0, 5.0, 3.0, 7.0]);

        let mut samples: Vec<f64> = (0..16).map(|v| v as f64).collect();
        _bit_reverse(&mut samples);
        // samples.iter().enumerate().for_each(|(i, v)| println!("{} = {} = {:04b}", i, v.clone() as usize,v.clone() as usize));
        assert_eq!(samples, vec![0.0, 8.0, 4.0, 12.0, 2.0, 10.0, 6.0, 14.0, 1.0, 9.0, 5.0, 13.0, 3.0, 11.0, 7.0, 15.0]);
    }

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

        // let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0], vec![0.0], 2, 1024.0), 1024.0)?;
        // let r = find_frequency_in_spectrum(spectrum, None);
        // // println!("\n{:?}", r);
        // assert_eq!(r.len(), 1);
        // assert_eq!(r[0].0, 5.0);

        // let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0, 10.0], vec![0.0, 10.0], 2, 1024.0), 1024.0)?;
        // let r = find_frequency_in_spectrum(spectrum, None);
        // // println!("{:?}", r);
        // assert_eq!(r.len(), 2);
        // assert_eq!(r[0].0, 5.0);
        // assert_eq!(r[1].0, 10.0);

        // let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0, 10.0, 7000.0], vec![0.0, 10.0, 10000.0], 2, 16.0*1024.0), 16.0*1024.0)?;
        // let r = find_frequency_in_spectrum(spectrum, None);
        // // println!("{:?}", r);
        // assert_eq!(r.len(), 3);
        // assert_eq!(r[0].0, 5.0);
        // assert_eq!(r[1].0, 10.0);
        // assert_eq!(r[2].0, 7000.0);

        // let spectrum = calc_spectrum_by_fft(&mock_cosine(vec![5.0, 10.0, 7000.0], vec![0.0, 10.0, 10000.0], 2, 16.0*1024.0), 16.0*1024.0)?;
        // let r = find_frequency_in_spectrum(spectrum, None);
        // // println!("{:?}", r);
        // assert_eq!(r.len(), 3);
        // assert_eq!(r[0].0, 5.0);
        // assert_eq!(r[1].0, 10.0);
        // assert_eq!(r[2].0, 7000.0);

        Ok(())
    }

}


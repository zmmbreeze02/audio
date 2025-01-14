use num_complex::Complex;
use thiserror::Error;
use std::f64::consts::PI;



#[derive(Error, Debug)]
pub enum FFTError {
    #[error("Size of input is not enough {0}")]
    NotEnoughSamples(usize),
    #[error("Size of input is not power of two {0}")]
    NotPowerOfTwo(usize),
    #[error("unknown data store error")]
    Unknown,
}

pub fn fft(input: &[f64]) -> Result<Vec<Complex<f64>>, FFTError> {
    let len = input.len();
    if len <= 1 {
        return Err(FFTError::NotEnoughSamples(len));
    }
    if !len.is_power_of_two() {
        return Err(FFTError::NotPowerOfTwo(len));
    }

    // Turn into Complex Vector
    let mut output: Vec<Complex<f64>> = input.into_iter().map(|s| Complex::new(*s, 0.0)).collect();

    _bit_reverse(&mut output);
    _butterflies(&mut output);

    Ok(output)
}

pub fn ifft(mut input: Vec<Complex<f64>>) -> Result<Vec<Complex<f64>>, FFTError> {
    let len = input.len();
    if len <= 1 {
        return Err(FFTError::NotEnoughSamples(len));
    }
    if !len.is_power_of_two() {
        return Err(FFTError::NotPowerOfTwo(len));
    }

    // conjugate 
    input.iter_mut().for_each(|c| c.im = -c.im );

    // fft
    _bit_reverse(&mut input);
    _butterflies(&mut input);

    // conjugate and divide by length 
    input.iter_mut().for_each(|c| {
        c.re = c.re / len as f64;
        c.im = -c.im / len as f64;
    });

    Ok(input)
}

/**
 * Use Gold Rader bit reversal algorithm to swap the elements of input.
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
fn _bit_reverse(input: &mut Vec<Complex<f64>>) {
    let len = input.len();
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
            input.swap(index, reserved_index);
        }
    }
}

/**
 * Use butterflies calculation to calc the FFT result.
 */
fn _butterflies(input: &mut Vec<Complex<f64>>) {
    let n = input.len();
    let max_stage = n.ilog2();

    for stage in 1..=max_stage {
        // Length of each block
        let l = 2usize.pow(stage);
        // Split input into blocks
        // for example when n=8, the blocks are:
        // stage=1,l=2 => 0,1|2,3|4,5|6,7
        // stage=2,l=4 => 0,1,2,3|4,5,6,7
        // stage=3,l=8 => 0,1,2,3,4,5,6,7
        for block_head_index in (0..n).step_by(l) {
            // half of L
            let half_l = l/2;
            // iterate through each parts
            let block_tail_index = block_head_index + half_l;
            for i in block_head_index..block_tail_index {
                let twiddle = _calc_twiddle(i % half_l, l);
                let tmp = input[i];
                input[i] = input[i] + input[i + half_l] * twiddle;
                input[i + half_l] = tmp - input[i + half_l] * twiddle;
            }
        }
    }
}


// result = W_N^k = e^{-j * 2\pi * k / N}
//        = cos(-2\pi * k / N) + i * sin(-2\pi * k / N)
fn _calc_twiddle(k: usize, len: usize) -> Complex<f64> {
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


// calc spectrum as `Vec<(frequency, complex_result)>`
pub fn calc_spectrum_by_fft(input: &[f64], sample_rate: f64) -> Result<Vec<(f64, Complex<f64>)>, FFTError> {
    let sample_count = input.len() as f64;
    let spectrum = fft(input)?
        .into_iter()
        .enumerate()
        .map(|(k, c)| (k as f64 * sample_rate / sample_count, c))
        .collect();
    Ok(spectrum)
}


#[cfg(test)]
mod tests {
    use num_complex::Complex;
    use super::{calc_spectrum_by_fft, fft, ifft, FFTError, _bit_reverse};
    use super::super::mock::{mock_sine, mock_cosine, find_frequency_in_spectrum};


    #[test]
    fn test_bit_reverse() {
        let mut input:Vec<Complex<f64>> = (0..4).map(|v| Complex::new(v as f64, 0.0)).collect();
        _bit_reverse(&mut input);
        let result:Vec<f64> = input.into_iter().map(|v| v.re).collect();
        assert_eq!(result, vec![0.0, 2.0, 1.0, 3.0]);

        let mut input:Vec<Complex<f64>> = (0..8).map(|v| Complex::new(v as f64, 0.0)).collect();
        _bit_reverse(&mut input);
        let result:Vec<f64> = input.into_iter().map(|v| v.re).collect();
        assert_eq!(result, vec![0.0, 4.0, 2.0, 6.0, 1.0, 5.0, 3.0, 7.0]);

        let mut input:Vec<Complex<f64>> = (0..16).map(|v| Complex::new(v as f64, 0.0)).collect();
        _bit_reverse(&mut input);
        let result:Vec<f64> = input.into_iter().map(|v| v.re).collect();
        assert_eq!(result, vec![0.0, 8.0, 4.0, 12.0, 2.0, 10.0, 6.0, 14.0, 1.0, 9.0, 5.0, 13.0, 3.0, 11.0, 7.0, 15.0]);
    }

    #[test]
    fn test_calc_spectrum_by_fft() -> Result<(), FFTError> {
        let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0], vec![0.0], 2, 1024.0), 1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].0, 5.0);

        let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0, 10.0], vec![0.0, 10.0], 2, 1024.0), 1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);

        let spectrum = calc_spectrum_by_fft(&mock_sine(vec![5.0, 10.0, 7000.0], vec![0.0, 10.0, 10000.0], 2, 16.0*1024.0), 16.0*1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        assert_eq!(r.len(), 3);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);
        assert_eq!(r[2].0, 7000.0);

        let spectrum = calc_spectrum_by_fft(&mock_cosine(vec![5.0, 10.0, 7000.0], vec![0.0, 10.0, 10000.0], 2, 16.0*1024.0), 16.0*1024.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        assert_eq!(r.len(), 3);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);
        assert_eq!(r[2].0, 7000.0);

        Ok(())
    }

    #[test]
    fn test_ifft() -> Result<(), FFTError> {
        let sample = mock_sine(vec![5.0], vec![0.0], 2, 1024.0);
        // println!("sample: {:?}\n", sample);
        let r = fft(&sample)?;
        let r = ifft(r)?;
        assert_eq!(r.len(), 2048);

        // imaginary should (almost) equal to 0
        let wrong_imag = r.clone().into_iter().map(|v| v.im).filter(|v| v.abs() > 1e-15).collect::<Vec<f64>>();
        // println!("imagination: {:?}\n", wrong_imag);
        assert_eq!(wrong_imag.len(), 0);

        // real should be (almost) equal to sample
        let r = r.into_iter().map(|v| v.re).collect::<Vec<f64>>();
        // println!("real: {:?}\n", r);
        let diff = r.iter().zip(sample.iter()).map(|(a, b)| a - b).filter(|v| v.abs() > 1e-15).collect::<Vec<f64>>();
        // println!("diff: {:?}\n", diff);
        assert_eq!(diff.len(), 0);

        Ok(())
    }


}


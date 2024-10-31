use num_complex::Complex;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FFTError {
    #[error("Size of samples is not power of two {0}")]
    NotPowerOfTwo(usize),
    #[error("unknown data store error")]
    Unknown,
}

pub fn fft(input: &[Complex<f64>]) -> Result<Vec<Complex<f64>>, FFTError> {
    let len = input.len();
    if !len.is_power_of_two() {
        return Err(FFTError::NotPowerOfTwo(len));
    }

    let input_even: &
}




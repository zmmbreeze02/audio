use num_complex::Complex;
use thiserror::Error;
use std::f64::consts::PI;

#[derive(Error, Debug)]
pub enum DFTError {
    // #[error("Size of samples is not power of two {0}")]
    // NotPowerOfTwo(usize),
    #[error("unknown data store error")]
    Unknown,
}

// X[k] = \sum_{n=0}^{N-1} x[n]*e^{-i2\pi\frac{k*n}{N}}
// \theta = -2\pi\frac{k*n}{N} 
// e^{\theta} = cos(\theta) + i*sin(\theta)
pub fn dft(samples: &[f64]) -> Result<Vec<Complex<f64>>, DFTError> {
    let len = samples.len();
    let len_f64 = len as f64;
    let mut result = vec![Complex::ZERO; len];

    for k in 0..len {
        let k_f64 = k as f64;
        let mut x_k = Complex::new(0.0, 0.0);

        for n in 0..len {
            let n_f64 = n as f64;
            let theta = -2.0 * PI * k_f64 * n_f64 / len_f64;
            x_k = x_k + Complex::new(samples[n], 0.0) * Complex::new(theta.cos(), theta.sin());
        }

        result[k] = x_k;
    }

    Ok(result)
}
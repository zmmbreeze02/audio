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

// calc spectrum as `Vec<(frequency, complex_result)>`
pub fn calc_spectrum_by_dft(samples: &[f64], sample_rate: f64) -> Result<Vec<(f64, Complex<f64>)>, DFTError> {
    let sample_count = samples.len() as f64;
    let spectrum = dft(samples)?
        .into_iter()
        .enumerate()
        .map(|(k, c)| (k as f64 * sample_rate / sample_count, c))
        .collect();
    Ok(spectrum)
}


#[cfg(test)]
mod tests {
    use super::{calc_spectrum_by_dft, dft, DFTError};
    use super::super::mock::{mock_sine, mock_cosine, find_frequency_in_spectrum};

    #[test]
    fn test_calc_spectrum_by_dft() -> Result<(), DFTError> {
        let spectrum = calc_spectrum_by_dft(&mock_sine(vec![5.0], vec![0.0], 2, 1000.0), 1000.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        // println!("{:?}", r);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].0, 5.0);

        let spectrum = calc_spectrum_by_dft(&mock_sine(vec![5.0, 10.0], vec![0.0, 10.0], 2, 1000.0), 1000.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        // println!("{:?}", r);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);

        let spectrum = calc_spectrum_by_dft(&mock_sine(vec![5.0, 10.0, 7000.0], vec![0.0, 10.0, 10000.0], 2, 16000.0), 16000.0)?;
        let r = find_frequency_in_spectrum(spectrum, None);
        // println!("{:?}", r);
        assert_eq!(r.len(), 3);
        assert_eq!(r[0].0, 5.0);
        assert_eq!(r[1].0, 10.0);
        assert_eq!(r[2].0, 7000.0);

        Ok(())
    }

}

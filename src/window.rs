/// Reference: https://github.com/numpy/numpy/blob/main/numpy/lib/_function_base_impl.py#L3267
use std::f64::consts::PI;

/**
 * Generate symmetric indices.
 * Equivalent to `np.arange(1-size, size, 2)` in Python
 */
fn calc_with_symmetric_indices<F>(size: usize, callback: F) -> Vec<f64>
where
    F: Fn(isize) -> f64,
{
    if size <= 0 {
        return vec![];
    }
    if size == 1 {
        return vec![1.0];
    }

    let start = 1 - size as isize;
    let end = size as isize;
    let step = 2;

    (start..end)    // [start, end)
        .step_by(step)
        .map(callback)
        .collect()
}

/**
 * Return the Hanning window.
 * w(n) = 0.5 - 0.5\\cos\\left(\\frac{2\\pi{n}}{M-1}\\right)
 *        \\qquad 0 \\leq n \\leq M-1
 */
pub fn hanning(size: usize) -> Vec<f64> {
    calc_with_symmetric_indices(size, |n| {
        0.5 + 0.5 * (PI * n as f64 / (size - 1) as f64).cos()
    })
}

/**
 * Return the Hamming window.
 * w(n) = 0.54 - 0.46\cos\left(\frac{2\pi{n}}{M-1}\right)
               \qquad 0 \leq n \leq M-1
 */
pub fn hamming(size: usize) -> Vec<f64> {
    calc_with_symmetric_indices(size, |n| {
        0.54 + 0.46 * (PI * n as f64 / (size - 1) as f64).cos()
    })
}

/**
 * Return the Blackman window.
 * w(n) = 0.42 - 0.5 \cos(2\pi n/M) + 0.08 \cos(4\pi n/M)
 */
pub fn blackman(size: usize) -> Vec<f64> {
    calc_with_symmetric_indices(size, |n| {
        0.42 + 0.5 * (PI * n as f64 / (size - 1) as f64).cos() + 0.08 * (2.0 * PI * n as f64 / (size - 1) as f64).cos()
    })
}

/**
 * Return the Bartlett window.
 *  w(n) = \frac{2}{M-1} \left(
              \frac{M-1}{2} - \left|n - \frac{M-1}{2} \right|
              \right)
 */
pub fn bartlett(size: usize) -> Vec<f64> {
    calc_with_symmetric_indices(size, |n| {
        if n <= 0 {
            1.0 + n as f64 / (size - 1) as f64
        } else {
            1.0 - n as f64 / (size - 1) as f64
        }
    })
}




#[cfg(test)]
mod tests {
    use super::{hanning, hamming, blackman, bartlett};

    #[test]
    fn test() {
        assert_eq!(bartlett(12), vec![0.0, 0.18181818181818177, 0.36363636363636365, 0.5454545454545454, 0.7272727272727273, 0.9090909090909091, 0.9090909090909091, 0.7272727272727273, 0.5454545454545454, 0.36363636363636365, 0.18181818181818177, 0.0]);
        assert_eq!(hanning(12), vec![0.0, 0.07937323358440945, 0.29229249349905684, 0.5711574191366425, 0.8274303669726426, 0.9797464868072487, 0.9797464868072487, 0.8274303669726426, 0.5711574191366425, 0.29229249349905684, 0.07937323358440945, 0.0]);
        assert_eq!(hamming(12), vec![0.08000000000000002, 0.15302337489765672, 0.3489090940191323, 0.6054648256057111, 0.8412359376148312, 0.9813667678626689, 0.9813667678626689, 0.8412359376148312, 0.6054648256057111, 0.3489090940191323, 0.15302337489765672, 0.08000000000000002]);
        assert_eq!(blackman(12), vec![-1.3877787807814457e-17, 0.032606434624560324, 0.159903634783434, 0.4143979812474828, 0.7360451799107798, 0.9670467694337431, 0.9670467694337431, 0.7360451799107798, 0.4143979812474828, 0.159903634783434, 0.032606434624560324, -1.3877787807814457e-17]);
    }
}
use thiserror::Error;
use std::fmt::Display;

// biquad filter
// y[n] = b_0 * x[n] + b_1 * x[n-1] + b_2 * x[n-2] - a_1 * y[n-1] - a_2 * y[n-2]
//
//        b_0 + b_1 • z^(-1) + b_2 • z^(-2)
// H(z) = ---------------------------------
//         1 + a_1 • z^(-1) + a_2 • z^(-2)
//
#[derive(Error, Debug)]
pub struct BiquadFilter {
    b: (f64, f64, f64),
    a: (f64, f64),
}

impl BiquadFilter {
    pub fn process(&mut self, x: Vec<f64>) -> Vec<f64> {
        let mut y = Vec::new();
        let mut x1  = 0.0;
        let mut x2 = 0.0;
        let mut y1 = 0.0;
        let mut y2 = 0.0;
        for x0 in x {
            let y0 = self.b.0 * x0 + self.b.1 * x1 + self.b.2 * x2 - self.a.0 * y1 - self.a.1 * y2;
            x2 = x1;
            x1 = x0;
            y2 = y1;
            y1 = y0;
            y.push(y0);
        }
        y
    }
}

impl Display for BiquadFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BiquadFilter {{ b: {:?}, a: {:?} }}", self.b, self.a)
    }
}

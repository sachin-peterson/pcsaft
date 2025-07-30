use crate::structures::*;
use crate::constants::*;

// Hard-chain term
pub fn a_hc(model: &PCSAFTModel, V: f64, _T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> f64 {
    let (d, z0, z1, z2, z3, m_bar) = data;
    let m = &model.params.m;
    let sum_z: f64 = z.iter().sum();
    
    let c1 = 1.0 / (1.0-z3);
    let c2 = 3.0 * z2 / (1.0-z3).powi(2);
    let c3 = 2.0 * z2.powi(2) / (1.0-z3).powi(3);

    let a_hs = if z3 != &0.0 {
        bmcs_hs(z0, z1, z2, z3)
    } else {
        bmcs_hs_zero_v(model, V, z, d)
    };

    let mut res = 0.0;
    for i in 0..z.len() {
        let d_i = d[i];
        let g_hs_i = c1 + d_i*d_i/(d_i+d_i) * c2 + (d_i*d_i/(d_i+d_i)).powi(2) * c3;
        res += z[i] * (m[i]-1.0) * g_hs_i.ln();
    }

    m_bar*a_hs - res/sum_z
}

// Hard-sphere term (Boublik-Mansoori-Carnahan-Starling)
// Case I: Non-zero volume
pub fn bmcs_hs(z0: &f64, z1: &f64, z2: &f64, z3: &f64) -> f64 {
    let z3m1 = 1.0 - z3;
    let z3m1_sq = z3m1*z3m1;

    let term1 = 3.0*z1*z2 / z3m1;
    let term2 = z2.powi(3) / (z3*z3m1_sq);
    let term3 = (z2.powi(3)/z3.powi(2) - z0) * (1.0 - z3).ln();

    (term1 + term2 + term3) / z0
}

// Case II: Zero volume
pub fn bmcs_hs_zero_v(model: &PCSAFTModel, V: f64, z: &[f64], d: &[f64]) -> f64 {
    let mut z0v = 0.0;
    let mut z1v = 0.0;
    let mut z2v = 0.0;
    let mut z3v = 0.0;

    for i in 0..z.len() {
        let m_i = model.params.m[i];
        let d_i = d[i];
        let xS = z[i] * m_i;
        z0v += xS;
        z1v += xS * d_i;
        z2v += xS * d_i.powi(2);
        z3v += xS * d_i.powi(3);
    }
    
    let c = PI/6.0 * NA;
    let d0 = c * z0v;
    let d1 = c * z1v;
    let d2 = c * z2v;
    let d3 = c * z3v;
    let rho = 1.0/V;

    let z2 = rho * z2v;
    let z3 = rho * z3v;
    let dz3 = 1.0 - z3;
    let log_dz3 = dz3.ln();

    3.0 * d1/d0 * z2/dz3
        + d2.powi(2)*z2 / (d3*d0*dz3.powi(2))
        - log_dz3
        + d2.powi(3) / (d3.powi(2)*d0) * log_dz3
}

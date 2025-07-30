use crate::structures::*;
use crate::constants::*;

// Dispersion term
pub fn a_disp(model: &PCSAFTModel, V: f64, T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> f64 {
    let (_, _, _, _, _, m_bar) = data;
    let sum_z: f64 = z.iter().sum();
    let rho = PI*NA*sum_z / V;

    let (m2eo3_1, m2eo3_2) = m2eo3(model, V, T, z);

    -2.0*rho*i(model, V, T, z, data, 1)*m2eo3_1 - m_bar*rho*c1(model, V, T, z, data)*m2eo3_2*i(model, V, T, z, data, 2)
}

// Mixture parameters
pub fn m2eo3(model: &PCSAFTModel, _V: f64, T: f64, z: &[f64]) -> (f64, f64) {
    let m = &model.params.m;
    let sigma = &model.params.sigma;
    let epsilon = &model.params.epsilon;

    let mut v1 = 0.0;
    let mut v2 = 0.0;

    for i in 0..z.len() {
        for j in 0..z.len() {
            let cst = z[i]*z[j]*m[i]*m[j]*sigma[i][j].powi(3);
            let exp = epsilon[i][j] / T;
            v1 += cst * exp;
            v2 += cst * exp.powi(2);
        }
    }

    let sum_z: f64 = z.iter().sum();
    let k = 1.0 / sum_z.powi(2);
    (k*v1, k*v2) 
}

// Integral of the radial distribution function
pub fn i(_model: &PCSAFTModel, _V: f64, _T: f64, _z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64), n: usize) -> f64 {
    let (_, _, _, _, eta, m_bar) = data;
    
    // 1st or 2nd order correlation
    let corr = if n == 1 {
        &PCSAFT_CONSTS.aij
    } else if n == 2 {
        &PCSAFT_CONSTS.bij
    } else {
        panic!("Invalid n value")
    };

    let m1 = (m_bar - 1.0) / m_bar;
    let m2 = (m_bar - 1.0) / m_bar * (m_bar - 2.0) / m_bar;

    let mut res = 0.0;
    let mut eta_power = 1.0;

    for &(c1, c2, c3) in corr {
        let ki = c1 + m1*c2 + m2*c3;
        res += ki*eta_power;
        eta_power *= eta;
    }

    res
}

// Compressibility expression
pub fn c1(_model: &PCSAFTModel, _V: f64, _T: f64, _z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> f64 {
    let (_, _, _, _, eta, m_bar) = data;

    // Evaluate polynomial: 0 + 20η - 27η² + 12η³ - 2η⁴
    let poly = 20.0*eta - 27.0*eta.powi(2) + 12.0*eta.powi(3) - 2.0*eta.powi(4);

    let term1 = 1.0 + m_bar*(8.0*eta - 2.0*eta.powi(2))/(1.0-eta).powi(4);
    let term2 = (1.0-m_bar)*poly / ((1.0-eta)*(2.0-eta).powi(2));

    1.0 / (term1 + term2)
}

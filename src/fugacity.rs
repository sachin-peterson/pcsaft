use crate::structures::*;
use crate::constants::*;
use crate::newton::*;
use crate::model::*;

// STEP 1: Guess eta, calculate rho (number density)
// STEP 2: Calculate zeta values
// STEP 3: Estimate Z = 1 + Z_hc + Z_disp (+ Z_assoc)
// STEP 4: Compute P = z * rho * RT
// STEP 5: Compute residual Psys - Pcalc
// STEP 6: Root solve for correct eta
// STEP 7: With correct Z, compute fugacity coefficient

pub fn fugacity_coeff(k: usize, x: &[f64], Z: f64, T: f64, rho: f64, model: &PCSAFTModel) -> f64 {
    let V = 1.0 / rho;
    let n = x.len();
    let h = 1e-6;

    // ∂a_res/∂x_k
    let df_dxk = |xk_perturbed: f64| {
        let mut x_perturbed = x.to_vec();
        x_perturbed[k] = xk_perturbed;
        a_res(model, V, T, &x_perturbed)
    };
    let dadx_k = five_pt_stencil(&df_dxk, x[k], h);

    // ∑ x_j ∂a_res/∂x_j
    let mut sum = 0.0;
    for j in 0..n {
        let df_dxj = |xj_perturbed: f64| {
            let mut x_perturbed = x.to_vec();
            x_perturbed[j] = xj_perturbed;
            a_res(model, V, T, &x_perturbed)
        };
        let dadx_j = five_pt_stencil(&df_dxj, x[j], h);
        sum += x[j] * dadx_j;
    }

    // Helmholtz and Z contribution
    let a_res = a_res(model, V, T, x);
    let mu_res = a_res + (Z - 1.0) + dadx_k - sum;

    // Final fugacity coefficient
    (mu_res - Z.ln()).exp()
}

pub fn density_solver(x: &[f64], T: f64, P: f64, model: &PCSAFTModel) -> (f64, f64) {
    let tol = 1e-10;
    let dx = 1e-5;
    let n = x.len();

    // Physical bounds on eta
    let (eta_min, eta_max) = (1e-10, 0.7405);

    // Define residual function
    let f = |eta: f64| -> f64 {residual(model, T, P, x, 2, eta)}; 

    // Solve for eta s.t. Pcalc = Psys
    let eta_roots = root_solver(f, eta_min, eta_max, dx, tol);

    let eta_l = eta_roots[eta_roots.len() - 1];
    let eta_v = eta_roots[0];

    let rho_l = density(model, T, x, n, eta_l) * 1e30 / NA;
    let rho_v = density(model, T, x, n, eta_v) * 1e30 / NA;

   (rho_l, rho_v)
}

// Residual function Psys - Pcalc
pub fn residual(model: &PCSAFTModel, T: f64, P: f64, x: &[f64], n: usize, eta: f64) -> f64 {
    // Compute rho [Angstrom^-3]
    let rho = density(model, T, x, n, eta);

    // Compute zeta [Angstrom]
    let z0 = zeta_n(model, T, x, rho, 0);
    let z1 = zeta_n(model, T, x, rho, 1);
    let z2 = zeta_n(model, T, x, rho, 2);
    let z3 = zeta_n(model, T, x, rho, 3);

    // Estimate Z
    let Z = Z(model, T, x, n, rho, &(z0, z1, z2, z3));

    // Compute P
    let rho_SI = rho / NA * 1e30;
    let Pcalc = Z * rho_SI * R * T;

    // Residual
    P - Pcalc
}

// Density (Angstrom^-3)
pub fn density(model: &PCSAFTModel, T: f64, x: &[f64], n: usize, eta: f64) -> f64 {

    let sigma = &model.params.sigma;
    let epsilon = &model.params.epsilon;

    // Computing density
    let mut res = 0.0;
    for i in 0..n {
        let xi = x[i];
        let mi = &model.params.m[i];
        let di = diameter(sigma[i][i], epsilon[i][i], T);

        // Converting diameter to angstrom
        let di_angstrom = di * 1e10;

        res += xi * mi * di_angstrom.powi(3);
    }

    6.0/PI * eta/res
}

// Zeta(n), w/ n = {0, 1, 2, 3}
pub fn zeta_n(model: &PCSAFTModel, T: f64, x: &[f64], rho: f64, n: i32) -> f64 {

    let sigma = &model.params.sigma;
    let epsilon = &model.params.epsilon;
    
    let mut res = 0.0;
    for i in 0..x.len() {
        let xi = x[i];
        let mi = &model.params.m[i];
        let di = diameter(sigma[i][i], epsilon[i][i], T);

        // Converting diameter to angstrom
        let di_angstrom = di * 1e10;

        res += xi * mi * di_angstrom.powi(n);
    }
    PI/6.0 * rho * res
}

// Diameter (Angstrom)
pub fn diameter(sigma: f64, epsilon: f64, T: f64) -> f64 {
    sigma * (1.0 - 0.12*(-3.0*epsilon/T).exp())
}

// Compressibility factor
pub fn Z(model: &PCSAFTModel, T: f64, x: &[f64], n: usize, rho: f64, zeta: &(f64, f64, f64, f64)) -> f64 {
    1.0 + Z_hc(model, T, x, n, zeta) + Z_disp(model, T, x, n, rho, zeta)
}

// Hard-chain contribution
pub fn Z_hc(model: &PCSAFTModel, T: f64, x: &[f64], n: usize, zeta: &(f64, f64, f64, f64)) -> f64 {
    let m_bar = m_bar(model, x);
    let Z_hs = Z_hs(zeta);

    let mut res = 0.0;
    for i in 0..n {
        let xi = x[i];
        let mi = model.params.m[i];
        let g_hs_ii = g_hs_ij(model, T, zeta, i, i);
        let rho_dg_hs_ii = rho_dg_hs_ij(model, T, zeta, i, i);

        res += xi * (mi-1.0) * rho_dg_hs_ii / g_hs_ii;
    }

    m_bar * Z_hs - res
}

// Mean segment number
pub fn m_bar(model: &PCSAFTModel, x: &[f64]) -> f64 {
    let m = &model.params.m;
    x.iter().zip(m).map(|(xi, mi)| xi * mi).sum::<f64>() / x.iter().sum::<f64>()
}

// Hard-sphere contributino
pub fn Z_hs(zeta: &(f64, f64, f64, f64)) -> f64 {
    let (z0, z1, z2, z3) = zeta;

    let term1 = z3 / (1.0-z3);
    let term2 = 3.0 * z1 * z2 / (z0*(1.0-z3).powi(2));
    let term3 = z2.powi(3) * (3.0 - z3) / (z0*(1.0-z3).powi(3));

    term1 + term2 + term3
}

// Radial distribution function
pub fn g_hs_ij(model: &PCSAFTModel, T: f64, zeta: &(f64, f64, f64, f64), i: usize, j: usize) -> f64 {
    let (_, _, z2, z3) = zeta;

    let sigma_i = model.params.sigma[i][i];
    let epsilon_i = model.params.epsilon[i][i];
    let di = diameter(sigma_i, epsilon_i, T);
    let di_angstrom = di * 1e10;

    let sigma_j = model.params.sigma[j][j];
    let epsilon_j = model.params.epsilon[j][j];
    let dj = diameter(sigma_j, epsilon_j, T);
    let dj_angstrom = dj * 1e10;

    let dterm = di_angstrom * dj_angstrom / (di_angstrom + dj_angstrom);

    let term1 = 1.0 / (1.0-z3);
    let term2 = dterm * 3.0*z2 / (1.0-z3).powi(2);
    let term3 = dterm.powi(2) * 2.0*z2.powi(2) / (1.0-z3).powi(3);

    term1 + term2 + term3
}

// Radial distribution derivative
pub fn rho_dg_hs_ij(model: &PCSAFTModel, T: f64, zeta: &(f64, f64, f64, f64), i: usize, j: usize) -> f64 {
    let (_, _, z2, z3) = zeta;

    let sigma_i = model.params.sigma[i][i];
    let epsilon_i = model.params.epsilon[i][i];
    let di = diameter(sigma_i, epsilon_i, T);
    let di_angstrom = di * 1e10;

    let sigma_j = model.params.sigma[j][j];
    let epsilon_j = model.params.epsilon[j][j];
    let dj = diameter(sigma_j, epsilon_j, T);
    let dj_angstrom = dj * 1e10;

    let dterm = di_angstrom * dj_angstrom / (di_angstrom + dj_angstrom);

    let term1 = z3 / (1.0-z3).powi(2);
    let term2 = dterm * (
        3.0*z2 / (1.0-z3).powi(2) +
        6.0*z2*z3 / (1.0-z3).powi(3)
    );
    let term3 = dterm.powi(2) * (
        4.0*z2.powi(2) / (1.0-z3).powi(3) +
        6.0*z2.powi(2)*z3 / (1.0-z3).powi(4)
    );

    term1 + term2 + term3
}

// Dispersion term
pub fn Z_disp(model: &PCSAFTModel, T: f64, x: &[f64], _n: usize, rho: f64, zeta: &(f64, f64, f64, f64)) -> f64 {
    let (_, _, _, eta) = zeta;
    let m_bar = m_bar(model, x);

    let (m2eo3, m2e2o3) = bar_terms(model, x, T);

    let dI1 = dI(*eta, m_bar, 1);
    let dI2 = dI(*eta, m_bar, 2);

    let C1 = C1(m_bar, *eta);
    let C2 = C2(m_bar, *eta, C1);

    let I2 = I(*eta, m_bar, 2);

    let term1 = -2.0 * PI * rho * dI1 * m2eo3;
    let term2 = -PI * rho * m_bar * (C1 * dI2 + C2 * eta * I2) * m2e2o3;

    term1 + term2
}

// Derivative of I
pub fn dI(eta: f64, m_bar: f64, corr: usize) -> f64 {
     
     // 1st or 2nd order correlation
    let corr = if corr == 1 {
        &PCSAFT_CONSTS.aij
    } else if corr == 2 {
        &PCSAFT_CONSTS.bij
    } else {
        panic!("Invalid n value")
    };

    let m1 = (m_bar - 1.0) / m_bar;
    let m2 = (m_bar - 1.0) / m_bar * (m_bar - 2.0) / m_bar;

    let mut res = 0.0;
    for (j, &(c0, c1, c2)) in corr.iter().enumerate() {
        let kj = c0 + m1 * c1 + m2 * c2;
        res += kj * (j as f64 + 1.0) * eta.powi(j as i32);
    }

    res
}

// Segment average term
pub fn bar_terms(model: &PCSAFTModel, x: &[f64], T: f64) -> (f64, f64) {
    let n = x.len();
    
    let m = &model.params.m;
    let sigma = &model.params.sigma;
    let epsilon = &model.params.epsilon;

    let mut v1 = 0.0;
    let mut v2 = 0.0;

    for i in 0..n {
        for j in 0..n {
            let cst = x[i]*x[j]*m[i]*m[j]*sigma[i][j].powi(3);
            let exp = epsilon[i][j] / T;
            v1 += cst * exp;
            v2 += cst * exp.powi(2);
        }
    }

    let sum_z: f64 = x.iter().sum();
    let k = 1.0 / sum_z.powi(2);

    (k*v1, k*v2) 
}

// Compressibility expression
pub fn C1(m_bar: f64, eta: f64) -> f64 {
    let term1 = (8.0*eta - 2.0*eta.powi(2)) / (1.0-eta).powi(4);
    let term2 = (20.0*eta - 27.0*eta.powi(2) + 12.0*eta.powi(3) - 2.0*eta.powi(4)) / ((1.0-eta)*(2.0-eta)).powi(2);

    1.0 + m_bar*term1 + (1.0-m_bar)*term2
}

// Derivative of C1 above
pub fn C2(m_bar: f64, eta: f64, C1: f64) -> f64 {
    let term1 = (-4.0*eta.powi(2) + 20.0*eta + 8.0) / (1.0-eta).powi(5);
    let term2 = (2.0*eta.powi(3) + 12.0*eta.powi(2) - 48.0*eta + 40.0) / ((1.0-eta)*(2.0-eta)).powi(3);

    -C1.powi(2) * (m_bar*term1 + (1.0-m_bar)*term2)
}

// I function
pub fn I(eta: f64, m_bar: f64, corr: usize) -> f64 {
    
    // 1st or 2nd order correlation
    let corr = if corr == 1 {
        &PCSAFT_CONSTS.aij
    } else if corr == 2 {
        &PCSAFT_CONSTS.bij
    } else {
        panic!("Invalid n value")
    };

    let m1 = (m_bar - 1.0) / m_bar;
    let m2 = (m_bar - 1.0) / m_bar * (m_bar - 2.0) / m_bar;

    let mut res = 0.0;
    for (j, &(c0, c1, c2)) in corr.iter().enumerate() {
        let kj = c0 + m1 * c1 + m2 * c2;
        res += kj * eta.powi(j as i32);
    }

    res
}

use crate::structures::*;
use crate::constants::*;
use crate::mixing::*;
use crate::a_assoc::*;
use crate::a_disp::*;
use crate::a_hc::*;

// Residual Helmholtz energy
// a_res = a_hc + a_disp + a_assoc
pub fn a_res(model: &PCSAFTModel, V: f64, T: f64, z: &[f64]) -> f64 {
    let data: (Vec<f64>, f64, f64, f64, f64, f64) = prep_data(model, V, T, z);

    let a_hc = a_hc(model, V, T, z, &data);
    let a_disp = a_disp(model, V, T, z, &data);
    let a_assoc = a_assoc(model, V, T, z, &data);

    a_hc + a_disp + a_assoc
}

// Model builder
pub fn build_model(components: Vec<PCSAFTComponent>, alpha: f64, atol: f64, rtol: f64, max_iters: usize, mixing: &str) -> PCSAFTModel {
    let n = components.len();

    // Scalar component parameters
    let mw: Vec<f64> = components.iter().map(|c| c.mw).collect();
    let m: Vec<f64> = components.iter().map(|c| c.m).collect();
    let sigma_vec: Vec<f64> = components.iter().map(|c| c.sigma * 1e-10).collect();
    let epsilon_vec: Vec<f64> = components.iter().map(|c| c.epsilon).collect();
    let bondvol_diag: Vec<f64> = components.iter().map(|c| c.bondvol).collect();
    let epsilon_assoc_diag: Vec<f64> = components.iter().map(|c| c.epsilon_assoc).collect();

    // Association site structure
    let site_labels: Vec<Vec<String>> = components.iter().map(|c| c.sites.clone()).collect();
    let mut flat_sites = Vec::new();
    let mut p = Vec::with_capacity(site_labels.len() + 1);
    let mut count = 0;
    for comp_sites in &site_labels {
        p.push(count);
        count += comp_sites.len();
        for _ in comp_sites {
            flat_sites.push(1);
        }
    }
    p.push(count);

    let site_param = SiteParam {
        components: (0..n).map(|i| format!("comp_{}", i)).collect(),
        sites: site_labels.clone(),
        n_sites: PackedVector { v: flat_sites, p },
    };

    // Full association matrix
    let mut outer_indices = Vec::new();
    let mut inner_indices = Vec::new();
    let mut bondvol_values = Vec::new();
    let mut epsilon_assoc_values = Vec::new();

    for i in 0..n {
        let n_sites = site_labels[i].len();
        for a in 0..n_sites {
            for b in 0..n_sites {
                outer_indices.push((i, i));
                inner_indices.push((a, b));
                bondvol_values.push(bondvol_diag[i]);
                epsilon_assoc_values.push(epsilon_assoc_diag[i]);
            }
        }
    }

    let bondvol = AssocMatrix::from_parts(bondvol_values, outer_indices.clone(), inner_indices.clone());
    let epsilon_assoc = AssocMatrix::from_parts(epsilon_assoc_values, outer_indices, inner_indices);

    // Construct base model
    let mut model = PCSAFTModel {
        params: PCSAFTParams {
            mw,
            m,
            sigma: vec![vec![0.0; n]; n],
            epsilon: vec![vec![0.0; n]; n],
            bondvol,
            epsilon_assoc,
            kij: None,
            lij: None,
        },
        assoc: AssocOptions {
            sites: site_param,
            alpha,
            atol,
            rtol,
            max_iters,
            combining: mixing.to_string(),
        },
    };

    recombine_saft(&mut model, sigma_vec, epsilon_vec, None, None);

    model

}

// Restructuring model
pub fn recombine_saft(model: &mut PCSAFTModel, sigma: Vec<f64>, epsilon: Vec<f64>, k: Option<&[Vec<f64>]>, l: Option<&[Vec<f64>]>) {    
    let (sigma_mix, epsilon_mix) = lorentz_berthelot(&sigma, &epsilon, k, l);

    model.params.sigma = sigma_mix;
    model.params.epsilon = epsilon_mix;

    // Update association scheme
    if assoc_pair_length(model) == 0 {return;}

    let sigma = &model.params.sigma;
    let bondvol = &model.params.bondvol;
    let epsilon_assoc = &model.params.epsilon_assoc;
    let sites = &model.assoc.sites.sites;
    let assoc_options = &model.assoc;

    let (new_bondvol, new_epsilon_assoc) = assoc_mix(bondvol, epsilon_assoc, sigma, sites, assoc_options);

    model.params.bondvol = new_bondvol;
    model.params.epsilon_assoc = new_epsilon_assoc;
}

// Data handler
pub fn prep_data(model: &PCSAFTModel, V: f64, T: f64, z: &[f64])-> (Vec<f64>, f64, f64, f64, f64, f64) {
    let d = ck_diameter(model, V, T, z, 0.12, 3.00);
    let (z0, z1, z2, z3) = zeta(model, V, z, &d);
    let m = &model.params.m;
    let m_bar = z.iter().zip(m).map(|(zi, mi)| zi * mi).sum::<f64>() / z.iter().sum::<f64>();
    
    (d, z0, z1, z2, z3, m_bar)
}

// Chen-Kregleswski diameter
pub fn ck_diameter(model: &PCSAFTModel, _V: f64, T: f64, z: &[f64], k1: f64, k2: f64) -> Vec<f64> {
    let mut d = vec![0.0; z.len()];
    
    for i in 0..z.len() {
        let sigma = model.params.sigma[i][i];
        let epsilon = model.params.epsilon[i][i];
        d[i] = sigma * (1.0 - k1 * (-k2*epsilon / T).exp());
    }

    d
}

// Atomic packing fractions
pub fn zeta(model: &PCSAFTModel, V: f64, z: &[f64], d: &[f64]) -> (f64, f64, f64, f64) {
    let mut z0 = 0.0;
    let mut z1 = 0.0;
    let mut z2 = 0.0;
    let mut z3 = 0.0;

    for i in 0..z.len() {
        let m_i = model.params.m[i];
        let d_i = d[i];
        let xS = z[i] * m_i;
        z0 += xS;
        z1 += xS * d_i;
        z2 += xS * d_i.powi(2);
        z3 += xS * d_i.powi(3);
    }

    let c = PI/6.0 * NA/V;
    (c * z0, c * z1, c * z2, c * z3)
}

use crate::structures::*;

// Lorentz-Berthelot mixing rules
pub fn lorentz_berthelot(sigma: &[f64], epsilon: &[f64], k: Option<&[Vec<f64>]>, l: Option<&[Vec<f64>]>) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = sigma.len();
    let mut sigma_ij = vec![vec![0.0; n]; n];
    let mut epsilon_ij = vec![vec![0.0;n]; n];

    for i in 0..n {
        for j in 0..n {
            let kij = k.map_or(0.0, |mat| mat[i][j]);
            let lij = l.map_or(0.0, |mat| mat[i][j]);

            sigma_ij[i][j] = 0.5 * (sigma[i]+sigma[j]) + lij;
            epsilon_ij[i][j] = (epsilon[i]*epsilon[j]).sqrt() * (1.0-kij)
        }
    }

    (sigma_ij, epsilon_ij)
}

// Mixing rules for association
pub fn assoc_mix(bondvol: &AssocMatrix<f64>, epsilon_assoc: &AssocMatrix<f64>, sigma: &Vec<Vec<f64>>, sites: &Vec<Vec<String>>, assoc_options: &AssocOptions) -> (AssocMatrix<f64>, AssocMatrix<f64>) {
    let method = assoc_options.combining.as_str();

    // ----- CASE 1: Use zero_mix -----
    if method == "ElliottRuntime" || method == "EsdRuntime" {
        let bv_mixed = zero_mix(bondvol, sites);
        let eps_mixed = zero_mix(epsilon_assoc, sites);
        (bv_mixed, eps_mixed)

    // ----- CASE 2: Use mixing rules -----
    } else if method == "Elliott" || method == "Esd" {
        let bv_mixed = bondvol_mix(bondvol, sigma, sites);
        let eps_mixed = epsilon_assoc_mix(epsilon_assoc, sites);
        (bv_mixed, eps_mixed)

    } else {
        panic!("Unknown association combining rule: {}", method);
    }
}

// Zero mixing rule: fill missing off-diagonal values using geometric mean and reset to zero if nonzero.
pub fn zero_mix(param: &AssocMatrix<f64>, sites: &Vec<Vec<String>>) -> AssocMatrix<f64> {
    if param.values.is_empty() {
        return param.clone();
    }

    let mut extended = assoc_extend_matrix(param, sites);
    let sentinel = -124.0;

    for (k, &(i, j)) in extended.outer_indices.iter().enumerate() {
        let (a, b) = extended.inner_indices[k];
        if extended.values[k] == 0.0 && valid_site_comb(sites, i, j, a, b) {
            let dij = (param.get(i, i, a, b) * param.get(j, j, a, b)).sqrt();
            if dij != 0.0 {
                extended.values[k] = sentinel;
            }
        }
    }

    // Convert sentinel back to 0
    for val in &mut extended.values {
        if *val == sentinel {
            *val = 0.0;
        }
    }

    extended
}

// Mixing rule for bond volume using geometric mean and sigma scaling.
pub fn bondvol_mix(bondvol: &AssocMatrix<f64>, sigma: &Vec<Vec<f64>>, sites: &Vec<Vec<String>>) -> AssocMatrix<f64> {
    if bondvol.values.is_empty() {
        return bondvol.clone();
    }

    let mut extended = assoc_extend_matrix(bondvol, sites);

    for (k, &(i, j)) in extended.outer_indices.iter().enumerate() {
        let (a, b) = extended.inner_indices[k];
        if extended.values[k] == 0.0 && valid_site_comb(sites, i, j, a, b) {
            let v = bondvol.get(i, i, a, b) * bondvol.get(j, j, a, b);
            let s = (sigma[i][i] * sigma[j][j]).sqrt() / sigma[i][j];
            extended.values[k] = v.sqrt() * s.powi(3);
        }
    }

    extended
}

// Mixing rule for association energy using arithmetic mean of diagonal terms.
pub fn epsilon_assoc_mix(epsilon_assoc: &AssocMatrix<f64>, sites: &Vec<Vec<String>>) -> AssocMatrix<f64> {
    if epsilon_assoc.values.is_empty() {
        return epsilon_assoc.clone();
    }

    let mut extended = assoc_extend_matrix(epsilon_assoc, sites);

    for (k, &(i, j)) in extended.outer_indices.iter().enumerate() {
        let (a, b) = extended.inner_indices[k];
        if extended.values[k] == 0.0 && valid_site_comb(sites, i, j, a, b) {
            extended.values[k] = 0.5 * (epsilon_assoc.get(i, i, a, b) + epsilon_assoc.get(j, j, a, b));
        }
    }

    extended
}


// ---------- HELPER METHODS ----------

// Check if site indices (a, b) are valid for components i and j.
pub fn valid_site_comb(sites: &Vec<Vec<String>>, i: usize, j: usize, a: usize, b: usize) -> bool {
    if sites[i].is_empty() || sites[j].is_empty() {
        return false;
    }
    a < sites[i].len() && b < sites[j].len()
}

/// Expand association matrix to include all valid (i, j, a, b) site combinations.
pub fn assoc_extend_matrix(mat: &AssocMatrix<f64>, sites: &Vec<Vec<String>>) -> AssocMatrix<f64> {
    if mat.values.is_empty() {
        return mat.clone();
    }

    let comps = sites.len();
    let mut idx = Vec::new();

    for i in 0..comps {
        for j in 0..=i {
            let la = sites[i].len();
            let lb = sites[j].len();
            if la != 0 && lb != 0 {
                for a in 0..la {
                    let start = if i == j { a } else { 0 };
                    for b in start..lb {
                        idx.push((i, j, a, b));
                    }
                }
            }
        }
    }

    idx.sort();
    let mut extended_vals = vec![0.0; idx.len()];
    for (k, &(i, j, a, b)) in idx.iter().enumerate() {
        extended_vals[k] = mat.get(i, j, a, b);
    }
    let outer_indices: Vec<(usize, usize)> = idx.iter().map(|&(i, j, _, _)| (i, j)).collect();
    let inner_indices: Vec<(usize, usize)> = idx.iter().map(|&(_, _, a, b)| (a, b)).collect();

    AssocMatrix::from_parts(extended_vals, outer_indices, inner_indices)
}

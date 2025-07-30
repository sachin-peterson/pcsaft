use crate::structures::*;
use crate::constants::*;
use crate::utils::*;
use crate::solver::*;

// Association term
pub fn a_assoc(model: &PCSAFTModel, V: f64, T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> f64 {
    
    // Number of sites
    let n = assoc_pair_length(model);
    
    // No association sites
    if n == 0 {return 0.0}

    // 1 association site
    if n == 1 {return compute_a_assoc_exact1(model, V, T, z, data)};

    // 2+ association sites
    let (X, delta) = assoc_fraction(model, V, T, z, data);
    compute_a_assoc(model, V, T, z, &X, delta)
}

// Radial distribution function
pub fn g_hs(_model: &PCSAFTModel, i: usize, j: usize, data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> f64 {
    let (d, _z0, _z1, z2, z3, _m_bar) = data;
    let d_i = d[i];
    let d_j = d[j];

    let c1 = 1.0 / (1.0-z3);
    let c2 = 3.0 * z2 / (1.0-z3).powi(2);
    let c3 = 2.0 * z2.powi(2) / (1.0-z3).powi(3);

    c1 + d_i*d_j/(d_i+d_j) * c2 + (d_i*d_j/(d_i+d_j)).powi(2) * c3
}


// ---------- 1 SITE ----------

// Association energy for 1 site
pub fn compute_a_assoc_exact1(model: &PCSAFTModel, V: f64, T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> f64 {
    let (xia, xjb, i, j, a, b, _n_sites_flat, _idxs, _delta) = X_exact1(model, V, T, z, data);

    let n = &model.assoc.sites.n_sites;
    let nia = n.get(i, a) as f64;
    let njb = n.get(j, b) as f64;

    let mut res = z[i] * nia * (xia.ln() - 0.5*xia + 0.5);
    if i != j || a != b {
        res += z[j] * njb * (xjb.ln() - 0.5*xjb + 0.5);
    }

    res / z.iter().sum::<f64>()
}

// Association fraction for 1 site
pub fn X_exact1(model: &PCSAFTModel, V: f64, T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> (f64, f64, usize, usize, usize, usize, usize, Vec<usize>, f64) {
    let kappa = assoc_structure(model);
    let (i, j) = kappa.outer_indices[0];
    let (a, b) = kappa.inner_indices[0];

    let delta = delta_ijab(model, V, T, z, i, j, a, b, data);

    let sites = &model.assoc.sites.n_sites;
    let idxs = sites.p.clone();
    let n = (sites.v).len();

    let rho = NA/V;
    let z_i = z[i];
    let z_j = z[j];

    let nia = sites.get(i, a) as f64;
    let njb = sites.get(j, b) as f64;

    let kia = nia * z_i * rho * delta;
    let kjb = njb * z_j * rho * delta;

    let a_ = kia;
    let b_ = 1.0 - kia - kjb;
    let c_ = -1.0;

    let xia = -2.0*c_ / (b_ + (b_*b_ - 4.0*a_*c_).sqrt());
    let xjb = (1.0 - kia*xia) / (1.0 - (kia*xia).powi(2));

    (xia, xjb, i, j, a, b, n, idxs, delta)
}

// Association strength for 1 site
pub fn delta_ijab(model: &PCSAFTModel, _V: f64, T: f64, _z: &[f64], i: usize, j: usize, a: usize, b: usize, data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> f64 {
    let epsilon_assoc = &model.params.epsilon_assoc;
    let kappa = &model.params.bondvol;
    let sigma = &model.params.sigma;

    let mut kappa_ijab = 0.0;
    let mut epsilon_ijab = 0.0;
    let mut found = false;

    for idx in 0..kappa.values.len() {
        if kappa.outer_indices[idx] == (i, j) && kappa.inner_indices[idx] == (a, b) {
            kappa_ijab = kappa.values[idx];
            epsilon_ijab = epsilon_assoc.values[idx];
            found = true;
            break;
        }
    }

    if !found || kappa_ijab == 0.0 {
        return 0.0;
    }

    let g_ij = g_hs(model, i, j, data);

    g_ij * sigma[i][j].powi(3) * (epsilon_ijab / T).exp_m1() * kappa_ijab
}

// Association fraction & strength for 1 site 
pub fn X_and_delta_exact1(model: &PCSAFTModel, V: f64, T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> (PackedVecVec<f64>, AssocMatrix<f64>) {
    let (xia, xjb, i, j, a, b, n, idxs, delta_ijab) = X_exact1(model, V, T, z, data);
    let X = pack_X_exact1(xia, xjb, i, j, a, b, n, idxs);

    let mut delta = assoc_similar(model);
    delta.values[0] = delta_ijab;

    (X, delta)
}


// ---------- 2+ SITES ----------

// Association energy for 2+ sites
pub fn compute_a_assoc(model: &PCSAFTModel, _V: f64, _T: f64, z: &[f64], X: &PackedVecVec<f64>, _delta: AssocMatrix<f64>) -> f64 {
    let sites = &model.assoc.sites;
    let n = &sites.n_sites;

    let mut res = 0.0;
    for i in 0..z.len() {
        let n_i = get_component_sites(n, i);
        let z_i = z[i];
        if z_i == 0.0 || n_i.is_empty() {continue;}

        let x_i = &X.v[i];
        let mut res_inner = 0.0;

        for a in 0..n_i.len() {
            let n_ia = n_i[a] as f64;
            let x_ia = x_i[a];
            res_inner += n_ia * (x_ia.ln() - 0.5*x_ia + 0.5);
        }

        res += res_inner*z_i;
    }

    res / z.iter().sum::<f64>()

}

// Association fraction for 2+ sites
pub fn assoc_fraction(model: &PCSAFTModel, V: f64, T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> (PackedVecVec<f64>, AssocMatrix<f64>) {
    let nn: usize = assoc_pair_length(model);
    let zero_matrix = assoc_similar(model);
    
    // If only 1 association site
    if nn == 1 {
        let (xia, xjb, i, j, a, b, n, idxs, _delta) = X_exact1(model, V, T, z, data);
        let X = pack_X_exact1(xia, xjb, i, j, a, b, n, idxs);
        return (X, zero_matrix)
    }

    // General case (>1 sites)
    let (X, delta) = compute_X_and_delta(model, V, T, z, data);
    (X, delta)
}

// Association strength for 2+ sites
pub fn delta(model: &PCSAFTModel, _V: f64, T: f64, _z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> AssocMatrix<f64> {
    let epsilon_assoc = &model.params.epsilon_assoc;
    let kappa = &model.params.bondvol;
    let sigma = &model.params.sigma;

    // Shape compatible zero matrix
    let mut delta = assoc_similar(model);

    for idx in 0..delta.values.len() {
        let (i, j) = delta.outer_indices[idx];
        let (a, b) = delta.inner_indices[idx];

        let mut found = false;
        let mut kappa_ijab = 0.0;
        let mut epsilon_ijab = 0.0;

        for k in 0..kappa.values.len() {
            if kappa.outer_indices[k] == (i, j) && kappa.inner_indices[k] == (a, b) {
                kappa_ijab = kappa.values[k];
                epsilon_ijab = epsilon_assoc.values[k];
                found = true;
                break;
            }
        }

        if !found || kappa_ijab == 0.0 {
            delta.values[idx] = 0.0;
            continue;
        }

        let g_ij = g_hs(model, i, j, data);
        let exp_term = (epsilon_ijab / T).exp_m1();

        delta.values[idx] = g_ij * sigma[i][j].powi(3) * exp_term * kappa_ijab;
    }
    delta
}

// Association fraction & strength for 2+ sites
pub fn compute_X_and_delta(model: &PCSAFTModel, V: f64, T: f64, z: &[f64], data: &(Vec<f64>, f64, f64, f64, f64, f64)) -> (PackedVecVec<f64>, AssocMatrix<f64>) {
    let nn: usize = assoc_pair_length(model);
    if nn == 1 {return X_and_delta_exact1(model, V, T, z, data);}

    let options = &model.assoc;
    let delta = delta(model, V, T, z, data);
    let K = assoc_site_matrix(model, V, T, z, data, &delta);
    let idxs = model.assoc.sites.n_sites.p.clone();
    let Xsol = assoc_matrix_solve(&K, &options);
    let packed_X = PackedVecVec::new(idxs, Xsol);

    (packed_X, delta)
}


// ---------- MATRIX METHODS ----------

// Association strength matrix, K
pub fn assoc_site_matrix(model: &PCSAFTModel, V: f64, _T: f64, z: &[f64], _data: &(Vec<f64>, f64, f64, f64, f64, f64), delta: &AssocMatrix<f64>) -> Vec<Vec<f64>> {
    let sites = &model.assoc.sites.n_sites;
    let p = &sites.p;
    let n = &sites.v;
    let rho = NA/V;

    let delta_vals = &delta.values;

    let nn = n.len();
    // Allocate dense matrix K
    let mut K = vec![vec![0.0; nn]; nn];

    let ii = &delta.outer_indices;
    let aa = &delta.inner_indices;
    let idxs = ii.len();

    // Populating matrix
    for i in 0..z.len() {
        let site_start = p[i];
        let site_end = p[i+1];
        for a in 0..(site_end-site_start) {
            let ia = compute_idx(p, i, a);
            for idx in 0..idxs {
                let ij = ii[idx];
                let ab = aa[idx];
                if issite(i, a, ij, ab) {
                    let j = complement_idx(i, ij);
                    let b = complement_idx(a, ab);
                    let jb = compute_idx(p, j, b);
                    let njb = n[jb] as f64;
                    let zj = z[j];
                    if zj != 0.0 {K[ia][jb] = rho*njb*zj*delta_vals[idx]}
                }
            }
        }
    }
    K
}

// Solving the association matrix
pub fn assoc_matrix_solve(K: &Vec<Vec<f64>>, options: &AssocOptions) -> Vec<f64> {
    let alpha = options.alpha;
    let atol = options.atol;
    let rtol = options.rtol;
    let max_iters = options.max_iters;

    assoc_matrix_solve_full(K, alpha, atol, rtol, max_iters)
}

// Find vector s.t. Ax .* x + x - 1 = 0
pub fn assoc_matrix_solve_full(K: &Vec<Vec<f64>>, alpha: f64, atol: f64, rtol: f64, max_iters: usize) -> Vec<f64> {
    let n = K.len();
    let mut X0 = vec![0.0; n];
    let success = assoc_matrix_X0(K, &mut X0);
    if success {return X0;}

    let mut Xsol = X0.clone();
    let n_iters = 5*n;

    // Successive substitution
    for _ in 0..n_iters {
        partial_matmul(fx,&mut Xsol, K, &X0, alpha);
        let (converged, finite) = check_convergence(&Xsol, &X0, atol, rtol, false);
        if converged {
            if !finite {
                Xsol.fill(f64::NAN);
            }
            return Xsol
        }
        X0.copy_from_slice(&Xsol);
    }

    // Nonlinear Newton
    let mut H = vec![vec![0.0; n]; n];
    let mut pivot = vec![0; n];
    let mut dX = Xsol.clone();
    let mut KX = vec![0.0; n];

    for _ in n_iters..max_iters {
        matrix_vector_mul(&mut KX, K, &Xsol);
        for i in 0..n {
            for j in 0..n {
                H[i][j] = -K[i][j];
            }
            H[i][i] -= (1.0+KX[i]) / Xsol[i];
        }

        // Gradient
        for i in 0..n {
            dX[i] = 1.0 / Xsol[i] - 1.0 - KX[i];
        }

        // LU solver
        let (A, pivots, _) = lu_solver(&mut H, &mut pivot);
        ldiv(&A, pivots, &mut dX);
        X0.copy_from_slice(&Xsol);

        for k in 0..n {
            let Xk = Xsol[k];
            let dXk = dX[k];
            let X_newton = Xk - dXk;

            // SS step
            if !(0.0..=1.0).contains(&X_newton) {
                Xsol[k] = 0.5 * (Xk + X0[k]);
            } else {
                Xsol[k] = X_newton;
            }
        }

        let (converged, finite) = check_convergence(&Xsol, &X0, atol, rtol, false);
        if converged {
            if !finite {
                Xsol.fill(f64::NAN);
            }
            return Xsol;
        }
    }

    Xsol.fill(f64::NAN);
    Xsol
}

// Solving Ax * x + x - 1 = 0
pub fn assoc_matrix_X0(K: &Vec<Vec<f64>>, X: &mut Vec<f64>) -> bool {
    let success = false;

    let n = K.len();
    if n != K[0].len() {
        panic!("Matrix K must be square");
    }

    if n == 1{
        let k = K[0][0];
        X[0] = 0.5 * (-1.0 + (1.0 + 4.0 * k).sqrt());
        return true;
    } 
    
    if check_antidiagonal2x2(K) {
        X_exact2(K, X);
        return true;
    }

    if check_antidiagonal4x4(K) {
        // Views for K11, K12, K21, K22
        let k11 = vec![vec![K[0][0], K[0][1]], vec![K[1][0], K[1][1]]];
        let k12 = vec![vec![K[0][2], K[0][3]], vec![K[1][2], K[1][3]]];
        let k21 = vec![vec![K[2][0], K[2][1]], vec![K[3][0], K[3][1]]];
        let k22 = vec![vec![K[2][2], K[2][3]], vec![K[3][2], K[3][3]]];

        X_exact2(&k11, &mut X[0..2]);
        X_exact2(&k22, &mut X[2..4]);

        if is_zero_block(&k12) && is_zero_block(&k21) {
            return true;
        } else {
            return false;
        }
    }

    let (kmin, kmax) = nonzero_extrema(K);
    let f = if kmax > 1.0 { 1.0 / kmin } else { 1.0 - kmin };

    for x in X.iter_mut() {
        *x = f.min(1.0);
    }

    success
}


// ---------- HELPER METHODS ----------

// No. of site combinations
pub fn assoc_pair_length(model: &PCSAFTModel) -> usize {
    let val = assoc_structure(model);
    val.values.len()
}

// Structure of association matrix
pub fn assoc_structure(model: &PCSAFTModel) -> AssocMatrix<usize> {
     let matrix = &model.params.bondvol;
     let l = matrix.values.len();
     
     AssocMatrix { 
        values: (0..l).collect(), 
        outer_indices: matrix.outer_indices.clone(), 
        inner_indices: matrix.inner_indices.clone(),
        outer_size: matrix.outer_size,
        inner_size: matrix.inner_size 
    }
}

// Zero-initialized association matrix
pub fn assoc_similar(model: &PCSAFTModel) -> AssocMatrix<f64> {
    let sites = &model.assoc.sites.sites;
    let mut values = Vec::new();
    let mut outer_indices = Vec::new();
    let mut inner_indices = Vec::new();

    for (i, sites_i) in sites.iter().enumerate() {
        for (j, sites_j) in sites.iter().enumerate() {
            for (a, _) in sites_i.iter().enumerate() {
                for (b, _) in sites_j.iter().enumerate() {
                    values.push(0.0);
                    outer_indices.push((i, j));
                    inner_indices.push((a, b));
                }
            }
        }
    }

    let outer_size = (sites.len(), sites.len());
    let inner_size = (
        sites.iter().map(|s| s.len()).max().unwrap_or(0),
        sites.iter().map(|s| s.len()).max().unwrap_or(0),
    );

    AssocMatrix {
        values,
        outer_indices,
        inner_indices,
        outer_size,
        inner_size,
    }
}

// Pack solution for a 1 site pair
pub fn pack_X_exact1(xia: f64, xjb: f64, i: usize, j: usize, a: usize, b: usize, n: usize, idxs: Vec<usize>) -> PackedVecVec<f64> {
    let xsol = vec![1.0; n];
    let mut packed = PackedVecVec::new(idxs, xsol);
    packed.v[j][b] = xjb;
    packed.v[i][a] = xia;
    packed
}

// List of site counts for component i
pub fn get_component_sites(n: &PackedVector<usize>, i: usize) -> &[usize] {
    let start = n.p[i];
    let end = if i + 1 < n.p.len() {
        n.p[i + 1]
    } else {
        n.v.len()
    };
    &n.v[start..end]
}

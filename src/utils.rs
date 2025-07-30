// Complement index
pub fn complement_idx(i: usize, ij: (usize, usize)) -> usize {
    let (i1, i2) = ij;
    if i1 == i {i2} else {i1}
}

// Compute index
pub fn compute_idx(idxs: &[usize], i: usize, a: usize) -> usize {
    return idxs[i] + a;
}

// Check site validity
pub fn issite(i: usize, a: usize, ij: (usize, usize), ab: (usize, usize)) -> bool {
    let ia = (i, a);
    let ia1 = (ij.0, ab.0);
    let ia2 = (ij.1, ab.1);
    
    ia == ia1 || ia == ia2
}

// Check 2x2 antidiagonal
pub fn check_antidiagonal2x2(x: &[Vec<f64>]) -> bool {
    if x.len() != 2 || x[0].len() != 2 {
        return false;
    }

    let x11 = x[0][0];
    let x12 = x[0][1];
    let x21 = x[1][0];
    let x22 = x[1][1];

    x11 == 0.0 && x22 == 0.0 && x12 >= 0.0 && x21 >= 0.0
}

// Check 4x4 antidigonals
pub fn check_antidiagonal4x4(x: &[Vec<f64>]) -> bool {
    if x.len() != 4 || x[0].len() != 4 {
        return false;
    }

    check_antidiagonal2x2(&[vec![x[0][0], x[0][1]], vec![x[1][0], x[1][1]]]) &&
    check_antidiagonal2x2(&[vec![x[0][2], x[0][3]], vec![x[1][2], x[1][3]]]) &&
    check_antidiagonal2x2(&[vec![x[2][0], x[2][1]], vec![x[3][0], x[3][1]]]) &&
    check_antidiagonal2x2(&[vec![x[2][2], x[2][3]], vec![x[3][2], x[3][3]]])
}

// Analytical solution for 2x2 association matrix
pub fn X_exact2(K: &[Vec<f64>], X: &mut [f64]) {
    let k1 = K[0][1];
    let k2 = K[1][0];

    let a = k2;
    let b = 1.0 - k2 + k1;
    let c = -1.0;
    let denom = b + (b * b - 4.0 * a * c).sqrt();

    let x1 = -2.0 * c / denom;
    let x1k = k2 * x1;
    let x2 = (1.0 - x1k) / (1.0 - x1k * x1k);

    X[0] = x1;
    X[1] = x2;
}

// Smallest and largest nonzero matrix values
pub fn nonzero_extrema(K: &[Vec<f64>]) -> (f64, f64) {
    let mut min_val = 0.0;
    let mut max_val = 0.0;

    for row in K {
        for &k in row {
            if k > max_val {
                max_val = k;
            }
            if min_val == 0.0 {
                min_val = k;
            } else if k != 0.0 {
                min_val = min_val.min(k);
            }
        }
    }

    (min_val, max_val)
}

// Check if a matrix is entirely zeroes
pub fn is_zero_block(K: &[Vec<f64>]) -> bool {
    K.iter().flatten().all(|&x| x == 0.0)
}

use crate::newton::*;

// Function Ax * x + x - 1 = 0
pub fn fx(alpha: f64, kx: f64, x: f64) -> f64 {
    alpha/(1.0+kx) + (1.0-alpha)*x
}

// Partial matrix multiplication
pub fn partial_matmul<F>(fx: F, out: &mut [f64], K: &[Vec<f64>], x: &[f64], alpha: f64) where F: Fn(f64, f64, f64) -> f64 {
    let n = x.len();
    for i in 0..n {
        let mut kxi = 0.0;
        let ki = &K[i];

        for j in 0..i {
            kxi += ki[j] * out[j];
        }

        for j in 0..n {
            kxi += ki[j] * x[j];
        }

        out[i] = fx(alpha, kxi, x[i]);
    }
}

// Matrix vector multiplication
pub fn matrix_vector_mul(out: &mut Vec<f64>, K: &Vec<Vec<f64>>, x: &Vec<f64>) {
    let n = x.len();
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += K[i][j] * x[j];
        }
        out[i] = sum;
    }
}

// Lower-upper solver, via NL Newton w/ successive substitution
pub fn lu_solver(A: &mut Vec<Vec<f64>>, pivots: &mut Vec<usize>) -> (Vec<Vec<f64>>, Vec<usize>, usize) {
    let n = A.len();
    let m = A[0].len();
    let min_mn = n.min(m);
    let mut info = 0;

    // Pivoting
    for k in 0..min_mn {
        let mut kp = k;
        let mut a_max = A[k][k].abs();
        for i in (k+1)..n {
            let abs_i = A[i][k].abs();
            if abs_i > a_max {
                kp = i;
                a_max = abs_i;
            }
        }
        pivots[k] = kp;

        if A[kp][k] != 0.0 {
            if kp != k {
                A.swap(k, kp);
            }

            // Scaling
            let Akk_inv = 1.0 / A[k][k];
            for i in k+1..m {
                A[i][k] *= Akk_inv;
            }

        } else if info == 0 {
            info = k+1;
        }

        // Update rest
        for j in (k+1)..m {
            for i in (k+1)..n {
                A[i][j] -= A[i][k]*A[k][j];
            }
        }
    }

    // Return LU result
    (A.clone(), pivots.clone(), info)
}

// Solving the linear system using LU decomposition
pub fn ldiv(A: &Vec<Vec<f64>>, piv: Vec<usize>, b: &mut [f64]) {
    let n = b.len();
    let L = &A;
    let U = &A;

    // Apply pivoting
    for i in 0..n {
        if piv[i] != i {
            b.swap(i, piv[i]);
        }
    }

    // Forward substitution: solve L*y = b
    for i in 0..n {
        for j in 0..i {
            b[i] -= L[i][j] * b[j];
        }
    }

    // Backward substitution: solve U*x = y
    for i in (0..n).rev() {
        for j in (i + 1)..n {
            b[i] -= U[i][j] * b[j];
        }
        b[i] /= U[i][i];
    }
}

// Checking convergence for NL Newton method
pub fn check_convergence(xsol: &[f64], xi: &[f64], atol: f64, rtol: f64, lognorm: bool) -> (bool, bool) {
    // Check for NaNs or infinities
    if xi.iter().any(|&x| !x.is_finite()) {
        return (true, false);
    }

    // Check for identity xi == xsol
    if xi.iter().zip(xsol.iter()).all(|(&a, &b)| a == b) {
        return (true, true);
    }

    // Compute dx and norm
    let mut norm_x = 0.0 as f64;
    let mut norm_xi = 0.0 as f64;

    let n = xi.len();

    for i in 0..n {
        let dx = if lognorm {
            xi[i] - xsol[i]
        } else {
            xi[i] / xsol[i] - 1.0
        };

        let abs_dx = dx.abs();
        if abs_dx > norm_x {
            norm_x = abs_dx;
        }

        if !lognorm {
            let abs_xi = xi[i].abs();
            if abs_xi > norm_xi {
                norm_xi = abs_xi;
            }
        }
    }

    if norm_x.abs() < f64::max(atol, norm_xi*rtol) {
        (true, true)
    } else {
        (false, false)
    }
}


// ROOT SOLVING SCHEME FOR P(ALPHA)

// Main solver for P(alpha)
pub fn root_solver_p_alpha<F>(fun: F,xi: f64, xf: f64, dx: f64, tol: f64, z: &[f64], k: &[f64]) -> f64 where F: Fn(f64) -> f64 {
    // Offset initial bound
    let mut xl = xi - dx / 2.0;
    let mut xr = xl + dx;

    // Offset final bound
    let xft = xf + dx / 2.0;

    // Initialize storage
    let mut r = 0.0;

    // Looping until the end of the interval
    while xl <= xft {
        let mut yl = fun(xl);
        let mut yr = fun(xr);

        // Bracketing a root or singularity
        while ((yr * yl) > 0.0 || (yr * yl).is_infinite() || (yr * yl).is_nan()) && xr <= xft {
            xl = xr;
            xr = xl + dx;

            yl = yr;
            yr = fun(xr);
        }

        // Linear interpolation
        let xg = (xl * yr - xr * yl) / (yr - yl);

        // Newton-Raphson method
        let (mut xn, mut nr_fail) = newton_raphson_p_alpha(&fun, tol, xg, (xl, xr), z, k);

        // Bisection (if NR fails)
        if nr_fail {
            (xn, nr_fail) = bisection(&fun, tol, (xl, xr), (yl, yr));
        }

        // Store root if valid
        if !nr_fail && xi <= xn && xn <= xf {
            r += xn;
        }

        // Step forward
        xl = xn + dx / 1000.0;
        xr = xl + dx;
    }
    r
}

// Newton-Raphson for P(alpha)
pub fn newton_raphson_p_alpha<F>(fun: F, tol: f64, mut xg: f64, (xl, xr): (f64, f64), z: &[f64], k: &[f64]) -> (f64, bool) where F: Fn(f64) -> f64 {
    // Initializing
    let mut nr_fail = false;
    let mut i = 0;
    let mut check = 1.0;
    let mut root = xg;

    // Looping until tolerance is met
    while tol < check {

        // Update iteration
        i += 1;

        // Evaluate function at xg
        let y = fun(xg);

        // Evaluate derivative at xg
        let der = p_alpha_der(xg, z, k);

        // Approximate root
        root = xg - y / der;

        // Method check
        if root < xl || root > xr || i > 10 || der == 0.0 {
            nr_fail = true;
            break;
        }

        // Compute error on root
        if root.abs() <= 1.0 {
            check = (xg - root).abs();
        } else {
            check = (1.0 - xg / root).abs();
        }
        
        // Update guess
        xg = root;
    }
    (root, nr_fail)
}

// P(alpha) function
pub fn p_alpha(alpha: f64, z: &[f64], k: &[f64]) -> f64 {
    let n = z.len();
    let mut sum = 0.0;

    for i in 0..n {
        sum += z[i] * (k[i] - 1.0) / (1.0 + alpha * (k[i] - 1.0));
    }
    sum
}

// P(alpha) derivative
pub fn p_alpha_der(alpha: f64, z: &[f64], k: &[f64]) -> f64 {
    let n = z.len();
    let mut sum = 0.0;

    for i in 0..n {
        let denominator = 1.0 + alpha * (k[i] - 1.0);
        sum += z[i] * (k[i] - 1.0).powi(2) / denominator.powi(2);
    }
    sum  *= -1.0;
    sum

}

// Root solver
pub fn root_solver<F>(fun: F, xi: f64, xf: f64, dx: f64, tol: f64) -> Vec<f64> where F: Fn(f64) -> f64 {
    // Offset initial bound
    let mut xl = xi - dx / 2.0;
    let mut xr = xl + dx;

    // Offset final bound
    let xft = xf + dx / 2.0;

    // Initialize storage
    let mut roots = Vec::new();

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
        let (mut xn, mut nr_fail) = newton_raphson(&fun, tol, xg, (xl, xr));

        // Bisection (if NR fails)
        if nr_fail {
            (xn, nr_fail) = bisection(&fun, tol, (xl, xr), (yl, yr));
        }

        // Store root if valid
        if !nr_fail && xi <= xn && xn <= xf {
            roots.push(xn);
        }

        // Step forward
        xl = xn + dx / 1000.0;
        xr = xl + dx;
    }
    roots
}

// Newton-Raphson
pub fn newton_raphson<F>(fun: F, tol: f64, mut xg: f64, (xl, xr): (f64, f64)) -> (f64, bool) where F: Fn(f64) -> f64 {
    // Initializing
    let mut nr_fail = false;
    let mut i = 0;
    let mut check = 1.0;
    let mut root = xg;
    let h = 1e-6;

    // Looping until tolerance is met
    while tol < check {

        // Update iteration
        i += 1;

        // Evaluate function at xg
        let y = fun(xg);

        // Evaluate derivative at xg
        let der = five_pt_stencil(&fun, xg, h);

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

// Bisection function
pub fn bisection<F>(fun: F, tol: f64, (mut xl, mut xr): (f64, f64), (mut yl, mut yr): (f64, f64)) -> (f64, bool) where F: Fn(f64) -> f64 {
    // Initializing
    let mut check  = 1.0;
    let mut singularity = false;
    let mut dy = 1e12;

    // Method
    while tol < check {
        let xmid = (xr + xl) / 2.0;
        let ymid = fun(xmid);

        if (yl * ymid) > 0.0 {
            xl = xmid;
            yl = ymid;
        } else {
            xr = xmid;
            yr = ymid;
        }

        // Singularity check
        let dy_new = (yl - yr).abs();
        if dy_new > dy {
            singularity = true;
            break;
        }
        dy = dy_new;

        // Tolerance check
        if (xr).abs() <= 1.0 {
            check = (xr - xl).abs()
        } else {
            check = (1.0 - xl / xr).abs()
        }
    }
    (xr, singularity)
}

// Five-point stencil
pub fn five_pt_stencil<F>(fun: &F, x: f64, h: f64) -> f64 where F: Fn(f64) -> f64 {
    let f_m2h = fun(x - 2.0*h);
    let f_m1h = fun(x - h);
    let f_p1h = fun(x + h);
    let f_p2h = fun(x + 2.0*h);

    (-f_p2h + 8.0*f_p1h - 8.0*f_m1h + f_m2h) / (12.0*h)
}

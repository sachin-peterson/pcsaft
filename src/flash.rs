use crate::structures::*;
use crate::solver::*;
use crate::fugacity::*;
use crate::constants::*;

pub fn flash_pcsaft(model: &PCSAFTModel, T: f64, P: f64, z: &[f64], k: &mut [f64], n: usize, tol: f64) -> (Vec<f64>, Vec<f64>, f64) {

    // Initializing
    let smax = 2;

    let mut alpha = 0.0;
    let mut x = vec![0.0; n];
    let mut y = vec![0.0; n];
    let mut theta_old = 0.0;

    for s in 0..smax {
        let mut p0 = 0.0;
        let mut p1 = 0.0;

        for i in 0..n {
            p0 += z[i] * (k[i] - 1.0);
            p1 += z[i] * (1.0 - 1.0 / k[i]);
        }
    
        // Number of phases
        let mut n_phase = 1;

        // Case I, II, or III
        if p0 * p1 <= 0.0 {
            let fun = |alpha: f64| p_alpha(alpha, &z, &k);
            alpha = root_solver_p_alpha(fun, 0.0, 1.0, 0.01, tol, &z, &k);            
            n_phase = 2;
        } else if p0 <= 0.0 {
            alpha = 0.0;
        } else if p1 >= 0.0 {
            alpha = 1.0;
        } else {
            panic!("System will not converge!")
        }

        // Calculating xi and yi
        for i in 0..n {
            x[i] = z[i] / (1.0 + alpha * (k[i] - 1.0));
            y[i] = k[i] * x[i];
        }

        // Normalizing
        let sum_x = x.iter().fold(0f64, |sum, i| sum + (*i as f64));
        let sum_y = y.iter().fold(0f64, |sum, i| sum + (*i as f64));
        for i in 0..n {
            x[i] /= sum_x;
            y[i] /= sum_y;
        }

        // Calculating densities
        let (rho_l, _) = density_solver(&x, T, P, model);
        let (_, rho_v) = density_solver(&y, T, P, model);

        // Compressibility factors
        let Zv = P / (rho_v * R * T);
        let Zl = P / (rho_l * R * T);

        // Calculating fugacities
        let mut fl = vec![0.0; n];
        let mut fv = vec![0.0; n];
        for i in 0..n {
            fl[i] = x[i] * P * fugacity_coeff(i, &x, Zl, T, rho_l, model);
            fv[i] = y[i] * P * fugacity_coeff(i, &y, Zv, T, rho_v, model);
        }

        println!("--- Iteration {} ---", s);
        println!("T:     {}", T);
        println!("z:     {:?}", z);
        println!("k:     {:?}", k);
        println!("alpha: {:.6}", alpha);
        println!("x:     {:?}", x);
        println!("y:     {:?}", y);
        println!("rho_l: {:.6} mol/m³", rho_l);
        println!("rho_v: {:.6} mol/m³", rho_v);
        println!("zl:    {}", Zl);
        println!("zv:    {}", Zv);
        println!("fl:    {:?}", fl);
        println!("fv:    {:?}", fv);

        // Calculating theta
        let mut theta: f64 = 0.0;
        for i in 0..n {
            theta += y[i] * f64::ln(fv[i]/fl[i])
        }

        // Check flash convergence
        if theta.abs() <= 1e-7 {
            break;
        } else if s != 0 && (theta - theta_old).abs() <= 1e-13 {
            if theta > 0.0 {
                println!("Feed is liquid-like, check if alpha = 0.");
                break;
            } else {
                println!("Feed is vapor-like, check if alpha = 1.");
                break;
            }
        } else {
            for i in 0..n {
                k[i] *= fl[i] / fv[i];
            }
        }

        // Update K-factors (1 phase present)
        if n_phase != 2 {
            for i in 0..n {
                k[i] *= theta.exp();
            }
        }

        // Check number of iterations
        if s >= smax {
            println!("Maximum number of iterations reached!")
        }

        // Update theta_old
        theta_old = theta;
    }

    (x, y, alpha)
}



// CORRECT VALUES
// 
// z:     [0.001, 0.999]
// k:     [576.4479837058954, 0.962353435956696]
// alpha: 0.024827
// x:     [6.54170357401373e-5, 0.9999345829642599]
// y:     [0.03770951835241864, 0.9622904816475814]
// rho_l: 3391.708075 mol/m³
// rho_v: 57.327034 mol/m³
// zl: 0.015518100054578654
// zv: 0.9181159661909821
// fl:    [1584.8677440112124, 191662.9374656585]
// fv:    [8484.95680307676, 184741.43875682115]

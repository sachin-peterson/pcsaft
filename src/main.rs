#![allow(non_snake_case)]
#![allow(dead_code)]

mod constants;
mod structures;

mod utils;
mod solver;
mod newton;

mod a_hc;
mod a_disp;
mod a_assoc;

mod model;
mod mixing;
mod references;

mod fugacity;
mod flash;

fn main() {

    // Components
    let _water = structures::PCSAFTComponent {
        mw: 18.015,
        m: 1.065,
        sigma: 3.0007,
        epsilon: 366.51,
        bondvol: 0.034868,
        epsilon_assoc: 2500.7,
        sites: vec!["H".into(), "O".into()]
    };

    let methane = structures::PCSAFTComponent {
        mw: 16.043,
        m: 1.0,
        sigma: 3.7039,
        epsilon: 150.03,
        bondvol: 0.0,
        epsilon_assoc: 0.0,
        sites: vec![]
    };

    let decane = structures::PCSAFTComponent {
        mw: 142.29,
        m: 4.6627,
        sigma: 3.8384,
        epsilon: 243.87,
        bondvol: 0.0,
        epsilon_assoc: 0.0,
        sites: vec![]
    };

    let components = vec![methane, decane];
    let n = components.len();

    // Association parameters
    let atol = 1e-8;
    let rtol = 1e-8;
    let alpha = 0.5;
    let iters = 50;
    let mixing = "Elliott";

    // Building the model
    let model = model::build_model(components, alpha, atol, rtol, iters, mixing);
    

    // TESTING TESTING TESTING    
    let T: f64 = 477.590;
    let P: f64 = 209000.0;
    
    // Methane
    let tc1 = 190.6;
    let pc1 = 4.599e6;
    let w1 = 0.012;

    // n-Decane
    let tc2 = 617.7;
    let pc2 = 2.110e6;
    let w2 = 0.492;

    // Vectorizing
    let tc: Vec<f64> = vec![tc1, tc2];
    let pc: Vec<f64> = vec![pc1, pc2];
    let w: Vec<f64> = vec![w1, w2];

    // Initial K's
    let mut k: Vec<f64> = vec![0.0; n];
    for i in 0..n {
        let lnk = (pc[i] / P).ln() + 5.37 * (1.0 + w[i]) * (1.0 - tc[i] / T);
        k[i] = lnk.exp();
    }

    // Feed fraction vector
    let z1 = 0.001;
    let z = vec![z1, 1.0 - z1];

    // Looking at residual function
    // let eta_vec: Vec<f64> = (0..1000)
    // .map(|i| 1e-6 + i as f64 * (0.49 - 1e-6) / 999.0)
    // .collect();

    // for eta in eta_vec {
    //     let res = fugacity::residual(&model, T, P, &x, n, eta);
    //     println!("res: {:?}", res)
    // }

    // TESTING FLASH
    flash::flash_pcsaft(&model, T, P, &z, &mut k, n, 1e-8);

}

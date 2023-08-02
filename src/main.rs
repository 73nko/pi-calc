use rand::Rng;
use rug::{ops::Pow, Float};
use std::f64::consts::SQRT_2;

const PRECISION: u32 = 1_000_000;

fn compute_pi_gauss_legendre() -> Float {
    let mut a: Float = Float::with_val(PRECISION, 1.0);
    let mut b: Float = Float::with_val(PRECISION, 1.0 / SQRT_2);
    let mut t: Float = Float::with_val(PRECISION, 0.25);
    let mut p: Float = Float::with_val(PRECISION, 1.0);

    for _ in 0..PRECISION {
        let a_next: Float = Float::with_val(PRECISION, Float::with_val(PRECISION, &a + &b) / 2.0);
        b = Float::with_val(PRECISION, Float::with_val(PRECISION, &a * &b).sqrt());
        t -= p.clone() * Float::with_val(PRECISION, &a - &a_next).pow(2);
        a = a_next;
        p *= 2.0;
    }

    Float::with_val(PRECISION, &a + &b).pow(2) / (4.0 * t)
}

fn compute_pi_montecarlo() -> f64 {
    let mut rng = rand::thread_rng();
    let mut hits = 0;

    for _ in 0..PRECISION {
        let x: f64 = rng.gen();
        let y: f64 = rng.gen();
        if x * x + y * y <= 1.0 {
            hits += 1;
        }
    }

    4.0 * (hits as f64) / (PRECISION as f64)
}

fn main() {
    let pi = compute_pi_gauss_legendre();
    println!("Pi is approximately using Gauss {}", pi);

    let pi = compute_pi_montecarlo();
    println!("Pi is approximately using Monte Carlo {}", pi);
}

#![allow(unused)]
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use sliding_blocks::{examples, solve_puzzle};

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Examples");
    group.sample_size(5);
    for (name, puzzle) in [
        ("ROYAL_ESCAPE", examples::ROYAL_ESCAPE),
        ("DIABOLICAL_BOX", examples::DIABOLICAL_BOX),
        ("GARBAGE_DISPOSAL", examples::GARBAGE_DISPOSAL),
    ] {
        group.bench_with_input(BenchmarkId::new(name, name), &puzzle, |b, p| {
            b.iter(|| solve_puzzle(p.start, p.goal))
        });
        /*group.bench_with_input(BenchmarkId::new("Iterative", i), i,
            |b, i| b.iter(|| fibonacci_fast(*i)));
        */
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

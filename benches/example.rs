#![allow(unused)]
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, SamplingMode};
use sliding_blocks::{examples, solve_puzzle};

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Examples");
    group.sample_size(10);
    group.sampling_mode(SamplingMode::Flat);
    for puzzle in examples::ALL_EXAMPLES {
        group.bench_with_input(BenchmarkId::new("lib_bfs", puzzle.name), &puzzle, |b, p| {
            b.iter(|| solve_puzzle(p.start, p.goal))
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

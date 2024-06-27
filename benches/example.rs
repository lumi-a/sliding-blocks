#![allow(unused)]
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use sliding_blocks::{examples, solve_puzzle_astar, solve_puzzle_lib_bfs};

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Examples");
    group.sample_size(10);
    for (name, puzzle) in [
        ("ROYAL_ESCAPE", examples::ROYAL_ESCAPE),
        ("DIABOLICAL_BOX", examples::DIABOLICAL_BOX),
        ("GARBAGE_DISPOSAL", examples::GARBAGE_DISPOSAL),
    ] {
        group.bench_with_input(BenchmarkId::new("lib_bfs", name), &puzzle, |b, p| {
            b.iter(|| solve_puzzle_lib_bfs(p.start, p.goal))
        });
        group.bench_with_input(BenchmarkId::new("astar", name), &puzzle, |b, p| {
            b.iter(|| solve_puzzle_astar(p.start, p.goal))
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

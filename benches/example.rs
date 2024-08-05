#![allow(unused)]
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, SamplingMode};
use sliding_blocks::{examples, solve_puzzle};

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Examples");
    group.sample_size(10);
    group.sampling_mode(SamplingMode::Flat);
    group.bench_function("Multi-Goal", |b| {
        b.iter(|| {
            examples::MULTI_GOAL_EXAMPLES
                .iter()
                .filter(|p| p.name != "🎩 The Time Machine")
                .for_each(|p| {
                    solve_puzzle(p.start, p.goal);
                })
        })
    });
    group.bench_function("Single-Goal", |b| {
        b.iter(|| {
            examples::SINGLE_GOAL_EXAMPLES.iter().for_each(|p| {
                solve_puzzle(p.start, p.goal);
            })
        })
    });
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

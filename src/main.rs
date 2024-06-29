use sliding_blocks::{examples, solve_puzzle};

fn main() {
    for puzzle in examples::ALL_EXAMPLES {
        print!("{} ", puzzle.name);
        solve_puzzle(puzzle.start, puzzle.goal);
    }
}

use sliding_blocks::{examples, solve_puzzle};

fn main() {
    for puzzle in examples::ALL_EXAMPLES {
        let computed = solve_puzzle(puzzle.start, puzzle.goal).len() - 1;
        let solution = puzzle.min_moves;
        println!(
            "{} {}{}{}",
            puzzle.name,
            computed,
            if solution == computed { "==" } else { "!=" },
            solution
        );
    }
}

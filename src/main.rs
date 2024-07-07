use sliding_blocks::{examples, solve_puzzle};

fn main() {
    for puzzle in examples::ALL_EXAMPLES {
        let computed = solve_puzzle(puzzle.start, puzzle.goal);
        let solution = puzzle.moves;
        println!(
            "{} {}{}{}",
            puzzle.name,
            computed,
            if solution == computed { "==" } else { "!=" },
            solution
        );
    }
}

use sliding_blocks::{examples, solve_puzzle_minmoves};

fn main() {
    for puzzle in examples::ALL_EXAMPLES {
        let computed = solve_puzzle_minmoves(puzzle.start, puzzle.goal)
            .unwrap()
            .unwrap_or(0);
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

use sliding_blocks::{examples, solve_puzzle_minmoves};

fn main() {
    for puzzle in examples::MULTI_GOAL_EXAMPLES
        .iter()
        .filter(|p| p.name != "🎩 The Time Machine")
    {
        let computed = solve_puzzle_minmoves(puzzle.start, puzzle.goal)
            .unwrap()
            .unwrap_or(0);
        let solution = puzzle.min_moves;
        println!(
            "{} {} {} {}",
            if solution == computed { "=" } else { "!=" },
            computed,
            solution,
            puzzle.name,
        );
    }
    println!("---");
    for puzzle in examples::SINGLE_GOAL_EXAMPLES {
        let computed = solve_puzzle_minmoves(puzzle.start, puzzle.goal)
            .unwrap()
            .unwrap_or(0);
        let solution = puzzle.min_moves;
        println!(
            "{} {} {} {}",
            if solution == computed { "=" } else { "!=" },
            computed,
            solution,
            puzzle.name,
        );
    }
}

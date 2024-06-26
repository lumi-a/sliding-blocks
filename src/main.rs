use sliding_blocks::{examples, solve_puzzle_lib_bfs, solve_puzzle_own_bfs};

fn main() {
    solve_puzzle_own_bfs(examples::ROYAL_ESCAPE.start, examples::ROYAL_ESCAPE.goal);
    solve_puzzle_lib_bfs(examples::ROYAL_ESCAPE.start, examples::ROYAL_ESCAPE.goal);
    solve_puzzle_own_bfs(
        examples::DIABOLICAL_BOX.start,
        examples::DIABOLICAL_BOX.goal,
    );
    solve_puzzle_lib_bfs(
        examples::DIABOLICAL_BOX.start,
        examples::DIABOLICAL_BOX.goal,
    );
    solve_puzzle_own_bfs(
        examples::GARBAGE_DISPOSAL.start,
        examples::GARBAGE_DISPOSAL.goal,
    );
    solve_puzzle_lib_bfs(
        examples::GARBAGE_DISPOSAL.start,
        examples::GARBAGE_DISPOSAL.goal,
    );
}

use super::{examples, solve_puzzle_minmoves};

#[test]
fn test_examples() {
    for puzzle in examples::ALL_EXAMPLES {
        assert_eq!(
            puzzle.min_moves,
            solve_puzzle_minmoves(&puzzle.start, &puzzle.goal)
        );
    }
}

use sliding_blocks::examples::ALL_EXAMPLES;
use sliding_blocks::{solve_puzzle, solve_puzzle_minmoves, SolvePuzzleError};

#[test]
fn test_examples() {
    for puzzle in ALL_EXAMPLES {
        assert_eq!(
            Ok(Some(puzzle.min_moves)),
            solve_puzzle_minmoves(puzzle.start, puzzle.goal)
        );
    }
}

#[test]
fn test_edgecases() {
    // Puzzle with exactly one block, which is a non-goalblock
    assert_eq!(
        solve_puzzle(".a.", "..."),
        Ok(Some(vec![".a.".to_string()]))
    );
    // Puzzle with exactly two blocks, neither of which are goalblocks
    assert_eq!(
        solve_puzzle(".ab", "..."),
        Ok(Some(vec![".ab".to_string()]))
    );
    // Feasible puzzle with exactly one block, which is a goal block
    assert_eq!(
        solve_puzzle(
            "
            a .....
            . .   .
            ... . .
                ...
        ",
            "
            . .....
            . .   .
            ... a .
                ...
        "
        ),
        Ok(Some(vec![
            "a .....
. .   .
... . .
    ..."
            .to_string(),
            ". .....
. .   .
... a .
    ..."
            .to_string()
        ]))
    );
    // Infeasible puzzle with exactly one block, which is a goal block
    assert_eq!(solve_puzzle("a .", ". a"), Ok(None));
    // Trivial puzzle with exactly one block, which is a goal block
    assert_eq!(solve_puzzle("a.", "a."), Ok(Some(vec!["a.".to_string()])));
    // Trivial puzzle with exactly two blocks, both being goalblocks
    assert_eq!(
        solve_puzzle("a.b", "a.b"),
        Ok(Some(vec!["a.b".to_string()]))
    );
    // Trivial puzzle with exactly two blocks, exactly one of which being goalblocks
    assert_eq!(
        solve_puzzle("a.b", "a.."),
        Ok(Some(vec!["a.b".to_string()]))
    );
    // Puzzle without any blocks
    assert_eq!(
        solve_puzzle("...", "..."),
        Ok(Some(vec!["...".to_string()]))
    );
    assert_eq!(
        solve_puzzle("...", ".."),
        Err(SolvePuzzleError::MismatchedBounds)
    );
    // Impossible puzzle
    assert_eq!(solve_puzzle("a.b", "..a"), Ok(None));
    // Simple one-move puzzle
    assert_eq!(
        solve_puzzle("a.", ".a"),
        Ok(Some(vec!["a.".to_string(), ".a".to_string()]))
    );
    // Simple one-move puzzle with redundant nongoal-block
    assert_eq!(
        solve_puzzle("a. b", ".a ."),
        Ok(Some(vec!["a. b".to_string(), ".a b".to_string()]))
    );
    // Simple one-move puzzle with redundant goal-block
    assert_eq!(
        solve_puzzle("a. b", ".a b"),
        Ok(Some(vec!["a. b".to_string(), ".a b".to_string()]))
    );
    // Puzzle with multiple goalblocks
    assert_eq!(
        solve_puzzle("a.b.", "..ab"),
        Ok(Some(vec![
            "a.b.".to_string(),
            "a..b".to_string(),
            "..ab".to_string()
        ]))
    );
    // Puzzle with one goalblock
    assert_eq!(
        solve_puzzle(
            "
            cb...
              a
        ",
            "
            a....
              .
        "
        ),
        Ok(Some(vec![
            "cb...
  a  "
            .to_string(),
            "c...b
  a  "
            .to_string(),
            "...cb
  a  "
            .to_string(),
            "a..cb
  .  "
            .to_string()
        ]))
    );
}

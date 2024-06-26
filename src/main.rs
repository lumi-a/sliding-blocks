use ::sliding_blocks::solve_puzzle;

fn main() {
    let puzzle = (
        "
      tt
      tt
    ......
    .ppoo.
     ypog
     yygg
      bb
      ..
    ",
        "
      ..
      ..
    ......
    ......
     ....
     ....
      tt
      tt
    ",
    );
    solve_puzzle(puzzle.0, puzzle.1);
}

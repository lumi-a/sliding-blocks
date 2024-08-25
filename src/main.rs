use sliding_blocks::solve_puzzle;

fn main() {
    let start = "
        A.b
         .
    ";
    let end = "
        ..A
         .
    ";
    println!("{:#?}", solve_puzzle(start, end));
}

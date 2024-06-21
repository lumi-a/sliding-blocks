use std::collections::BTreeSet;
use std::collections::HashMap;

#[derive(PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
struct Coor(u8, u8);

type CoorSet = BTreeSet<Coor>;

fn build_charmap(s: &str) -> HashMap<char, CoorSet> {
    s.lines()
        .enumerate()
        .flat_map(|(y, l)| {
            l.chars().enumerate().filter_map(move |(x, c)| {
                if !c.is_whitespace() {
                    Some((c, Coor(x as u8, y as u8)))
                } else {
                    None
                }
            })
        })
        .fold(HashMap::new(), |mut charmap, (c, coor)| {
            charmap.entry(c).or_insert_with(CoorSet::new).insert(coor);
            charmap
        })
}

pub fn solve_puzzle(start: &str, end: &str) {
    println!("{:?}", build_charmap(start));
}

fn main() {
    let puzzle = (
        "


    #t#CCvv
    A.#DFFv

     #.#DDrr
    #.BEEuu
    #s#E
    ",
        "
    #.#D
    B.#DD
    #.#CC
    #.AEE
    #.#EFF
    ",
    );
    solve_puzzle(puzzle.0, puzzle.1);
}

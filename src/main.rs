use std::collections::BTreeSet;
use std::collections::HashMap;

#[derive(PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
struct Coor(u8, u8);

type CoorSet = BTreeSet<Coor>;

pub fn solve_puzzle(start: &str, end: &str) {
    fn build_charmap(s: &str) -> HashMap<char, CoorSet> {
        let enumerated: Vec<(char, Coor)> = s
            .lines()
            .enumerate()
            .flat_map(|(y, l)| {
                l.chars()
                    .enumerate()
                    .filter(|x| !x.1.is_whitespace())
                    .map(move |(x, c)| (c, Coor(x as u8, y as u8)))
            })
            .collect();

        let mut charmap: HashMap<char, CoorSet> = HashMap::new();
        for (c, coor) in enumerated {
            charmap.entry(c).or_insert_with(CoorSet::new).insert(coor);
        }
        charmap
    }

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

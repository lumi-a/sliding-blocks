use std::cmp::min;
use std::collections::BTreeSet;
use std::collections::HashMap;

type CoorComponent = u8;
#[derive(PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
struct Coor(CoorComponent, CoorComponent);

type CoorSet = BTreeSet<Coor>;

const BOUNDS_CHAR: char = '.';

fn build_charmap(s: &str) -> HashMap<char, CoorSet> {
    let mut min_x = CoorComponent::MAX;
    let mut min_y = CoorComponent::MAX;

    let mut temp_coords = Vec::new();

    for (y, l) in s.lines().enumerate() {
        for (x, c) in l.chars().enumerate() {
            if !c.is_whitespace() {
                let x = x as CoorComponent;
                let y = y as CoorComponent;
                min_x = min(min_x, x);
                min_y = min(min_y, y);
                temp_coords.push((c, Coor(x, y)));
                if c != BOUNDS_CHAR {
                    temp_coords.push((BOUNDS_CHAR, Coor(x, y)));
                }
            }
        }
    }

    let shift = Coor(min_x, min_y);
    let mut charmap: HashMap<char, CoorSet> = HashMap::new();
    for (c, coor) in temp_coords {
        let shifted_coor = Coor(coor.0 - shift.0, coor.1 - shift.1);
        charmap
            .entry(c)
            .or_insert_with(CoorSet::new)
            .insert(shifted_coor);
    }

    charmap
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

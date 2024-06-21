use std::cmp::min;
use std::collections::BTreeSet;
use std::collections::HashMap;

type CoorComponent = u8;
#[derive(PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
struct Coor(CoorComponent, CoorComponent);

type CoorSet = BTreeSet<Coor>;

#[derive(PartialEq, Eq)]
struct Shape(BTreeSet<Coor>);
struct Offset(CoorComponent, CoorComponent);

const BOUNDS_CHAR: char = '.';

fn string_to_charcoorsmap(s: &str) -> HashMap<char, CoorSet> {
    // TODO: Change type of this to actually return a HashMap<char, (Shape, Offset)>, or something.
    // Shapekey extraction should be done later, otherwise we'd already have to return
    // a (Vec<Shape>, Vec<Set<Offset>>) right now, and that's just too abstract nonsense at this point.
    // Oh, of course then also rename the function, and variables that later call this function
    // (the current name is horrible)

    let mut min_x = CoorComponent::MAX;
    let mut min_y = CoorComponent::MAX;

    let mut temp_coords = Vec::new();

    for (y, l) in s.lines().enumerate() {
        for (x, c) in l.chars().enumerate() {
            if !c.is_whitespace() {
                let x = x as CoorComponent;
                let y = y as CoorComponent;
                // TODO: Check that x and y fit into the CoorComponent type.
                // Doing so would mean we'd have to return a Result instead.
                // Currently, this doesn't even panic, but continues innocently (yet wrongly)
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

pub fn solve_puzzle(start: &str, goal: &str) {
    let start_charcoorsmap = string_to_charcoorsmap(start);
    let goal_charcoorsmap = string_to_charcoorsmap(goal);

    // TODO: Handle this gracefully rather than panicking
    assert_eq!(
        start_charcoorsmap
            .get(&BOUNDS_CHAR)
            .unwrap_or(&CoorSet::new()),
        goal_charcoorsmap
            .get(&BOUNDS_CHAR)
            .unwrap_or(&CoorSet::new()),
        "The start and goal must have the same bounds."
    );
}

fn main() {
    let puzzle = (
        "
    #.#D
    B.#DD
    #.#CC
    #.AEE
    #.#EFF
    ",
        "
    #...
    D....
    DD...
    E....
    E.....
    ",
    );
    solve_puzzle(puzzle.0, puzzle.1);
}

use std::cmp::min;
use std::collections::BTreeSet;
use std::collections::HashMap;

type Coor = u8;
#[derive(PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
struct Coordinates(Coor, Coor);

type CoordinatesSet = BTreeSet<Coordinates>;

#[derive(PartialEq, Eq, Debug)]
struct Shape(BTreeSet<Coordinates>);
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Offset(Coor, Coor);

type Offsets = BTreeSet<Offset>;
struct Shapekey(Vec<Shape>);
struct Blockstate(Vec<Offsets>);

const BOUNDS_CHAR: char = '.';

fn string_to_charcoorsmap(s: &str) -> HashMap<char, CoordinatesSet> {
    // TODO: Rename the function, and variables that later call this function
    // (the current name is horrible)

    let mut min_x = Coor::MAX;
    let mut min_y = Coor::MAX;

    let mut temp_coords = Vec::new();

    for (y, l) in s.lines().enumerate() {
        for (x, c) in l.chars().enumerate() {
            if !c.is_whitespace() {
                let x = x as Coor;
                let y = y as Coor;
                // TODO: Check that x and y fit into the Coor type.
                // Doing so would mean we'd have to return a Result instead.
                // Currently, this doesn't even panic, but continues innocently (yet wrongly)
                // TODO: Then also let x, y be usize here, and only later convert
                // them to Coor after subtraction
                min_x = min(min_x, x);
                min_y = min(min_y, y);
                temp_coords.push((c, Coordinates(x, y)));
                if c != BOUNDS_CHAR {
                    temp_coords.push((BOUNDS_CHAR, Coordinates(x, y)));
                }
            }
        }
    }

    let shift = Coordinates(min_x, min_y);
    let mut charmap: HashMap<char, CoordinatesSet> = HashMap::new();
    for (c, coor) in temp_coords {
        let shifted_coor = Coordinates(coor.0 - shift.0, coor.1 - shift.1);
        charmap
            .entry(c)
            .or_insert_with(CoordinatesSet::new)
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
            .unwrap_or(&CoordinatesSet::new()),
        goal_charcoorsmap
            .get(&BOUNDS_CHAR)
            .unwrap_or(&CoordinatesSet::new()),
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

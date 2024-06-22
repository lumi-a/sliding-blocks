use std::cmp::{max, min};
use std::collections::BTreeSet;
use std::collections::HashMap;

// TODO: Perhaps it's better to abstract most of these into structs
type Coor = u8;
type Coordinates = (Coor, Coor);

type CoordinatesSet = BTreeSet<Coordinates>;
type CharToCoors = HashMap<char, CoordinatesSet>;

type Shape = BTreeSet<Coordinates>;
type Bounds = Shape;
type Offset = (Coor, Coor);

type Shapekey = Vec<Shape>;
type Offsets = Vec<Offset>;
type Blockstate = Vec<Offsets>; // TODO: Perhaps this is better done on the stack, e.g. with https://crates.io/crates/arrayvec

const BOUNDS_CHAR: char = '.';

fn string_to_chartocoors(s: &str) -> (CharToCoors, Coor, Coor) {
    // TODO: Rename the function, and variables that later call this function
    // (the current name is horrible)

    let mut min_x = Coor::MAX;
    let mut min_y = Coor::MAX;
    let mut max_x = Coor::MIN;
    let mut max_y = Coor::MIN;

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
                max_x = max(max_x, x);
                max_y = max(max_y, y);
                temp_coords.push((c, (x, y)));
                if c != BOUNDS_CHAR {
                    temp_coords.push((BOUNDS_CHAR, (x, y)));
                }
            }
        }
    }
    // TODO: Doesn't handle the case where the puzzle is empty
    let width = max_x - min_x + 1;
    let height = max_y - min_y + 1;

    let shift = (min_x, min_y);
    let mut chartocoors: CharToCoors = CharToCoors::new();
    for (c, coor) in temp_coords {
        let shifted_coor = (coor.0 - shift.0, coor.1 - shift.1);
        chartocoors
            .entry(c)
            .or_insert_with(CoordinatesSet::new)
            .insert(shifted_coor);
    }

    (chartocoors, width, height)
}

fn extract_shapekey(
    start_chartocoors: &CharToCoors,
    goal_chartocoors: &CharToCoors,
    width: Coor,
    height: Coor,
) -> (Bounds, Shapekey, Offsets) {
    let mut shapestooffsets: HashMap<Shape, Vec<Offset>> = HashMap::new();
    // TODO: Handle empty strings gracefully
    let bounds = (start_chartocoors.remove(&BOUNDS_CHAR).unwrap());

    for (c, start_coords) in start_chartocoors.iter() {
        // Extract min-x and min-y.
        // Assumes that start_coords is nonempty. TODO: Is that a misassumption?
        let mut min_x = Coor::MAX;
        let mut min_y = Coor::MAX;
        for coor in start_coords {
            min_x = min(min_x, coor.0);
            min_y = min(min_y, coor.1);
        }
        let shape = (start_coords
            .iter()
            .map(|coor| (coor.0 - min_x, coor.1 - min_y))
            .collect(),);
        shapestooffsets
            .entry(shape)
            .or_insert_with(Vec::new)
            .push((min_x, min_y));
    }
    let ughy: Vec<Shape> = shapestooffsets.keys().collect();
    // let shapekey = Shapekey(shapestooffsets.keys().collect());
    // let offsets = Offsets(shapestooffsets.values().collect());

    (bounds, shapekey, offsets)
}

fn print_puzzle(shapekey: &Shapekey, blockstate: &Blockstate, width: Coor, height: Coor) {}

pub fn solve_puzzle(start: &str, goal: &str) {
    let (start_chartocoors, width, height) = string_to_chartocoors(start);
    let (goal_chartocoors, goal_width, goal_height) = string_to_chartocoors(goal);

    // TODO: Handle this gracefully rather than panicking
    assert_eq!(
        start_chartocoors
            .get(&BOUNDS_CHAR)
            .unwrap_or(&CoordinatesSet::new()),
        goal_chartocoors
            .get(&BOUNDS_CHAR)
            .unwrap_or(&CoordinatesSet::new()),
        "The start and goal must have the same bounds."
    );

    // TODO: Handle this gracefully rather than panicking
    assert_eq!(width, goal_width, "start_width and goal_width don't match. This should never happen, as bounds are already asserted to be the same.");
    assert_eq!(height, goal_height, "start_height and goal_height don't match. This should never happen, as bounds are already asserted to be the same.");
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

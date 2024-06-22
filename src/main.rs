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

fn string_to_offsetshapes(s: &str) -> Vec<(Offset, Shape)> {
    let mut global_min_x = Coor::MAX;
    let mut global_min_y = Coor::MAX;

    let mut temp_coords: HashMap<char, (Coor, Coor, CoordinatesSet)> = HashMap::new();
    fn insert_temp_coords(
        c: char,
        x: Coor,
        y: Coor,
        temp_coords: &mut HashMap<char, (Coor, Coor, CoordinatesSet)>,
    ) {
        let existing =
            temp_coords
                .entry(c)
                .or_insert((Coor::MAX, Coor::MAX, CoordinatesSet::new()));
        existing.0 = min(existing.0, x);
        existing.1 = min(existing.1, y);
        existing.2.insert(Coordinates(x, y));
    }

    for (y, l) in s.lines().enumerate() {
        for (x, c) in l.chars().enumerate() {
            if !c.is_whitespace() {
                // TODO: Check that x and y fit into the CoorComponent type.
                // Doing so would mean we'd have to return a Result instead.
                // Currently, this doesn't even panic, but continues innocently (yet wrongly)
                // TODO: Then also let x, y be usize here, and only later convert
                // them to Coor after subtraction
                let x = x as Coor;
                let y = y as Coor;
                global_min_x = min(global_min_x, x);
                global_min_y = min(global_min_y, y);

                insert_temp_coords(c, x, y, &mut temp_coords);
                if c != BOUNDS_CHAR {
                    insert_temp_coords(BOUNDS_CHAR, x, y, &mut temp_coords);
                }
            }
        }
    }

    temp_coords
        .into_values()
        .map(|(min_x, min_y, coors)| {
            (
                Offset(min_x - global_min_x, min_y - global_min_y),
                Shape(
                    coors
                        .into_iter()
                        .map(|Coordinates(x, y)| Coordinates(x - min_x, y - min_y))
                        .collect(),
                ),
            )
        })
        .collect()
}

pub fn solve_puzzle(start: &str, goal: &str) {
    let start_offsetshapes = string_to_offsetshapes(start);
    let goal_offsetshapes = string_to_offsetshapes(goal);
    // TODO: Maybe check that start and goal have the same bounds?

    println!("{start_offsetshapes:?} {goal_offsetshapes:?}");
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

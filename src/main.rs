use colored::{self, Colorize};
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
type Offsets = BTreeSet<Offset>;
type Blockstate = Vec<Offsets>; // TODO: Perhaps this is better done on the stack, e.g. with https://crates.io/crates/arrayvec
type Coortable<T> = Vec<Vec<T>>;

fn intersect_coortables(a: &Coortable<bool>, b: &Coortable<bool>) -> Coortable<bool> {
    a.iter()
        .zip(b)
        .map(|(row_a, row_b)| {
            row_a.iter()
                .zip(row_b)
                .map(|(&elem_a, &elem_b)| elem_a && elem_b)
                .collect()
        })
        .collect()
}

// Nonintersectionkey[ShapeA][CoordinatesA][ShapeB][CoordinatesB] == true iff:
//   (ShapeA offset by CoordinatesA) ∩ (ShapeB offset by CoordinatesB) == ∅.
// To save on memory, we always assume that (ShapeA offset by CoordinatesA) is in bounds.
// This assumption will be satisfied when using the Nonintersectionkey in the algorithm.
type Nonintersectionkey = Vec<Coortable<Vec<Coortable<bool>>>>;

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
) -> (Bounds, Shapekey, Blockstate) {
    let mut shapestooffsets: HashMap<Shape, Offsets> = HashMap::new();
    // TODO: Handle empty strings gracefully
    // TODO: Also maybe don't use clone, but I'm not into ownership enough to think through how to handle this
    let bounds: Shape = start_chartocoors.get(&BOUNDS_CHAR).unwrap().clone();

    for (c, start_coords) in start_chartocoors.iter() {
        if c == &BOUNDS_CHAR {
            continue;
        }
        // Extract min-x and min-y.
        // Assumes that start_coords is nonempty. TODO: Is that a misassumption?
        let mut min_x = Coor::MAX;
        let mut min_y = Coor::MAX;
        for coor in start_coords {
            min_x = min(min_x, coor.0);
            min_y = min(min_y, coor.1);
        }
        let shape = start_coords
            .iter()
            .map(|coor| (coor.0 - min_x, coor.1 - min_y))
            .collect();
        shapestooffsets
            .entry(shape)
            .or_insert_with(BTreeSet::new)
            .insert((min_x, min_y));
    }
    let shapekey: Shapekey = shapestooffsets.keys().map(|shape| shape.clone()).collect();
    let blockstate: Blockstate = shapestooffsets
        .values()
        .map(|offsets| offsets.clone())
        .collect();

    (bounds, shapekey, blockstate)
}

fn build_nonintersectionkey(
    bounds: &Bounds,
    shapekey: &Shapekey,
    width: Coor,
    height: Coor,
) -> Nonintersectionkey {
    // Brace yourselves

    let mut nik = Nonintersectionkey::new();
    for shape_a in shapekey {
        let mut nik_a: Coortable<Vec<Coortable<bool>>> = Coortable::new();
        for ya in 0..height {
            let mut nik_ay: Vec<Vec<Coortable<bool>>> = Vec::new();
            for xa in 0..width {
                let mut nik_ayx: Vec<Coortable<bool>> = Vec::new();

                // TODO: Extract into shift-function
                let shifted_a: CoordinatesSet =
                    shape_a.iter().map(|(x, y)| (x + xa, y + ya)).collect();
                if shifted_a.is_subset(bounds) {
                    // Let the fun begin
                    for shape_b in shapekey {
                        let mut nik_ayx_b: Coortable<bool> = Coortable::new();
                        for yb in 0..height {
                            let mut nik_ayx_by: Vec<bool> = Vec::new();
                            for xb in 0..width {
                                // TODO: Extract into shift-function
                                let shifted_b: CoordinatesSet =
                                    shape_b.iter().map(|(x, y)| (x + xb, y + yb)).collect();

                                let nik_ayx_byx: bool = shifted_b.is_subset(bounds)
                                    && shifted_b.is_disjoint(&shifted_a);
                                nik_ayx_by.push(nik_ayx_byx);
                            }
                            nik_ayx_b.push(nik_ayx_by);
                        }
                        nik_ayx.push(nik_ayx_b);
                    }
                }
                nik_ay.push(nik_ayx);
            }
            nik_a.push(nik_ay);
        }
        nik.push(nik_a);
    }
    nik
}

fn get_neighboring_blockstates(
    blockstate: &Blockstate,
    nonintersectionkey: &Nonintersectionkey,
    width: Coor, height: Coor
) -> Vec<Blockstate> {
    // TODO: Create current nonintersections using dynamic programming:
    // In the end,
    //  `left_nonintersection[shape][i] ∩ right_nonintersection[shape][?-i]
    //   ∩ (nik[shape][x1][y1][shape]∩…∩(nik[shape][x(i-1)][y(i-1)][shape]))
    //   ∩ (nik[shape][x(i+1)][y(i+1)][shape]∩…∩(nik[shape][x?][y?][shape]))
    //  `
    // will describe exactly the positions that block `i` of shape `shape` is allowed to move to.
    // The latter ugly thing can also be implemented using dynamic programming.
    // TODO: Maybe this is faster using Bitvecs rather than Coortable<Bool>
    // TODO: Capacity can be calculated ahead of time.
    // TODO: This not only assumes non-empty bounds, but also that at least one block exists (or does it?)
    /*
    let all_coors : Coortable<bool> = vec![vec![true; height as usize]; width as usize];
    let mut left_nonintersection: Vec<Coortable<bool>> = vec![all_coors.clone()];
    let mut right_nonintersection: Vec<Coortable<bool>> = vec![all_coors.clone()];
    for (shape_ix, shape_offsets) in blockstate.iter().enumerate() {
        let nik_shape = &nonintersectionkey[shape_ix];
        for (x,y) in shape_offsets {
            let new_coortable = intersect_coortables(
                left_nonintersection.last().unwrap(), nik_elem
            );
            left_nonintersection.push(
                new_coortable
            )
        }
    }
    */
    
    let mut neighboring_blockstates: Vec<Blockstate> = Vec::new();
    for (shape_ix, shape_offsets) in blockstate.iter().enumerate() {

        for offset in shape_offsets.iter().enumerate() {
            let mutated_offsets = mutate(offset);
            for mutated_offset in mutated_offsets {
                let mutated_shape_offsets = // replace shape_offsets[ix] with mutated_offset
                neighboring_blockstates.push(mutated_shape_offsets)
            }
        }
    }
    neighboring_blockstates
}

fn print_puzzle(
    bounds: &Bounds,
    shapekey: &Shapekey,
    blockstate: &Blockstate,
    width: Coor,
    height: Coor,
) {
    // Create vec of blocks:
    let mut blocks: Vec<CoordinatesSet> = Vec::new();
    for (shape, offsets) in shapekey.iter().zip(blockstate.iter()) {
        for offset in offsets {
            let block: CoordinatesSet = shape
                .iter()
                .map(|coor| (coor.0 + offset.0, coor.1 + offset.1))
                .collect();
            blocks.push(block);
        }
    }
    blocks.sort(); // ensures consistent indices
    let blocks = blocks;

    const IN_BLOCK: &str = "██";
    const IN_BOUNDS: &str = "  ";
    const OUT_OF_BOUNDS: &str = "░░";
    println!("{}", OUT_OF_BOUNDS.repeat(width as usize + 2));
    for y in 0..height {
        print!("{}", OUT_OF_BOUNDS);
        for x in 0..width {
            // Find block that contains (x, y)
            let option_block_ix: Option<usize> =
                blocks.iter().position(|block| block.contains(&(x, y)));
            match option_block_ix {
                Some(block_ix) => {
                    let r = ((block_ix + 1) * 7573 % 256) as u8;
                    let g = ((block_ix + 1) * 6841 % 256) as u8;
                    let b = ((block_ix + 1) * 5953 % 256) as u8;
                    print!("{}", IN_BLOCK.truecolor(r, g, b));
                }
                None => {
                    if bounds.contains(&(x, y)) {
                        print!("{}", IN_BOUNDS);
                    } else {
                        print!("{}", OUT_OF_BOUNDS);
                    }
                }
            }
        }
        println!("{}", OUT_OF_BOUNDS);
    }
    println!("{}", OUT_OF_BOUNDS.repeat(width as usize + 2));
}

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

    println!("{},{}", width, height);

    let (bounds, shapekey, blockstate) =
        extract_shapekey(&start_chartocoors, &goal_chartocoors, width, height);

    // TODO: Remove this
    print_puzzle(&bounds, &shapekey, &blockstate, width, height);

    let nonintersectionkey = build_nonintersectionkey(&bounds, &shapekey, width, height);
}

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
      ..
      ..
    ",
    );
    solve_puzzle(puzzle.0, puzzle.1);
}

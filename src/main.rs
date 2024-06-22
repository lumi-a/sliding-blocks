use colored::{self, Colorize};
use std::cmp::{max, min};
use std::collections::hash_map::Entry::Occupied;
use std::collections::hash_map::Entry::Vacant;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::{BTreeSet, VecDeque};

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
            row_a
                .iter()
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
        for xa in 0..width {
            let mut nik_ax: Vec<Vec<Coortable<bool>>> = Vec::new();
            for ya in 0..height {
                let mut nik_axy: Vec<Coortable<bool>> = Vec::new();

                // TODO: Extract into shift-function
                let shifted_a: CoordinatesSet =
                    shape_a.iter().map(|(x, y)| (x + xa, y + ya)).collect();
                if shifted_a.is_subset(bounds) {
                    // Let the fun begin
                    for shape_b in shapekey {
                        let mut nik_axy_b: Coortable<bool> = Coortable::new();
                        for xb in 0..width {
                            let mut nik_axy_bx: Vec<bool> = Vec::new();
                            for yb in 0..height {
                                // TODO: Extract into shift-function
                                let shifted_b: CoordinatesSet =
                                    shape_b.iter().map(|(x, y)| (x + xb, y + yb)).collect();

                                let nik_axy_bxy: bool = shifted_b.is_subset(bounds)
                                    && shifted_b.is_disjoint(&shifted_a);
                                nik_axy_bx.push(nik_axy_bxy);
                            }
                            nik_axy_b.push(nik_axy_bx);
                        }
                        nik_axy.push(nik_axy_b);
                    }
                }
                nik_ax.push(nik_axy);
            }
            nik_a.push(nik_ax);
        }
        nik.push(nik_a);
    }
    nik
}

fn get_neighboring_blockstates(
    blockstate: &Blockstate,
    nonintersectionkey: &Nonintersectionkey,
    width: Coor,
    height: Coor,
) -> Vec<Blockstate> {
    // TODO: Create current nonintersections using dynamic programming:
    // In the end,
    //  `left_nonintersection[shape][i] ∩ right_nonintersection[shape][?-i]
    //   ∩ [(nik[shape][x1][y1][shape]∩…∩(nik[shape][x(i-1)][y(i-1)][shape]))
    //     ∩ (nik[shape][x(i+1)][y(i+1)][shape]∩…∩(nik[shape][x?][y?][shape]))
    //   ]`
    // will describe exactly the positions that block `i` of shape `shape` is allowed to move to.
    // The latter ugly thing can also be implemented using dynamic programming.
    // TODO: Maybe this is faster using Bitvecs rather than Coortable<Bool>
    // TODO: Capacity can be calculated ahead of time.
    // TODO: This not only assumes non-empty bounds, but also that at least TWO blocks exist!
    //  (otherwise we don't get to snack on the free bounds-check)
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
    */ // for now, just this bad implementation of bfs:
    // TODO: For this, the order of the shapekeys actually affects performance, so
    // we should make sure that larger shapes get sorted first
    let dfs = |moving_shape_ix: usize,
               moving_offset: Offset,
               trimmed_shape_offsets: &Offsets|
     -> Vec<Offset> {
        let is_legal = |offsetty: Offset| -> bool {
            for (shape_ix, shape_offsets) in blockstate.iter().enumerate() {
                if shape_ix == moving_shape_ix {
                    continue;
                }
                for offset in shape_offsets {
                    // TODO: How bad are these "as usize" conversions?
                    // If they're really bad, I might just end up using usize as the type for Coor
                    if !nonintersectionkey[shape_ix][offset.0 as usize][offset.1 as usize]
                        [moving_shape_ix][offsetty.0 as usize][offsetty.1 as usize]
                    {
                        return false;
                    }
                }
            }
            for offset in trimmed_shape_offsets {
                if !nonintersectionkey[moving_shape_ix][offset.0 as usize][offset.1 as usize]
                    [moving_shape_ix][offsetty.0 as usize][offsetty.1 as usize]
                {
                    return false;
                }
            }
            true
        };
        let mut legal_offsets: Vec<Offset> = Vec::new();
        let mut seen_offsets: BTreeSet<Offset> = BTreeSet::new();
        seen_offsets.insert(moving_offset);
        let mut stack: Vec<Offset> = vec![moving_offset];
        while !stack.is_empty() {
            let offset = stack.pop().unwrap();

            // TODO: Maybe this can be made faster by not going back in
            // the direction we just came from

            // TODO: This is AWFUL code.
            // And no, we can't just extract the array into match-arms and for-loop
            // over the match-result, because arrays of different length are
            // different types
            // Should we extend the intersection-vector to reach one more cell? Hmm
            let inserty = |new_offset: Offset| {
                if new_offset.0 < width
                    && new_offset.1 < height
                    && seen_offsets.insert(new_offset)
                    && is_legal(new_offset)
                {
                    legal_offsets.push(new_offset);
                    stack.push(new_offset);
                }
            };
            // TODO: I hate this
            if offset.0 > 0 {
                if offset.1 > 0 {
                    [
                        (offset.0, offset.1 + 1),
                        (offset.0 + 1, offset.1),
                        (offset.0, offset.1 - 1),
                        (offset.0 - 1, offset.1),
                    ]
                    .map(inserty);
                } else {
                    [
                        (offset.0, offset.1 + 1),
                        (offset.0 + 1, offset.1),
                        (offset.0 - 1, offset.1),
                    ]
                    .map(inserty);
                }
            } else {
                if offset.1 > 0 {
                    [
                        (offset.0, offset.1 + 1),
                        (offset.0 + 1, offset.1),
                        (offset.0, offset.1 - 1),
                    ]
                    .map(inserty);
                } else {
                    [(offset.0, offset.1 + 1), (offset.0 + 1, offset.1)].map(inserty);
                }
            }
        }

        legal_offsets
    };

    // It's okay to gather all these into a vector rather than a set,
    // because all neighbors WILL be unique.
    let mut neighboring_blockstates: Vec<Blockstate> = Vec::new();
    for (shape_ix, shape_offsets) in blockstate.iter().enumerate() {
        for offset in shape_offsets.iter() {
            let mut trimmed_shape_offsets = shape_offsets.clone();
            trimmed_shape_offsets.remove(offset);
            for mutated_offset in dfs(shape_ix, *offset, &trimmed_shape_offsets) {
                let mut mutated_shape_offsets = trimmed_shape_offsets.clone();
                mutated_shape_offsets.insert(mutated_offset);
                let mut new_blockstate = blockstate.clone();
                new_blockstate[shape_ix] = mutated_shape_offsets;
                neighboring_blockstates.push(new_blockstate);
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
    // TODO: Colors aren't the best choice here, because colors WILL change when
    //  blocks move, giving the illusion of some blocks having changed shapes

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

fn puzzle_bfs_with_path_reconstruction<F, G>(
    blockstate: &Blockstate,
    get_neighbors: F,
    goal_check: G,
) -> Option<Vec<Blockstate>>
where
    F: Fn(&Blockstate) -> Vec<Blockstate>,
    G: Fn(&Blockstate) -> bool,
{
    let mut queue: VecDeque<Blockstate> = VecDeque::new();
    // TODO: Better datastructures?
    // TODO: Definitely not memory efficient, we should store pointers instead
    let mut seen: HashMap<Blockstate, Option<Blockstate>> = HashMap::new();
    queue.push_back(blockstate.clone());
    seen.insert(blockstate.clone(), None);
    while !queue.is_empty() {
        let blockstate = queue.pop_front().unwrap();
        if goal_check(&blockstate) {
            let mut path: Vec<Blockstate> = vec![blockstate.clone()];
            let mut curr = blockstate.clone();
            while let Some(prev) = seen[&curr].clone() {
                path.push(prev.clone());
                curr = prev;
            }
            path.reverse();
            return Some(path);
        }
        for neighbor in get_neighbors(&blockstate) {
            match seen.entry(neighbor.clone()) {
                Occupied(_) => (),
                Vacant(entry) => {
                    entry.insert(Some(blockstate.clone()));
                    queue.push_back(neighbor);
                }
            }
        }
    }
    None
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

    let (bounds, shapekey, start_blockstate) =
        extract_shapekey(&start_chartocoors, &goal_chartocoors, width, height);

    // TODO: Remove this
    print_puzzle(&bounds, &shapekey, &start_blockstate, width, height);
    println!();
    println!();

    let nonintersectionkey = build_nonintersectionkey(&bounds, &shapekey, width, height);

    // TODO: Implement actual goal conditions
    // TODO: Implement A*
    // TODO: Add tests
    // TODO: delete this
    let mut big_ix = 0;
    let mut big_val = shapekey[0].len();
    for (ix, val) in shapekey.iter().enumerate() {
        let lenny = val.len();
        if lenny > big_val {
            big_ix = ix;
            big_val = lenny;
        }
    }

    let path = puzzle_bfs_with_path_reconstruction(
        &start_blockstate,
        |blockstate| get_neighboring_blockstates(&blockstate, &nonintersectionkey, width, height),
        |blockstate| blockstate[big_ix].iter().next().unwrap().1 > 5,
    )
    .unwrap();
    for blockstate in &path {
        print_puzzle(&bounds, &shapekey, &blockstate, width, height);
    }
    println!("{}", path.len());

    // TODO: Print path
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

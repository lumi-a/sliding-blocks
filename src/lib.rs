pub mod examples;
#[cfg(test)]
mod tests;

use bitvec::prelude::*;
use colored::{self, Colorize};
use itertools::Itertools;
use std::cmp::{max, min, Ordering};
use std::collections::{BTreeSet, HashMap};
use wasm_bindgen::prelude::*;

type Coor = u8;
type Width = Coor;
type Height = Coor;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Point(Coor, Coor);
// A "global" coordinate, in contrast to the Offset of a Shape
impl Point {
    // TODO: It'd be nice to have this using std::ops::Add, std::ops::Sub,
    // but those seem to kinda consume ownership? There are ways to avoid that,
    // I hear, but I don't know enough about ownership yet to understand those,
    // sorry.
    fn add(&self, other: &Point) -> Point {
        Point(self.0 + other.0, self.1 + other.1)
    }
    fn sub(&self, other: &Point) -> Point {
        Point(self.0 - other.0, self.1 - other.1)
    }
}
impl From<&Offset> for Point {
    fn from(offset: &Offset) -> Self {
        Self(offset.0, offset.1)
    }
}

type Points = BTreeSet<Point>;
type CharToPoints = HashMap<char, Points>;
type ReconstructionMap = HashMap<(Shape, Offset), char>;

type Shape = BTreeSet<Point>; // nonempty. min-x == 0, min-y == 0

// TODO: Make this a struct method?
fn get_extremes(coordinates_set: &Points) -> (Point, Point) {
    // Extract min-x and min-y.
    // Assumes that coordinatesSet is nonempty.
    // TODO: Is that a misassumption?
    let mut min_x = Coor::MAX;
    let mut min_y = Coor::MAX;
    let mut max_x = Coor::MIN;
    let mut max_y = Coor::MIN;
    for point in coordinates_set {
        max_x = min(max_x, point.0);
        max_y = min(max_y, point.1);
        min_x = min(min_x, point.0);
        min_y = min(min_y, point.1);
    }
    (Point(min_x, min_y), Point(max_x, max_y))
}
type Bounds = Shape; // min-x == 1, min-y == 1.

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Offset(Coor, Coor); // Should be >= (1, 1) for most offsets
impl Offset {
    #[inline]
    fn up(&self) -> Offset {
        Offset(self.0, self.1 + 1)
    }
    #[inline]
    fn down(&self) -> Offset {
        Offset(self.0, self.1 - 1)
    }
    #[inline]
    fn left(&self) -> Offset {
        Offset(self.0 - 1, self.1)
    }
    #[inline]
    fn right(&self) -> Offset {
        Offset(self.0 + 1, self.1)
    }
}
impl From<&Point> for Offset {
    fn from(point: &Point) -> Self {
        Self(point.0, point.1)
    }
}

type Offsets = BTreeSet<Offset>;
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Blockstate {
    nongoal_offsets: Vec<Offsets>,
    goal_offsets: Vec<Offset>,
}
type Shapekey = Vec<Shape>;
type GoalShapekeyKey = Vec<usize>; // Given an index in the Blockstate.goal_blocks vec, what is the index of its shape in the shapekey vec?
type GoalTargetOffsets = Vec<Offset>; // At what offset is a block in a goal-position?

// Nonintersectionkey[ShapeA][CoordinatesA][ShapeB][CoordinatesB] == true iff:
//   (ShapeA offset by CoordinatesA) ∩ (ShapeB offset by CoordinatesB) == ∅.
// To save on memory, we always assume that (ShapeA offset by CoordinatesA) is in bounds.
// This assumption will be satisfied when using the Nonintersectionkey in the algorithm.

type ShapevecForNik<T> = Vec<T>;
struct Nonintersectionkey(usize, ShapevecForNik<Vec<ShapevecForNik<BitVec>>>);
impl std::ops::Index<(usize, &Offset, usize, &Offset)> for Nonintersectionkey {
    type Output = bool;

    #[inline]
    fn index(
        &self,
        (shape_ix_a, Offset(xa, ya), shape_ix_b, Offset(xb, yb)): (usize, &Offset, usize, &Offset),
    ) -> &Self::Output {
        &self.1[shape_ix_a][*xa as usize + (*ya as usize) * self.0][shape_ix_b]
            [*xb as usize + (*yb as usize) * self.0]
    }
}

const BOUNDS_CHAR: char = '.';

fn build_nonintersectionkey(
    bounds: &Bounds,
    shapekey: &Shapekey,
    width: Width,
    height: Height,
) -> Nonintersectionkey {
    // TODO: This can be done faster by first storing the relative offsets for which the
    // nik is true, and then filling in the actual nik.

    // Brace yourselves

    let mut nikvec: ShapevecForNik<Vec<ShapevecForNik<BitVec>>> = Vec::new();
    for shape_a in shapekey {
        let mut nik_a: Vec<ShapevecForNik<BitVec>> = Vec::new();
        for ya in 0..=(height + 1) {
            for xa in 0..=(width + 1) {
                let mut nik_a_xy: ShapevecForNik<BitVec> = ShapevecForNik::new();
                let shift_a = Point(xa, ya);

                // TODO: Extract into shift-function
                let shifted_a: Points = shape_a.iter().map(|p| p.add(&shift_a)).collect();
                if shifted_a.is_subset(bounds) {
                    // Let the fun begin
                    for shape_b in shapekey {
                        let mut nik_a_xy_b: BitVec = BitVec::new();
                        for yb in 0..=(height + 1) {
                            for xb in 0..=(width + 1) {
                                let shift_b = Point(xb, yb);
                                // TODO: Extract into shift-function
                                let shifted_b: Points =
                                    shape_b.iter().map(|p| p.add(&shift_b)).collect();

                                let nik_a_xy_b_xy: bool = shifted_b.is_subset(bounds)
                                    && shifted_b.is_disjoint(&shifted_a);
                                nik_a_xy_b.push(nik_a_xy_b_xy);
                            }
                        }
                        nik_a_xy_b.shrink_to_fit();
                        nik_a_xy.push(nik_a_xy_b);
                    }
                }
                nik_a_xy.shrink_to_fit();
                nik_a.push(nik_a_xy);
            }
        }
        nik_a.shrink_to_fit();
        nikvec.push(nik_a);
    }
    nikvec.shrink_to_fit();
    Nonintersectionkey(width as usize + 2, nikvec)
}

fn get_neighboring_blockstates(
    blockstate: &Blockstate,
    nonintersectionkey: &Nonintersectionkey,
    goal_shapekey_key: &GoalShapekeyKey,
) -> Vec<Blockstate> {
    // TODO: I have no idea if `&dyn Fn(Offset) -> bool` is the right signature as I didn't learn about `&dyn` yet
    let dfs_general = |initial_offset: Offset, is_legal: &dyn Fn(&Offset) -> bool| {
        let mut legal_offsets: Vec<Offset> = Vec::new();
        let mut seen_offsets: BTreeSet<Offset> = BTreeSet::new(); // TODO: Different data structures?
        seen_offsets.insert(initial_offset.clone());

        // To eliminate backtracking:
        enum CameFrom {
            Up,
            Down,
            Left,
            Right,
        }

        let mut stack: Vec<(Offset, CameFrom)> = Vec::new();

        // Initial setup:
        for (new_offset, new_dir) in [
            (initial_offset.up(), CameFrom::Up),
            (initial_offset.down(), CameFrom::Down),
            (initial_offset.left(), CameFrom::Left),
            (initial_offset.right(), CameFrom::Right),
        ] {
            // No need to check if seen_offsets contains them, as we know it won't
            if is_legal(&new_offset) {
                seen_offsets.insert(new_offset.clone());
                legal_offsets.push(new_offset.clone());
                stack.push((new_offset, new_dir));
            }
        }

        while let Some((offset, dir)) = stack.pop() {
            for (new_offset, new_dir) in match dir {
                CameFrom::Up => [
                    (offset.up(), CameFrom::Up),
                    (offset.left(), CameFrom::Left),
                    (offset.right(), CameFrom::Right),
                ],
                CameFrom::Down => [
                    (offset.down(), CameFrom::Down),
                    (offset.left(), CameFrom::Left),
                    (offset.right(), CameFrom::Right),
                ],
                CameFrom::Left => [
                    (offset.up(), CameFrom::Up),
                    (offset.down(), CameFrom::Down),
                    (offset.left(), CameFrom::Left),
                ],
                CameFrom::Right => [
                    (offset.up(), CameFrom::Up),
                    (offset.down(), CameFrom::Down),
                    (offset.right(), CameFrom::Right),
                ],
            } {
                // TODO: Can we check seen_offsets membership without cloning?
                // TODO: Which of these checks is faster? And shouldn't we be using pathfinding::directed::dfs::dfs instead?
                if is_legal(&new_offset) && seen_offsets.insert(new_offset.clone()) {
                    legal_offsets.push(new_offset.clone());
                    stack.push((new_offset, new_dir));
                }
            }
        }

        legal_offsets
    };
    let dfs_nongoal = |movingshape_ix: usize,
                       moving_offset: Offset,
                       trimmed_movingshape_offsets: &Offsets|
     -> Vec<Offset> {
        let is_legal = |offsety: &Offset| -> bool {
            // TODO: This function assumes there are other blocks on the field, because
            // we currently use their nonintersectionkeys to additionally verify that
            // offsety is in-bounds
            for (shape_ix, shape_offsets) in blockstate.nongoal_offsets.iter().enumerate() {
                if shape_ix == movingshape_ix {
                    continue;
                }
                for offset in shape_offsets {
                    // TODO: How bad are these "as usize" conversions?
                    // If they're really bad, I might just end up using usize as the type for Coor
                    if !nonintersectionkey[(shape_ix, offset, movingshape_ix, offsety)] {
                        return false;
                    }
                }
            }
            for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
                let shape_ix = goal_shapekey_key[goalvec_ix];
                if !nonintersectionkey[(shape_ix, offset, movingshape_ix, offsety)] {
                    return false;
                }
            }
            for offset in trimmed_movingshape_offsets {
                if !nonintersectionkey[(movingshape_ix, offset, movingshape_ix, offsety)] {
                    return false;
                }
            }
            true
        };
        dfs_general(moving_offset, &is_legal)
    };

    let dfs_goal = |moving_goalvec_ix: usize, moving_offset: Offset| -> Vec<Offset> {
        let moving_shape_ix = goal_shapekey_key[moving_goalvec_ix];
        let is_legal = |offsety: &Offset| -> bool {
            // TODO: This function assumes there are other blocks on the field, because
            // we currently use their nonintersectionkeys to additionally verify that
            // offsety is in-bounds
            for (shape_ix, shape_offsets) in blockstate.nongoal_offsets.iter().enumerate() {
                for offset in shape_offsets {
                    // TODO: How bad are these "as usize" conversions?
                    // If they're really bad, I might just end up using usize as the type for Coor
                    if !nonintersectionkey[(shape_ix, offset, moving_shape_ix, offsety)] {
                        return false;
                    }
                }
            }
            for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
                if goalvec_ix == moving_goalvec_ix {
                    continue;
                }
                let shape_ix = goal_shapekey_key[goalvec_ix];
                if !nonintersectionkey[(shape_ix, offset, moving_shape_ix, offsety)] {
                    return false;
                }
            }
            true
        };
        dfs_general(moving_offset, &is_legal)
    };

    // It's okay to gather all these into a vector rather than a set,
    // because all neighbors WILL be unique.
    let mut neighboring_blockstates: Vec<Blockstate> = Vec::new();
    // Start with blockstate.goal_offsets first, to hopefully find the goal a little sooner
    for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
        // TODO: avoid .clone() here?
        for mutated_offset in dfs_goal(goalvec_ix, offset.clone()) {
            let mut new_blockstate = blockstate.clone();
            new_blockstate.goal_offsets[goalvec_ix] = mutated_offset;
            neighboring_blockstates.push(new_blockstate);
        }
    }
    for (shape_ix, offsets) in blockstate.nongoal_offsets.iter().enumerate() {
        for offset in offsets.iter() {
            let mut trimmed_shape_offsets = offsets.clone();
            trimmed_shape_offsets.remove(offset);
            // TODO: avoid .clone() here?
            for mutated_offset in dfs_nongoal(shape_ix, offset.clone(), &trimmed_shape_offsets) {
                let mut mutated_shape_offsets = trimmed_shape_offsets.clone();
                mutated_shape_offsets.insert(mutated_offset);
                let mut new_blockstate = blockstate.clone();
                new_blockstate.nongoal_offsets[shape_ix] = mutated_shape_offsets;
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
    goal_shapekey_key: &GoalShapekeyKey,
    width: Width,
    height: Height,
) {
    // TODO: Colors aren't the best choice here, because colors WILL change after
    // blocks move, giving the illusion of some blocks having changed shapes

    // Create vec of blocks:
    let mut blocks: Vec<Points> = Vec::new();
    for (shape, offsets) in shapekey.iter().zip(blockstate.nongoal_offsets.iter()) {
        for offset in offsets {
            // Just create shift-shape method already
            let block: Points = shape.iter().map(|p| p.add(&offset.into())).collect();
            blocks.push(block);
        }
    }
    for (shape_ix, offset) in goal_shapekey_key.iter().zip(blockstate.goal_offsets.iter()) {
        let block: Points = shapekey[*shape_ix]
            .iter()
            .map(|p| p.add(&offset.into()))
            .collect();
        blocks.push(block);
    }
    let blocks = blocks;

    const IN_BLOCK: &str = "██";
    const IN_BOUNDS: &str = "  ";
    const OUT_OF_BOUNDS: &str = "░░";
    for y in 0..=(height + 1) {
        for x in 0..=(width + 1) {
            // Find block that contains (x, y)
            let option_block_ix: Option<usize> =
                blocks.iter().position(|block| block.contains(&Point(x, y)));
            match option_block_ix {
                Some(block_ix) => {
                    let r = ((block_ix + 1) * 7573 % 256) as u8;
                    let g = ((block_ix + 1) * 6841 % 256) as u8;
                    let b = ((block_ix + 1) * 5953 % 256) as u8;
                    print!("{}", IN_BLOCK.truecolor(r, g, b));
                }
                None => {
                    if bounds.contains(&Point(x, y)) {
                        print!("{}", IN_BOUNDS);
                    } else {
                        print!("{}", OUT_OF_BOUNDS);
                    }
                }
            }
        }
        println!();
    }
}

fn puzzle_preprocessing(
    start: &str,
    goal: &str,
) -> (
    Bounds,
    Shapekey,
    Blockstate,
    Nonintersectionkey,
    GoalShapekeyKey,
    GoalTargetOffsets,
    ReconstructionMap,
    Width,
    Height,
) {
    fn string_to_chartocoors(s: &str) -> (CharToPoints, Width, Height) {
        let mut min_x = Coor::MAX;
        let mut min_y = Coor::MAX;
        let mut max_x = Coor::MIN;
        let mut max_y = Coor::MIN;

        let mut temp_points: Vec<(char, Point)> = Vec::new();

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
                    // TODO: Should min_x, min_y, max_x, max_y be handles using get_extremes instead?
                    //       Would make code less performant, but more legible
                    min_x = min(min_x, x);
                    min_y = min(min_y, y);
                    max_x = max(max_x, x);
                    max_y = max(max_y, y);
                    temp_points.push((c, Point(x, y)));
                    if c != BOUNDS_CHAR {
                        temp_points.push((BOUNDS_CHAR, Point(x, y)));
                    }
                }
            }
        }
        // TODO: Doesn't handle the case where the puzzle is empty
        let width: Width = max_x - min_x + 1;
        let height: Height = max_y - min_y + 1;
        // TODO: Doesn't handle the case where the puzzle is empty
        let shift = Point(min_x - 1, min_y - 1);

        let mut char_to_points: CharToPoints = CharToPoints::new();
        for (c, coor) in temp_points {
            let shifted_point = coor.sub(&shift); // !!!
            char_to_points.entry(c).or_default().insert(shifted_point);
        }

        (char_to_points, width, height)
    }
    let (start_chartocoors, width, height) = string_to_chartocoors(start);
    let (goal_chartocoors, goal_width, goal_height) = string_to_chartocoors(goal);

    // TODO: Handle this gracefully rather than panicking
    assert_eq!(
        start_chartocoors
            .get(&BOUNDS_CHAR)
            .unwrap_or(&Points::new()),
        goal_chartocoors.get(&BOUNDS_CHAR).unwrap_or(&Points::new()),
        "The start and goal must have the same bounds."
    );

    // TODO: Handle this gracefully rather than panicking
    assert_eq!(width, goal_width, "start_width and goal_width don't match. This should never happen, as bounds are already asserted to be the same.");
    assert_eq!(height, goal_height, "start_height and goal_height don't match. This should never happen, as bounds are already asserted to be the same.");

    fn extract_auxiliaries(
        start_chartopoints: &CharToPoints,
        goal_chartopoints: &CharToPoints,
    ) -> (
        Bounds,
        Shapekey,
        Blockstate,
        GoalShapekeyKey,
        GoalTargetOffsets,
        ReconstructionMap,
    ) {
        // TODO: Handle empty strings gracefully
        // TODO: Also maybe don't use clone, but I'm not into ownership enough to think through how to handle this
        let bounds: Shape = start_chartopoints.get(&BOUNDS_CHAR).unwrap().clone();

        let mut char_to_shape: HashMap<char, Shape> = HashMap::new();
        let mut shape_to_chars_and_offsets: HashMap<Shape, Vec<(char, Offset)>> = HashMap::new();
        let mut start_reconstruction_map: ReconstructionMap = ReconstructionMap::new();
        let mut goal_chars_startoffset_targetoffset: Vec<(char, Offset, Offset)> = Vec::new();
        for (c, start_points) in start_chartopoints.iter() {
            if c == &BOUNDS_CHAR {
                continue;
            }

            let (shape_min, _) = get_extremes(start_points);
            let shape: Shape = start_points
                .iter()
                .map(|point| point.sub(&shape_min))
                .collect();

            // This is only used to map goal-shapes to their indices in shapekey later
            char_to_shape.insert(*c, shape.clone());

            shape_to_chars_and_offsets
                .entry(shape.clone())
                .or_default()
                .push((*c, (&shape_min).into()));
            start_reconstruction_map.insert((shape, (&shape_min).into()), *c);

            if let Some(goal_points) = goal_chartopoints.get(c) {
                // TODO: Because we're currently in the start-loop, and won't separately
                // do a goal-loop, we'll silently ignore any characters in the goalstring
                // that aren't present in the startstring. If there are chars in the goalstring
                // that are not in the startstring, then that puzzle is not correctly posed,
                // which we should communicate instead.

                let (target_min, _) = get_extremes(goal_points);
                goal_chars_startoffset_targetoffset.push((
                    *c,
                    (&shape_min).into(),
                    (&target_min).into(),
                ));
            }
        }

        let mut raw_shapekey: Vec<(Shape, Vec<(char, Offset)>)> = shape_to_chars_and_offsets
            .iter()
            .map(|(shape, chars_and_offsets)| (shape.clone(), chars_and_offsets.clone()))
            .collect();

        raw_shapekey.sort_by(
            |(a_shape, a_chars_and_offsets), (b_shape, b_chars_and_offsets)| {
                let a_shape_only_for_goals = a_chars_and_offsets
                    .iter()
                    .all(|(c, _)| goal_chartopoints.get(c).is_some());
                let b_shape_only_for_goals = b_chars_and_offsets
                    .iter()
                    .all(|(c, _)| goal_chartopoints.get(c).is_some());

                // Sort shapes that are only for goals last
                if a_shape_only_for_goals && !b_shape_only_for_goals {
                    Ordering::Greater
                } else if !a_shape_only_for_goals && b_shape_only_for_goals {
                    Ordering::Less
                } else {
                    // If both or neither are only for goals, sort by size first
                    // (Idea being: If shapes are larger, we find intersections earlier)
                    // ((which probably won't matter anyway))
                    // (((but we need a total order)))
                    let a_size = a_shape.len();
                    let b_size = b_shape.len();
                    if a_size == b_size {
                        // And if sizes equal, just compare the shapes as sets
                        a_shape.cmp(b_shape)
                    } else {
                        b_size.cmp(&a_size)
                    }
                }
            },
        );

        let mut shapekey: Shapekey = raw_shapekey
            .iter()
            .map(|(shape, _)| shape.clone())
            .collect();
        shapekey.shrink_to_fit();

        // For all goal-blocks, now look up which index their shape in shapekey corresponds to
        // TODO: Should we sort this first?
        let mut goal_shapekey_key: GoalShapekeyKey = goal_chars_startoffset_targetoffset
            .iter()
            .map(|(c, _, _)| {
                raw_shapekey
                    .iter()
                    .position(|(shape, _)| shape == char_to_shape.get(c).unwrap())
                    .unwrap()
            })
            .collect();
        goal_shapekey_key.shrink_to_fit();

        let mut nongoal_offsets = raw_shapekey
            .iter()
            .map(|(_, chars_and_offsets)| {
                chars_and_offsets
                    .iter()
                    .filter(|(c, _)| goal_chartopoints.get(c).is_none())
                    .map(|(_, offset)| offset.clone())
                    .collect()
            })
            .filter(|offsets: &Offsets| !offsets.is_empty())
            .collect_vec();
        nongoal_offsets.shrink_to_fit();

        let mut goal_offsets = goal_chars_startoffset_targetoffset
            .iter()
            .map(|(_, start, _)| start.clone())
            .collect_vec();
        goal_offsets.shrink_to_fit();

        let blockstate: Blockstate = Blockstate {
            nongoal_offsets,
            goal_offsets,
        };
        let mut goal_target_offsets = goal_chars_startoffset_targetoffset
            .iter()
            .map(|(_, _, target)| target.clone())
            .collect_vec();
        goal_target_offsets.shrink_to_fit();

        (
            bounds,
            shapekey,
            blockstate,
            goal_shapekey_key,
            goal_target_offsets,
            start_reconstruction_map,
        )
    }

    let (
        bounds,
        shapekey,
        start_blockstate,
        goal_shapekey_key,
        goal_target_offsets,
        reconstruction_map,
    ) = extract_auxiliaries(&start_chartocoors, &goal_chartocoors);
    let nonintersectionkey = build_nonintersectionkey(&bounds, &shapekey, width, height);

    (
        bounds,
        shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        reconstruction_map,
        width,
        height,
    )
}

// TODO: Prove admissibility
fn misplaced_goalblocks_heuristic(
    blockstate: &Blockstate,
    goal_target_offsets: &GoalTargetOffsets,
) -> usize {
    // You might think that, rather than just counting the misplaced blocks,
    // we could also check their distance from their goal_target_offsets.
    // I implemented and benchmarked that, and it was worse than this.

    blockstate
        .goal_offsets
        .iter()
        .zip(goal_target_offsets)
        .filter(|(goal_offset, target_offset)| goal_offset != target_offset)
        .count()
}

fn solve_puzzle_path(start: &str, goal: &str) -> Vec<Blockstate> {
    let (
        _bounds,
        _shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        _start_reconstruction_map,
        _width,
        _height,
    ) = puzzle_preprocessing(start, goal);

    // If we have more than one goalblock, an astar heuristic helps speed things up.
    // Otherwise, default to usual bfs
    if goal_shapekey_key.len() > 1 {
        pathfinding::directed::astar::astar(
            &start_blockstate,
            |blockstate| {
                // TODO: More performant solution than using into_iter?
                get_neighboring_blockstates(blockstate, &nonintersectionkey, &goal_shapekey_key)
                    .into_iter()
                    .map(|blockstate| (blockstate, 1))
                    .collect_vec()
            },
            |blockstate| misplaced_goalblocks_heuristic(blockstate, &goal_target_offsets),
            |blockstate| blockstate.goal_offsets == goal_target_offsets,
        )
        .unwrap()
        .0
    } else {
        pathfinding::directed::bfs::bfs(
            &start_blockstate,
            |blockstate| {
                get_neighboring_blockstates(blockstate, &nonintersectionkey, &goal_shapekey_key)
            },
            |blockstate| blockstate.goal_offsets == goal_target_offsets,
        )
        .unwrap()
    }
}

fn solve_puzzle_minmoves(start: &str, goal: &str) -> usize {
    solve_puzzle_path(start, goal).len() - 1
}

#[wasm_bindgen]
pub fn solve_puzzle(start: &str, goal: &str) -> Vec<String> {
    let (
        bounds,
        shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        start_reconstruction_map,
        width,
        height,
    ) = puzzle_preprocessing(start, goal);

    // If we have more than one goalblock, an astar heuristic helps speed things up.
    // Otherwise, default to usual bfs
    let path = if goal_shapekey_key.len() > 1 {
        pathfinding::directed::astar::astar(
            &start_blockstate,
            |blockstate| {
                // TODO: More performant solution than using into_iter?
                get_neighboring_blockstates(blockstate, &nonintersectionkey, &goal_shapekey_key)
                    .into_iter()
                    .map(|blockstate| (blockstate, 1))
                    .collect_vec()
            },
            |blockstate| misplaced_goalblocks_heuristic(blockstate, &goal_target_offsets),
            |blockstate| blockstate.goal_offsets == goal_target_offsets,
        )
        .unwrap()
        .0
    } else {
        pathfinding::directed::bfs::bfs(
            &start_blockstate,
            |blockstate| {
                get_neighboring_blockstates(blockstate, &nonintersectionkey, &goal_shapekey_key)
            },
            |blockstate| blockstate.goal_offsets == goal_target_offsets,
        )
        .unwrap()
    };

    // Reconstruct path
    let reconstruction_map_to_string = |rm: &ReconstructionMap| -> String {
        let mut board: Vec<Vec<char>> = vec![vec![' '; width as usize]; height as usize];
        for offset in &bounds {
            board[(offset.1 - 1) as usize][(offset.0 - 1) as usize] = BOUNDS_CHAR;
        }

        for ((shape, offset), c) in rm.iter() {
            for Point(x, y) in shape.iter() {
                board[(y + offset.1 - 1) as usize][(x + offset.0 - 1) as usize] = *c;
            }
        }

        board
            .iter()
            .map(|row| row.iter().collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    };
    let mut reconstruction_map = start_reconstruction_map;
    let mut string_path = vec![reconstruction_map_to_string(&reconstruction_map)];
    let mut previous_blockstate = &path[0];
    for blockstate in path.iter().skip(1) {
        'check_single_block: {
            // nongoal-blocks
            for (shape_ix, (offsets, previous_offsets)) in blockstate
                .nongoal_offsets
                .iter()
                .zip(previous_blockstate.nongoal_offsets.iter())
                .enumerate()
            {
                if *offsets != *previous_offsets {
                    // TODO: Probably safe unwrap, but still, maybe handle this more gracefully
                    let old_offset = previous_offsets.difference(offsets).next().unwrap();
                    let new_offset = offsets.difference(previous_offsets).next().unwrap();
                    let shape = &shapekey[shape_ix];

                    // TODO: This is sooo ugly
                    let c = reconstruction_map
                        .remove(&(shape.clone(), old_offset.clone()))
                        .unwrap();

                    reconstruction_map.insert((shape.clone(), new_offset.clone()), c);

                    break 'check_single_block;
                }
            }

            // goal-blocks
            for (goalvec_ix, (offset, previous_offset)) in blockstate
                .goal_offsets
                .iter()
                .zip(previous_blockstate.goal_offsets.iter())
                .enumerate()
            {
                if *offset != *previous_offset {
                    // TODO: Probably safe unwrap, but still, maybe handle this more gracefully
                    let old_offset = previous_offset;
                    let new_offset = offset;
                    let shape = &shapekey[goal_shapekey_key[goalvec_ix]];

                    // TODO: This is sooo ugly
                    let c = reconstruction_map
                        .remove(&(shape.clone(), old_offset.clone()))
                        .unwrap();

                    reconstruction_map.insert((shape.clone(), new_offset.clone()), c);

                    break 'check_single_block;
                }
            }
        }
        // TODO: Maybe check if reconstruction_map changed at all?
        string_path.push(reconstruction_map_to_string(&reconstruction_map));

        previous_blockstate = blockstate;
    }

    string_path
}

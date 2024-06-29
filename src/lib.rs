pub mod examples;

use colored::{self, Colorize};
use std::cmp::{max, min, Ordering};
use std::collections::BTreeSet;
use std::collections::HashMap;

// TODO: Perhaps it's better to abstract most of these into structs
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
    // TODO: Should we #[inline] these?
    fn up(&self) -> Offset {
        Offset(self.0, self.1 + 1)
    }
    fn down(&self) -> Offset {
        Offset(self.0, self.1 - 1)
    }
    fn left(&self) -> Offset {
        Offset(self.0 - 1, self.1)
    }
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
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Blockstate {
    nongoal_offsets: Vec<Offsets>, // TODO: Perhaps this is better done on the stack, e.g. with https://crates.io/crates/arrayvec
    goal_offsets: Vec<Offset>,
}
type Shapekey = Vec<Shape>;
type GoalShapekeyKey = Vec<usize>; // Given an index in the Blockstate.goal_blocks vec, what is the index of its shape in the shapekey vec?
type GoalTargetOffsets = Vec<Offset>; // At what offset is a block in a goal-position?
type Offsettable<T> = Vec<Vec<T>>;

fn _intersect_coortables(a: &Offsettable<bool>, b: &Offsettable<bool>) -> Offsettable<bool> {
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
type Nonintersectionkey = Vec<Offsettable<Vec<Offsettable<bool>>>>;

const BOUNDS_CHAR: char = '.';

// TODO: This function has a bad name by now. I think it can also be moved into
// the pre-processing function?
fn extract_shapekey(
    start_chartopoints: &CharToPoints,
    goal_chartopoints: &CharToPoints,
) -> (
    Bounds,
    Shapekey,
    Blockstate,
    GoalShapekeyKey,
    GoalTargetOffsets,
) {
    // TODO: Handle empty strings gracefully
    // TODO: Also maybe don't use clone, but I'm not into ownership enough to think through how to handle this
    let bounds: Shape = start_chartopoints.get(&BOUNDS_CHAR).unwrap().clone();

    let mut char_to_shape: HashMap<char, Shape> = HashMap::new();
    let mut shape_to_chars_and_offsets: HashMap<Shape, Vec<(char, Offset)>> = HashMap::new();
    let mut goal_chars_startoffset_targetoffset: Vec<(char, Offset, Offset)> = Vec::new();
    for (c, start_points) in start_chartopoints.iter() {
        if c == &BOUNDS_CHAR {
            continue;
        }

        // TODO: Should get_extremes return two Points: (Point(min), Point(max)) instead?
        let (shape_min, _) = get_extremes(start_points);
        let shape: Shape = start_points
            .iter()
            .map(|point| point.sub(&shape_min))
            .collect();

        // This is only used to map goal-shapes to their indices in shapekey later
        char_to_shape.insert(*c, shape.clone());

        shape_to_chars_and_offsets
            .entry(shape)
            .or_default()
            .push((*c, (&shape_min).into()));
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

    let shapekey: Shapekey = raw_shapekey
        .iter()
        .map(|(shape, _)| shape.clone())
        .collect();
    // For all goal-blocks, now look up which index their shape in shapekey corresponds to
    // TODO: Should we sort this first?
    let goal_shapekey_key: GoalShapekeyKey = goal_chars_startoffset_targetoffset
        .iter()
        .map(|(c, _, _)| {
            raw_shapekey
                .iter()
                .position(|(shape, _)| shape == char_to_shape.get(c).unwrap())
                .unwrap()
        })
        .collect();

    let blockstate: Blockstate = Blockstate {
        nongoal_offsets: raw_shapekey
            .iter()
            .map(|(_, chars_and_offsets)| {
                chars_and_offsets
                    .iter()
                    .filter(|(c, _)| goal_chartopoints.get(c).is_none())
                    .map(|(_, offset)| offset.clone())
                    .collect()
            })
            .filter(|offsets: &Offsets| !offsets.is_empty())
            .collect(),
        goal_offsets: goal_chars_startoffset_targetoffset
            .iter()
            .map(|(_, start, _)| start.clone())
            .collect(),
    };
    let goal_target_offsets = goal_chars_startoffset_targetoffset
        .iter()
        .map(|(_, _, target)| target.clone())
        .collect();

    (
        bounds,
        shapekey,
        blockstate,
        goal_shapekey_key,
        goal_target_offsets,
    )
}

fn build_nonintersectionkey(
    bounds: &Bounds,
    shapekey: &Shapekey,
    width: Width,
    height: Height,
) -> Nonintersectionkey {
    // Brace yourselves

    let mut nik = Nonintersectionkey::new();
    for shape_a in shapekey {
        let mut nik_a: Offsettable<Vec<Offsettable<bool>>> = Offsettable::new();
        // TODO: Do we really have to include 0 and width+1?
        for xa in 0..=(width + 1) {
            let mut nik_ax: Vec<Vec<Offsettable<bool>>> = Vec::new();
            for ya in 0..=(height + 1) {
                let mut nik_axy: Vec<Offsettable<bool>> = Vec::new();
                let shift_a = Point(xa, ya);

                // TODO: Extract into shift-function
                let shifted_a: Points = shape_a.iter().map(|p| p.add(&shift_a)).collect();
                if shifted_a.is_subset(bounds) {
                    // Let the fun begin
                    for shape_b in shapekey {
                        let mut nik_axy_b: Offsettable<bool> = Offsettable::new();
                        for xb in 0..=(width + 1) {
                            let mut nik_axy_bx: Vec<bool> = Vec::new();
                            for yb in 0..=(height + 1) {
                                let shift_b = Point(xb, yb);
                                // TODO: Extract into shift-function
                                let shifted_b: Points =
                                    shape_b.iter().map(|p| p.add(&shift_b)).collect();

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
    goal_shapekey_key: &GoalShapekeyKey,
) -> Vec<Blockstate> {
    // TODO: Better-yet than the next todo, could we pass around the
    // nik-intersections that are described in the next todo?
    // Expressed in the context of Point-Sets rather than niks (and
    // a concept still worth exploring if the idea doesn't work for niks,
    // as it definitely *does* work for Point-Sets):
    // - First, take the union U of all block-coordinates
    // - To dfs a block, remove it from U, and test its legality by checking
    //   if it shifted doesn't intersect U
    //   (Should bounds be baked into U, or be kept in a separate union, or
    //    should we maybe do nonintersectionkeys but just for bounds?)
    // - And all this can be sped up in the future, because we can store
    //   the union U in the blockstate, so that we don't have to calculate it
    //   from scratch every time we call get_neighboring_blockstates!!

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
    */ // for now, just this bad implementation of dfs:
    // TODO: I have no idea if `&dyn Fn(Offset) -> bool` is the right signature as I didn't learn about `&dyn` yet
    let dfs_general = |initial_offset: Offset, is_legal: &dyn Fn(&Offset) -> bool| {
        let mut legal_offsets: Vec<Offset> = Vec::new();
        let mut seen_offsets: BTreeSet<Offset> = BTreeSet::new(); // TODO: Different data structures?
        seen_offsets.insert(initial_offset.clone());
        let mut stack: Vec<Offset> = vec![initial_offset];
        while let Some(offset) = stack.pop() {
            // TODO: Maybe this can be made faster by not going back in
            // the direction we just came from
            for new_offset in [offset.up(), offset.down(), offset.left(), offset.right()] {
                // TODO: Can we check seen_offsets membership without cloning?
                // TODO: Which of these checks is faster? And shouldn't we be using pathfinding::directed::dfs::dfs instead?
                if is_legal(&new_offset) && seen_offsets.insert(new_offset.clone()) {
                    legal_offsets.push(new_offset.clone());
                    stack.push(new_offset);
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
                    if !nonintersectionkey[shape_ix][offset.0 as usize][offset.1 as usize]
                        [movingshape_ix][offsety.0 as usize][offsety.1 as usize]
                    {
                        return false;
                    }
                }
            }
            for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
                let shapekey_ix = goal_shapekey_key[goalvec_ix];
                if !nonintersectionkey[shapekey_ix][offset.0 as usize][offset.1 as usize]
                    [movingshape_ix][offsety.0 as usize][offsety.1 as usize]
                {
                    return false;
                }
            }
            for offset in trimmed_movingshape_offsets {
                if !nonintersectionkey[movingshape_ix][offset.0 as usize][offset.1 as usize]
                    [movingshape_ix][offsety.0 as usize][offsety.1 as usize]
                {
                    return false;
                }
            }
            true
        };
        dfs_general(moving_offset, &is_legal)
    };

    let dfs_goal = |moving_goalvec_ix: usize, moving_offset: Offset| -> Vec<Offset> {
        let moving_shapekey_ix = goal_shapekey_key[moving_goalvec_ix];
        let is_legal = |offsetty: &Offset| -> bool {
            // TODO: This function assumes there are other blocks on the field, because
            // we currently use their nonintersectionkeys to additionally verify that
            // offsety is in-bounds
            for (shape_ix, shape_offsets) in blockstate.nongoal_offsets.iter().enumerate() {
                for offset in shape_offsets {
                    // TODO: How bad are these "as usize" conversions?
                    // If they're really bad, I might just end up using usize as the type for Coor
                    if !nonintersectionkey[shape_ix][offset.0 as usize][offset.1 as usize]
                        [moving_shapekey_ix][offsetty.0 as usize][offsetty.1 as usize]
                    {
                        return false;
                    }
                }
            }
            for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
                if goalvec_ix == moving_goalvec_ix {
                    continue;
                }
                let shapekey_ix = goal_shapekey_key[goalvec_ix];
                if !nonintersectionkey[shapekey_ix][offset.0 as usize][offset.1 as usize]
                    [moving_shapekey_ix][offsetty.0 as usize][offsetty.1 as usize]
                {
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
    for (shapekey_ix, offsets) in blockstate.nongoal_offsets.iter().enumerate() {
        for offset in offsets.iter() {
            let mut trimmed_shape_offsets = offsets.clone();
            trimmed_shape_offsets.remove(offset);
            // TODO: avoid .clone() here?
            for mutated_offset in dfs_nongoal(shapekey_ix, offset.clone(), &trimmed_shape_offsets) {
                let mut mutated_shape_offsets = trimmed_shape_offsets.clone();
                mutated_shape_offsets.insert(mutated_offset);
                let mut new_blockstate = blockstate.clone();
                new_blockstate.nongoal_offsets[shapekey_ix] = mutated_shape_offsets;
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
    for (shapekey_ix, offset) in goal_shapekey_key.iter().zip(blockstate.goal_offsets.iter()) {
        let block: Points = shapekey[*shapekey_ix]
            .iter()
            .map(|p| p.add(&offset.into()))
            .collect();
        blocks.push(block);
    }
    blocks.sort(); // ensures consistent indices
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
    Coor,
    Coor,
) {
    fn string_to_chartocoors(s: &str) -> (CharToPoints, Width, Height) {
        // TODO: Rename the function, and variables that later call this function
        // (the current name is horrible)

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

    let (bounds, shapekey, start_blockstate, goal_shapekey_key, goal_target_offsets) =
        extract_shapekey(&start_chartocoors, &goal_chartocoors);
    let nonintersectionkey = build_nonintersectionkey(&bounds, &shapekey, width, height);

    (
        bounds,
        shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        width,
        height,
    )
}

// TODO: Implement A*
// TODO: Add tests

pub fn solve_puzzle(start: &str, goal: &str) {
    let (
        bounds,
        shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        width,
        height,
    ) = puzzle_preprocessing(start, goal);

    if false {
        print_puzzle(
            &bounds,
            &shapekey,
            &start_blockstate,
            &goal_shapekey_key,
            width,
            height,
        );
    }

    let path = pathfinding::directed::bfs::bfs(
        &start_blockstate,
        |blockstate| {
            get_neighboring_blockstates(blockstate, &nonintersectionkey, &goal_shapekey_key)
        },
        |blockstate| blockstate.goal_offsets == goal_target_offsets,
    )
    .unwrap();

    println!("{}", path.len() - 1);

    assert!(path.len() < 1000); // TODO: Remove this
}

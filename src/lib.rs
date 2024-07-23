pub mod examples;

use bitvec::prelude::*;
use itertools::Itertools;
use std::cmp::{max, min, Ordering};
use std::collections::{BTreeMap, BTreeSet, HashMap};
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
type CharToPoints = BTreeMap<char, Points>;
type ReconstructionMap = HashMap<(Shape, Offset), char>;

type Shape = BTreeSet<Point>; // nonempty. min-x == 0, min-y == 0

// TODO: Make this a struct method?
fn get_extremes(coordinates_set: &Points) -> (Point, Point) {
    // Extract min-x and min-y.
    // Assumes that coordinatesSet is nonempty.
    let mut min_x = Coor::MAX;
    let mut min_y = Coor::MAX;
    let mut max_x = Coor::MIN;
    let mut max_y = Coor::MIN;
    for point in coordinates_set {
        max_x = max(max_x, point.0);
        max_y = max(max_y, point.1);
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
    fn up(&self) -> Self {
        Self(self.0, self.1 + 1)
    }
    #[inline]
    fn down(&self) -> Self {
        Self(self.0, self.1 - 1)
    }
    #[inline]
    fn left(&self) -> Self {
        Self(self.0 - 1, self.1)
    }
    #[inline]
    fn right(&self) -> Self {
        Self(self.0 + 1, self.1)
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
#[cfg_attr(test, derive(PartialEq, Debug))]
struct Nonintersectionkey {
    width: usize,
    nik: ShapevecForNik<Vec<ShapevecForNik<BitVec>>>,
}
impl std::ops::Index<(usize, &Offset, usize, &Offset)> for Nonintersectionkey {
    type Output = bool;

    #[inline]
    fn index(
        &self,
        (shape_ix_a, Offset(xa, ya), shape_ix_b, Offset(xb, yb)): (usize, &Offset, usize, &Offset),
    ) -> &Self::Output {
        &self.nik[shape_ix_a][*xa as usize + (*ya as usize) * self.width][shape_ix_b]
            [*xb as usize + (*yb as usize) * self.width]
    }
}
impl Nonintersectionkey {
    // TODO: Can (should) you rewrite the bounds-check into something else?
    fn abuse_this_datastructure_for_in_bounds_check(
        &self,
        shape_ix: usize,
        Offset(xa, ya): &Offset,
    ) -> bool {
        !self.nik[shape_ix][*xa as usize + (*ya as usize) * self.width].is_empty()
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum SolvePuzzleError {
    MismatchedBounds,
    MismatchedGoalShapes(char),
    GoalblockWithoutStartingblock(char),
    WidthTooLarge,
    HeightTooLarge,
    EmptyPuzzle,
}
impl std::fmt::Display for SolvePuzzleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // TODO: Better error messages.
            SolvePuzzleError::MismatchedBounds => write!(f, "Start-bounds don't match goal-bounds."),
            SolvePuzzleError::MismatchedGoalShapes(c) => write!(f, "The shape of the goal-block '{c}' in the start-configuration doesn't match its shape in the goal-configuration."),
            SolvePuzzleError::GoalblockWithoutStartingblock(c) => write!(f, "The block '{c}' is in the goal-configuration, but not in the start-configuration."),
            SolvePuzzleError::WidthTooLarge => write!(f, "The width of the puzzle is too large to fit into the data-type the solver uses. If you encounter this while trying to solve an actual puzzle, please file an issue."),
            SolvePuzzleError::HeightTooLarge => write!(f, "The height of the puzzle is too large to fit into the data-type the solver uses. If you encounter this while trying to solve an actual puzzle, please file an issue."),
            SolvePuzzleError::EmptyPuzzle => write!(f, "Please add some blocks."),
        }
    }
}
impl std::error::Error for SolvePuzzleError {}
impl From<SolvePuzzleError> for JsValue {
    fn from(val: SolvePuzzleError) -> Self {
        JsValue::from(val.to_string())
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
    Nonintersectionkey {
        width: width as usize + 2,
        nik: nikvec,
    }
}

// TODO: I have no idea if `&dyn Fn(Offset) -> bool` is the right signature as I didn't learn about `&dyn` yet
fn dfs_general(
    initial_offset: &Offset,
    is_legal: &dyn Fn(&Offset) -> bool,
) -> impl Iterator<Item = Offset> {
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

    // TODO: Refactor to not use a vector at all
    legal_offsets.into_iter()
}

// To keep track of which block we just moved,
// so that we don't move it again on the next turn.
#[derive(Clone)]
enum Justmoved {
    Nongoal(usize, Offset), // Shape index, offset
    Goal(usize),            // goalvec-ix
    Nothing,
}
#[derive(Clone)]
struct BlockstateJustmoved {
    blockstate: Blockstate,
    justmoved: Justmoved,
}
impl PartialEq for BlockstateJustmoved {
    fn eq(&self, other: &Self) -> bool {
        self.blockstate == other.blockstate
    }
}
impl Eq for BlockstateJustmoved {}
impl std::hash::Hash for BlockstateJustmoved {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.blockstate.hash(state);
    }
}

fn get_neighboring_blockstates(
    BlockstateJustmoved {
        blockstate,
        justmoved,
    }: &BlockstateJustmoved,
    nonintersectionkey: &Nonintersectionkey,
    goal_shapekey_key: &GoalShapekeyKey,
) -> Vec<BlockstateJustmoved> {
    fn dfs_nongoal<'a>(
        moving: impl Iterator<Item = (usize, &'a Offset, &'a Offsets)> + 'a,
        blockstate: &'a Blockstate,
        nonintersectionkey: &'a Nonintersectionkey,
        goal_shapekey_key: &'a GoalShapekeyKey,
    ) -> impl Iterator<Item = BlockstateJustmoved> + 'a {
        moving.flat_map(
            move |(moving_shape_ix, moving_offset, offsets_with_same_shape_ix)| {
                // TODO: Workhorse outside of closure?
                let mut trimmed_movingshape_offsets = offsets_with_same_shape_ix.clone();
                trimmed_movingshape_offsets.remove(moving_offset);
                let is_legal = |offsety: &Offset| -> bool {
                    for (shape_ix, shape_offsets) in blockstate.nongoal_offsets.iter().enumerate() {
                        if shape_ix == moving_shape_ix {
                            continue;
                        }
                        for offset in shape_offsets {
                            // TODO: How bad are these "as usize" conversions?
                            // If they're really bad, I might just end up using usize as the type for Coor
                            if !nonintersectionkey[(shape_ix, offset, moving_shape_ix, offsety)] {
                                return false;
                            }
                        }
                    }
                    for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
                        let shape_ix = goal_shapekey_key[goalvec_ix];
                        if !nonintersectionkey[(shape_ix, offset, moving_shape_ix, offsety)] {
                            return false;
                        }
                    }
                    for offset in &trimmed_movingshape_offsets {
                        if !nonintersectionkey[(moving_shape_ix, offset, moving_shape_ix, offsety)]
                        {
                            return false;
                        }
                    }
                    true
                };
                dfs_general(moving_offset, &is_legal).map(move |mutated_offset| {
                    // TODO: Workhorse outside of closure?
                    let mut mutated_trimmed_offsets = trimmed_movingshape_offsets.clone();
                    mutated_trimmed_offsets.insert(mutated_offset.clone());
                    let mut new_blockstate = blockstate.clone();
                    new_blockstate.nongoal_offsets[moving_shape_ix] = mutated_trimmed_offsets;
                    let justmoved = Justmoved::Nongoal(moving_shape_ix, mutated_offset);
                    BlockstateJustmoved {
                        blockstate: new_blockstate,
                        justmoved,
                    }
                })
            },
        )
    }

    fn dfs_goal<'a>(
        moving: impl Iterator<Item = (usize, &'a Offset)> + 'a,
        blockstate: &'a Blockstate,
        nonintersectionkey: &'a Nonintersectionkey,
        goal_shapekey_key: &'a GoalShapekeyKey,
    ) -> impl Iterator<Item = BlockstateJustmoved> + 'a {
        moving.flat_map(move |(moving_goalvec_ix, moving_offset)| {
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
            let justmoved = Justmoved::Goal(moving_goalvec_ix);
            dfs_general(moving_offset, &is_legal).map(move |mutated_offset| {
                let mut new_blockstate = blockstate.clone();
                new_blockstate.goal_offsets[moving_goalvec_ix] = mutated_offset;
                BlockstateJustmoved {
                    blockstate: new_blockstate,
                    justmoved: justmoved.clone(),
                }
            })
        })
    }

    // It's okay to gather all these into a vector rather than a set,
    // because all neighbors WILL be unique.
    let mut neighboring_blockstates: Vec<BlockstateJustmoved> = Vec::new();
    // Start with blockstate.goal_offsets first, to hopefully find the goal a little sooner

    // TODO: Lots of code-duplication over these branches. Should we export it into
    // two closures, one for adding the nongoal offsets, and one for adding the
    // the goal-offsets?
    match justmoved {
        Justmoved::Nongoal(justmoved_shape_ix, justmoved_offset) => {
            neighboring_blockstates.extend(dfs_goal(
                blockstate.goal_offsets.iter().enumerate(),
                blockstate,
                nonintersectionkey,
                goal_shapekey_key,
            ));
            neighboring_blockstates.extend(dfs_nongoal(
                blockstate
                    .nongoal_offsets
                    .iter()
                    .enumerate()
                    .flat_map(|(shape_ix, offsets)| {
                        offsets
                            .into_iter()
                            .filter(move |offset| {
                                // TODO: The .filter checks shape_ix != *justmoved_shape_ix for every single offset.
                                // Can we do better?
                                // And do we even need to? Branch-predictor should do some solid work here,
                                // and these checks really are cheap.
                                shape_ix != *justmoved_shape_ix || *offset != justmoved_offset
                            })
                            .map(move |offset| (shape_ix, offset, offsets))
                    }),
                &blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
            ));
        }
        Justmoved::Goal(moved_goalvec_ix) => {
            neighboring_blockstates.extend(dfs_goal(
                blockstate
                    .goal_offsets
                    .iter()
                    .enumerate()
                    .filter(|(goalvec_ix, _)| *goalvec_ix != *moved_goalvec_ix),
                blockstate,
                nonintersectionkey,
                goal_shapekey_key,
            ));
            neighboring_blockstates.extend(dfs_nongoal(
                blockstate
                    .nongoal_offsets
                    .iter()
                    .enumerate()
                    .flat_map(|(shape_ix, offsets)| {
                        offsets
                            .into_iter()
                            .map(move |offset| (shape_ix, offset, offsets))
                    }),
                &blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
            ));
        }
        Justmoved::Nothing => {
            neighboring_blockstates.extend(dfs_goal(
                blockstate.goal_offsets.iter().enumerate(),
                blockstate,
                nonintersectionkey,
                goal_shapekey_key,
            ));
            neighboring_blockstates.extend(dfs_nongoal(
                blockstate
                    .nongoal_offsets
                    .iter()
                    .enumerate()
                    .flat_map(|(shape_ix, offsets)| {
                        offsets
                            .into_iter()
                            .map(move |offset| (shape_ix, offset, offsets))
                    }),
                &blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
            ));
        }
    }
    neighboring_blockstates
}

fn _print_puzzle(
    bounds: &Bounds,
    shapekey: &Shapekey,
    blockstate: &Blockstate,
    goal_shapekey_key: &GoalShapekeyKey,
    width: Width,
    height: Height,
) {
    // Create vec of blocks:
    let mut blocks: Vec<Points> = Vec::new();
    for (shape, offsets) in shapekey.iter().zip(blockstate.nongoal_offsets.iter()) {
        for offset in offsets {
            // TODO: Just create shift-shape method already
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

    for y in 0..=(height + 1) {
        let mut line0: Vec<&str> = Vec::new();
        let mut line1: Vec<&str> = Vec::new();
        for x in 0..=(width + 1) {
            // Find block that contains (x, y)
            let option_block_ix: Option<usize> =
                blocks.iter().position(|block| block.contains(&Point(x, y)));
            match option_block_ix {
                Some(block_ix) => {
                    let t = y > 0 && blocks[block_ix].contains(&Point(x, y - 1));
                    let b = blocks[block_ix].contains(&Point(x, y + 1));
                    let l = x > 0 && blocks[block_ix].contains(&Point(x - 1, y));
                    let r = blocks[block_ix].contains(&Point(x + 1, y));
                    let tl = y > 0 && x > 0 && blocks[block_ix].contains(&Point(x - 1, y - 1));
                    let tr = y > 0 && blocks[block_ix].contains(&Point(x + 1, y - 1));
                    let bl = x > 0 && blocks[block_ix].contains(&Point(x - 1, y + 1));
                    let br = blocks[block_ix].contains(&Point(x + 1, y + 1));
                    line0.push(match (t, l, tl) {
                        (false, false, _) => "╭─",
                        (true, false, _) => "│ ",
                        (false, true, _) => "──",
                        (true, true, false) => "╯ ",
                        (true, true, true) => "  ",
                    });
                    line0.push(match (t, r, tr) {
                        (false, false, _) => "─╮",
                        (true, false, _) => " │",
                        (false, true, _) => "──",
                        (true, true, false) => " ╰",
                        (true, true, true) => "  ",
                    });
                    line1.push(match (b, l, bl) {
                        (false, false, _) => "╰─",
                        (true, false, _) => "│ ",
                        (false, true, _) => "──",
                        (true, true, false) => "╮ ",
                        (true, true, true) => "  ",
                    });
                    line1.push(match (b, r, br) {
                        (false, false, _) => "─╯",
                        (true, false, _) => " │",
                        (false, true, _) => "──",
                        (true, true, false) => " ╭",
                        (true, true, true) => "  ",
                    });
                }
                None => {
                    if bounds.contains(&Point(x, y)) {
                        line0.push("    ");
                        line1.push("    ");
                    } else {
                        line0.push("████");
                        line1.push("████");
                    }
                }
            }
        }
        println!("{}", line0.concat());
        println!("{}", line1.concat());
    }
}

#[cfg_attr(test, derive(PartialEq, Debug))]
struct Auxiliaries {
    bounds: Bounds,
    shapekey: Shapekey,
    start_blockstate: Blockstate,
    nonintersectionkey: Nonintersectionkey,
    goal_shapekey_key: GoalShapekeyKey,
    goal_target_offsets: GoalTargetOffsets,
    reconstruction_map: ReconstructionMap,
    width: Width,
    height: Height,
}

#[cfg_attr(test, derive(PartialEq, Debug))]
enum PreprocessingOutput {
    EmptyPuzzle,
    ProperPuzzle(Auxiliaries),
}

fn preprocess_proper_puzzle(
    start_chartopoints: &CharToPoints,
    goal_chartopoints: &CharToPoints,
    width: Width,
    height: Height,
) -> Result<PreprocessingOutput, SolvePuzzleError> {
    // Assumes that bounds exist in CharToPoints.

    // Check that bounds match.
    let start_bounds = &start_chartopoints[&BOUNDS_CHAR];
    let goal_bounds = &goal_chartopoints[&BOUNDS_CHAR];
    if start_bounds != goal_bounds {
        return Err(SolvePuzzleError::MismatchedBounds);
    }

    let bounds: Shape = start_chartopoints[&BOUNDS_CHAR].clone();

    let mut char_to_shape: HashMap<char, Shape> = HashMap::new();
    let mut shape_to_chars_and_offsets: HashMap<Shape, Vec<(char, Offset)>> = HashMap::new();
    let mut reconstruction_map: ReconstructionMap = ReconstructionMap::new();
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
        reconstruction_map.insert((shape, (&shape_min).into()), *c);

        if let Some(goal_points) = goal_chartopoints.get(c) {
            let (target_min, _) = get_extremes(goal_points);
            goal_chars_startoffset_targetoffset.push((
                *c,
                (&shape_min).into(),
                (&target_min).into(),
            ));
        }
    }
    let mut char_to_goalshape: HashMap<char, Shape> = HashMap::new();
    for (c, goal_points) in goal_chartopoints.iter() {
        // We already know the bounds match, so we don't need to care about those
        if c == &BOUNDS_CHAR {
            continue;
        }
        let (shape_min, _) = get_extremes(goal_points);
        let shape: Shape = goal_points
            .iter()
            .map(|point| point.sub(&shape_min))
            .collect();
        char_to_goalshape.insert(*c, shape.clone());
    }
    // Check that the start and goal shapes are the same
    for (c, goalshape) in char_to_goalshape.iter() {
        let startshape = char_to_shape
            .get(c)
            .ok_or(SolvePuzzleError::GoalblockWithoutStartingblock(*c))?;
        if startshape != goalshape {
            return Err(SolvePuzzleError::MismatchedGoalShapes(*c));
        }
    }

    let mut raw_shapekey: Vec<(Shape, Vec<(char, Offset)>)> = shape_to_chars_and_offsets
        .iter()
        .map(|(shape, chars_and_offsets)| (shape.clone(), chars_and_offsets.clone()))
        .collect();

    raw_shapekey.sort_unstable_by(
        |(a_shape, a_chars_and_offsets), (b_shape, b_chars_and_offsets)| {
            let a_shape_only_for_goals = a_chars_and_offsets
                .iter()
                .all(|(c, _)| goal_chartopoints.contains_key(c));
            let b_shape_only_for_goals = b_chars_and_offsets
                .iter()
                .all(|(c, _)| goal_chartopoints.contains_key(c));

            // Sort shapes that are only for goals last
            if a_shape_only_for_goals && !b_shape_only_for_goals {
                Ordering::Greater
            } else if !a_shape_only_for_goals && b_shape_only_for_goals {
                Ordering::Less
            } else {
                // If both or neither are only for goals, sort by size firste idea
                // being: If shapes are larger, we find intersections earlier
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

    // Sort goal-blocks by size, too.
    // This almost never matters, because most puzzles have at most one
    // goal-block anyway, and those that don't often have goal-blocks of
    // the same size.
    goal_chars_startoffset_targetoffset.sort_unstable_by(|(ca, _, _), (cb, _, _)| {
        let a_size = char_to_shape[ca].len();
        let b_size = char_to_shape[cb].len();
        if a_size == b_size {
            ca.cmp(cb)
        } else {
            b_size.cmp(&a_size)
        }
    });

    // For all goal-blocks, now look up which index their shape in shapekey corresponds to
    let mut goal_shapekey_key: GoalShapekeyKey = goal_chars_startoffset_targetoffset
        .iter()
        .map(|(c, _, _)| {
            raw_shapekey
                .iter()
                .position(|(shape, _)| *shape == char_to_shape[c])
                .unwrap()
        })
        .collect();
    goal_shapekey_key.shrink_to_fit();

    let mut nongoal_offsets = raw_shapekey
        .iter()
        .map(|(_, chars_and_offsets)| {
            chars_and_offsets
                .iter()
                .filter(|(c, _)| !goal_chartopoints.contains_key(c))
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

    let nonintersectionkey = build_nonintersectionkey(&bounds, &shapekey, width, height);

    Ok(PreprocessingOutput::ProperPuzzle(Auxiliaries {
        bounds,
        shapekey,
        start_blockstate: blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        reconstruction_map,
        width,
        height,
    }))
}

// TODO: This code is awful, and awfully named.
enum StringToCharToPointsResult {
    EmptyPuzzle,
    ProperPuzzle(CharToPoints, Width, Height),
}
fn string_to_chartopoints(s: &str) -> Result<StringToCharToPointsResult, SolvePuzzleError> {
    let mut min_x = usize::MAX;
    let mut min_y = usize::MAX;
    let mut max_x = usize::MIN;
    let mut max_y = usize::MIN;

    let mut char_annotated_coordinates: Vec<(char, (usize, usize))> = Vec::new();

    for (y, l) in s.lines().enumerate() {
        for (x, c) in l.chars().enumerate() {
            if !c.is_whitespace() {
                min_x = min(min_x, x);
                min_y = min(min_y, y);
                max_x = max(max_x, x);
                max_y = max(max_y, y);
                char_annotated_coordinates.push((c, (x, y)));
                if c != BOUNDS_CHAR {
                    char_annotated_coordinates.push((BOUNDS_CHAR, (x, y)));
                }
            }
        }
    }
    if char_annotated_coordinates.is_empty() {
        return Ok(StringToCharToPointsResult::EmptyPuzzle);
    }
    let width_usize: usize = max_x + 1 - min_x;
    let height_usize: usize = max_y + 1 - min_y;

    // Try converting to Coor, and otherwise throw SolvePuzzleError::WidthTooLarge
    let width: Width = width_usize
        .try_into()
        .map_err(|_| SolvePuzzleError::WidthTooLarge)?;
    let height: Height = height_usize
        .try_into()
        .map_err(|_| SolvePuzzleError::HeightTooLarge)?;

    let mut char_to_points: CharToPoints = CharToPoints::new();
    for (c, pair) in char_annotated_coordinates {
        let (x_usize, y_usize) = (pair.0 + 1 - min_x, pair.1 + 1 - min_y);
        let x: Coor = x_usize
            .try_into()
            .map_err(|_| SolvePuzzleError::WidthTooLarge)?;
        let y: Coor = y_usize
            .try_into()
            .map_err(|_| SolvePuzzleError::HeightTooLarge)?;
        let point = Point(x, y);
        char_to_points.entry(c).or_default().insert(point);
    }

    Ok(StringToCharToPointsResult::ProperPuzzle(
        char_to_points,
        width,
        height,
    ))
}

fn preprocessing(start: &str, goal: &str) -> Result<PreprocessingOutput, SolvePuzzleError> {
    let start_stc_result = string_to_chartopoints(start)?;
    let goal_stc_result = string_to_chartopoints(goal)?;

    match (start_stc_result, goal_stc_result) {
        (StringToCharToPointsResult::EmptyPuzzle, StringToCharToPointsResult::EmptyPuzzle) => {
            Ok(PreprocessingOutput::EmptyPuzzle)
        }
        (
            StringToCharToPointsResult::ProperPuzzle(start_chartopoints, width, height),
            StringToCharToPointsResult::ProperPuzzle(goal_chartopoints, _goal_width, _goal_height),
        ) => preprocess_proper_puzzle(&start_chartopoints, &goal_chartopoints, width, height),
        _ => Err(SolvePuzzleError::MismatchedBounds),
    }
}

fn misplaced_goalblocks_heuristic(
    blockstate: &Blockstate,
    goal_target_offsets: &GoalTargetOffsets,
) -> usize {
    // You might think that, rather than just counting the misplaced blocks,
    // we could also check their distance from their goal_target_offsets.
    // I implemented and benchmarked that, and it was worse than this.

    // TODO: Should we carry heuristic-stats throughout the A* nodes and
    // only update the counter if a goal-block enters or leaves its target
    // position?

    blockstate
        .goal_offsets
        .iter()
        .zip(goal_target_offsets)
        .filter(|(goal_offset, target_offset)| goal_offset != target_offset)
        .count()
}

fn solution_from_auxiliaries(
    Auxiliaries {
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        ..
    }: &Auxiliaries,
) -> Option<Vec<Blockstate>> {
    // If we have more than one goalblock, an astar heuristic helps speed things up.
    // Otherwise, default to usual bfs
    match goal_shapekey_key.len() {
        0 => Some(vec![start_blockstate.clone()]),
        1 => {
            // If this goalblock is also the only block in the puzzle, we can't just call
            // the usual BFS function, as that assumes we have at least two distinct
            // blocks in the puzzle to do its in-bounds-check.
            if start_blockstate.nongoal_offsets.is_empty() {
                if start_blockstate.goal_offsets == *goal_target_offsets {
                    return Some(vec![start_blockstate.clone()]);
                }
                let beginning_offset = start_blockstate.goal_offsets[0].clone();
                let is_legal = |offset: &Offset| -> bool {
                    nonintersectionkey.abuse_this_datastructure_for_in_bounds_check(0, offset)
                };
                let neighbors = dfs_general(&beginning_offset, &is_legal).collect_vec();
                if neighbors.contains(&goal_target_offsets[0]) {
                    let mut goal_blockstate = start_blockstate.clone();
                    goal_blockstate.goal_offsets.clone_from(goal_target_offsets);
                    Some(vec![start_blockstate.clone(), goal_blockstate])
                } else {
                    None
                }
            } else {
                pathfinding::directed::bfs::bfs(
                    &BlockstateJustmoved {
                        blockstate: start_blockstate.clone(),
                        justmoved: Justmoved::Nothing,
                    },
                    |blockstate| {
                        get_neighboring_blockstates(
                            blockstate,
                            nonintersectionkey,
                            goal_shapekey_key,
                        )
                    },
                    |BlockstateJustmoved { blockstate, .. }| {
                        blockstate.goal_offsets == *goal_target_offsets
                    },
                )
                .map(|path| {
                    path.into_iter()
                        .map(|BlockstateJustmoved { blockstate, .. }| blockstate)
                        .collect_vec()
                })
            }
        }
        _ => pathfinding::directed::astar::astar(
            &BlockstateJustmoved {
                blockstate: start_blockstate.clone(),
                justmoved: Justmoved::Nothing,
            },
            |blockstate| {
                // TODO: More performant solution than using into_iter?
                get_neighboring_blockstates(blockstate, nonintersectionkey, goal_shapekey_key)
                    .into_iter()
                    .map(|blockstate| (blockstate, 1))
                    .collect_vec()
            },
            |BlockstateJustmoved { blockstate, .. }| {
                misplaced_goalblocks_heuristic(blockstate, goal_target_offsets)
            },
            |BlockstateJustmoved { blockstate, .. }| {
                blockstate.goal_offsets == *goal_target_offsets
            },
        )
        .map(|path| {
            path.0
                .into_iter()
                .map(|BlockstateJustmoved { blockstate, .. }| blockstate)
                .collect_vec()
        }),
    }
}

fn solve_puzzle_path(start: &str, goal: &str) -> Result<Option<Vec<Blockstate>>, SolvePuzzleError> {
    match preprocessing(start, goal)? {
        PreprocessingOutput::ProperPuzzle(auxiliaries) => {
            Ok(solution_from_auxiliaries(&auxiliaries))
        }
        PreprocessingOutput::EmptyPuzzle => Err(SolvePuzzleError::EmptyPuzzle),
    }
}

pub fn solve_puzzle_minmoves(start: &str, goal: &str) -> Result<Option<usize>, SolvePuzzleError> {
    let maybe_path = solve_puzzle_path(start, goal)?;
    Ok(maybe_path.map(|path| path.len() - 1))
}

#[wasm_bindgen]
pub fn solve_puzzle(start: &str, goal: &str) -> Result<Option<Vec<String>>, SolvePuzzleError> {
    match preprocessing(start, goal)? {
        PreprocessingOutput::EmptyPuzzle => Err(SolvePuzzleError::EmptyPuzzle),
        PreprocessingOutput::ProperPuzzle(auxiliaries) => {
            Ok(solution_from_auxiliaries(&auxiliaries).map(|path| {
                let width = auxiliaries.width;
                let height = auxiliaries.height;
                let bounds = auxiliaries.bounds;
                let shapekey = auxiliaries.shapekey;
                let goal_shapekey_key = auxiliaries.goal_shapekey_key;

                // Reconstruct path
                let reconstruction_map_to_string = |rm: &ReconstructionMap| -> String {
                    let mut board: Vec<Vec<char>> =
                        vec![vec![' '; width as usize]; height as usize];
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
                let mut reconstruction_map = auxiliaries.reconstruction_map;
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
                                let old_offset =
                                    previous_offsets.difference(offsets).next().unwrap();
                                let new_offset =
                                    offsets.difference(previous_offsets).next().unwrap();
                                let shape = &shapekey[shape_ix];

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
                                let old_offset = previous_offset;
                                let new_offset = offset;
                                let shape = &shapekey[goal_shapekey_key[goalvec_ix]];

                                let c = reconstruction_map
                                    .remove(&(shape.clone(), old_offset.clone()))
                                    .unwrap();

                                reconstruction_map.insert((shape.clone(), new_offset.clone()), c);

                                break 'check_single_block;
                            }
                        }
                    }
                    string_path.push(reconstruction_map_to_string(&reconstruction_map));

                    previous_blockstate = blockstate;
                }

                string_path
            }))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_preprocessing() {
        assert_eq!(solve_puzzle("", ""), Err(SolvePuzzleError::EmptyPuzzle));
        assert_eq!(preprocessing("", ""), Ok(PreprocessingOutput::EmptyPuzzle));
        assert_eq!(
            preprocessing("          ", "    "),
            Ok(PreprocessingOutput::EmptyPuzzle)
        );
        assert_eq!(
            preprocessing(
                "
             	     
        ",
                "
             	     

             	     
            "
            ),
            Ok(PreprocessingOutput::EmptyPuzzle)
        );
        assert_eq!(
            preprocessing("a", "b"),
            Err(SolvePuzzleError::GoalblockWithoutStartingblock('b'))
        );
        let wide_str: &str = &"a".repeat(u8::MAX as usize + 1);
        let tall_str: &str = &"a
        "
        .repeat(u8::MAX as usize + 1);
        assert_eq!(
            preprocessing(wide_str, wide_str),
            Err(SolvePuzzleError::WidthTooLarge)
        );
        assert_eq!(
            preprocessing(tall_str, tall_str),
            Err(SolvePuzzleError::HeightTooLarge)
        );
        assert_eq!(
            preprocessing("aa", "a."),
            Err(SolvePuzzleError::MismatchedGoalShapes('a'))
        );
        assert_eq!(
            preprocessing("a.", "a"),
            Err(SolvePuzzleError::MismatchedBounds)
        );
    }
}

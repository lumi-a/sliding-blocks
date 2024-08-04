//! Solving sliding-block puzzles like the [`ones from the Professor Layton games`]
//!
//! [`ones from the Professor Layton games`](https://layton.fandom.com/wiki/Category:Sliding)

#![warn(
    clippy::all,
    clippy::pedantic,
    clippy::nursery,
    clippy::cargo,
    clippy::decimal_literal_representation,
    clippy::empty_enum_variants_with_brackets,
    clippy::empty_structs_with_brackets,
    clippy::exhaustive_enums,
    clippy::exhaustive_structs,
    clippy::get_unwrap,
    clippy::if_then_some_else_none,
    clippy::indexing_slicing,
    clippy::infinite_loop,
    clippy::integer_division,
    clippy::missing_assert_message,
    clippy::missing_asserts_for_indexing,
    clippy::missing_inline_in_public_items,
    clippy::multiple_inherent_impl,
    clippy::print_stdout,
    clippy::redundant_type_annotations,
    clippy::shadow_unrelated,
    clippy::shadow_same,
    clippy::shadow_reuse,
    clippy::single_call_fn,
    clippy::string_lit_chars_any,
    clippy::try_err,
    clippy::undocumented_unsafe_blocks,
    clippy::unneeded_field_pattern,
    clippy::unwrap_in_result,
    clippy::unwrap_used,
    clippy::use_debug,
    clippy::wildcard_enum_match_arm,
    clippy::default_numeric_fallback
)]

pub mod examples;

use bitvec::prelude::*;
use core::cmp::{max, min, Ordering};
use core::ops::{Add, Sub};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use vec_collections::{AbstractVecSet, VecSet};
use wasm_bindgen::prelude::*;

/// Type of a single coordinate. If changing this type, also change the type of `Offset`.
type Coor = u8;

/// Alias for `Coor`, for readability.
type Width = Coor;

/// Alias for `Coor`, for readability.
type Height = Coor;

/// A "global" coordinate, in contrast to the `Offset` of a `Shape`.
/// Over time, these two types diverged a lot, for what I hope (but
/// don't believe) to be good reasons. I might implement `Point`
/// as a wrapper around `Offset` at some point. (TODO)
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Copy)]
struct Point(Coor, Coor);
impl Add for Point {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0, self.1 + other.1)
    }
}
impl Sub for Point {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0, self.1 - other.1)
    }
}
impl From<Offset> for Point {
    fn from(offset: Offset) -> Self {
        Self(offset.x(), offset.y())
    }
}

/// A collection of Points. It's a `BTreeSet` rather than any
/// other set, to implement `Ord`.
type Points = BTreeSet<Point>;

/// Shift a collection of points by an `Offset`.
fn shift_points(points: &Points, offset: Offset) -> Points {
    points.iter().map(|point| *point + offset.into()).collect()
}

/// A map from a `char` to the cells it occupies. Used in pre-processing.
type CharToPoints = BTreeMap<char, Points>;

/// A map from a given `Shape` and its `Offset` to its `char`. Used in mutable form
/// for outputting the path the solution takes with the input-`char`s.
type ReconstructionMap = HashMap<(Shape, Offset), char>;

/// A shape. This must not be empty.
type Shape = BTreeSet<Point>;

/// Type-alias for `Shape` to describe bounds.
/// Its minimum x-value and minimum y-value will always be 1.
type Bounds = Shape;

/// Extract the min/max x/y values from `Points`.
fn get_points_dimensions(coordinates_set: &Points) -> (Point, Point) {
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

/// An offset (of a shape). We'll later assume that each shape has,
/// as long as its in-bounds, a min-x-offset of at least 1, and
/// min-y-offset of at least 1.
///
/// Implemented as a `u16` to support fast hashing, addition, and
/// other int-niceties.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Copy)]
struct Offset(u16);
impl Offset {
    /// Create an `Offset` from the given coordinates.
    #[inline]
    fn new(x: Coor, y: Coor) -> Self {
        Self(((u16::from(x)) << 8) | u16::from(y))
    }
    /// Extract x-coordinate from self.
    #[inline]
    const fn x(self) -> Coor {
        (self.0 >> 8) as Coor
    }
    /// Extract y-coordinate from self.
    #[inline]
    const fn y(self) -> Coor {
        (self.0 & 0xff) as Coor
    }
    /// Extract x-coordinate from self as usize.
    #[inline]
    const fn x_usize(self) -> usize {
        (self.0 >> 8) as usize
    }
    /// Extract y-coordinate from self as usize.
    #[inline]
    const fn y_usize(self) -> usize {
        (self.0 & 0xff) as usize
    }
    /// Return new offset shifted up by one cell.
    #[inline]
    const fn up(self) -> Self {
        Self(self.0 + 1)
    }
    /// Return new offset shifted down by one cell.
    #[inline]
    const fn down(self) -> Self {
        Self(self.0 - 1)
    }
    /// Return new offset shifted left by one cell.
    #[inline]
    const fn left(self) -> Self {
        Self(self.0 - 0x100)
    }
    /// Return new offset shifted right by one cell.
    #[inline]
    const fn right(self) -> Self {
        Self(self.0 + 0x100)
    }
}
impl From<&Point> for Offset {
    fn from(point: &Point) -> Self {
        Self::new(point.0, point.1)
    }
}
impl tinyset::Fits64 for Offset {
    unsafe fn from_u64(x: u64) -> Self {
        Self(x as u16)
    }
    fn to_u64(self) -> u64 {
        self.0 as u64
    }
}

/// A collection of `Offset`s, used for storing
/// offsets of the same shape in one blockstate.
/// TODO: Experiment with the initial capacity.
type Offsets = VecSet<[Offset; 16]>;

/// The blockstates, stored as a tuple:
/// - A `Vec`, each entry being a collection of non-goal-`Offsets`.
/// - A `Vec`, each entry being a goal-`Offset`.
/// By "non-goal-offset" we mean that the offset refers
/// to a block that does not have a specified target-position.
/// "goal-offsets" correspond to blocks with specified target-positions.
///
/// For an offset-collection in the `nongoal_offsets` (i.e. an entry in the
/// `nongoal_offsets` vec), all blocks in that offset-collection
/// have the same shape. This is the sole reason we store them
/// as a *set*: If two non-goal-blocks are in interchangable
/// positions, that defines equivalent blockstates.
/// Depending on the puzzle, this can lead to considerable speedups.
/// This cannot be applied to `goal_offsets`, as this won't lead
/// to equivalent blockstates (in particular, two goal-blocks of
/// the same shape might have different goal-positions, so a solved-blockstate
/// would then be equivalent to a not-yet-solved-blockstate, i.e. we'd
/// get messed up.)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Blockstate {
    nongoal_offsets: Vec<Offsets>,
    goal_offsets: Vec<Offset>,
}

/// Keeps track of shapes, for all `Blockstate`s. Always has
/// length at least `nongoal_offsets.len()`. All blocks corresponding
/// to entry `i` in `nongoal_offsets`have shape `shapekey[i]`.
type Shapekey = Vec<Shape>;

/// To also keep track of goal-block-shapes, we store what
/// indices of `Shapekey` we should look them up in.
/// If the shape of a goal-block is not shared with any
/// non-goal-block, that shape is appended to the end of
/// the shapekey.
/// For a goal-block at entry `i` in `goal_offsets`, its shape
/// is `shapekey[goal_shapekeykey[i]]`.
type GoalShapekeyKey = Vec<usize>;

/// Stores the target-offsets each goal-block needs to reach.
/// For a goal-block at entry `i` in `goal_offsets`, its target-offset
/// is `goal_target_offsets[i]`.
type GoalTargetOffsets = Vec<Offset>;

/// Type-alias to make the definition of `Nonintersectionkey`
/// more readable.
type ShapevecForNik<T> = Vec<T>;

/// So that we not always have to check whether two blocks of
/// two given shape and two given offsets intersect each other,
/// we calculate all that in pre-processing and store it in
/// what I called a Nonintersectionkey, a horrible name (TODO)
/// Its `width` keeps track of the width of the board, so that we
/// don't have to index over x and y separately.
/// Its `nik` attribute stores the values. For a block A of shape-index
/// `shape_a` at offset `(xa, ya)` and a block B of shape-index `shape_b` at
/// offset `(xb, yb)`:
/// - `nik[shape_a][xa+ya*width]` is a nonempty vec iff A is in-bounds,
///   and if so:
/// - `nik[shape_a][xa+ya*width][shape_b][xb+yb*width]` is true iff B
///   is in-bounds and does not intersect A.
///
/// The first assumption saves on memory and will always be satisfied
/// later on.
/// Storing whether B is in-bounds here speeds up checks later and means
/// we don't have to keep track of bounds that much, but it also means
/// that we can't rely on Nonintersectionkey for in-bounds-checking if
/// there is only a single block in the puzzle. But, of course, if there's
/// only a single block, the puzzle is rather easy to solve.
///
/// We only allow indices `(x,y)` between `(0,0)` and `(width+1, height+1)`.
/// This is done to avoid checking for edge-cases where a block is right at the
/// edge of the bounds, so that our addition doesn't overflow. It also
/// means we have to let in-bounds-blocks have a min-offset of `(1,1)`,
/// which is only a mild inconvenience.
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
        (shape_ix_a, a, shape_ix_b, b): (usize, &Offset, usize, &Offset),
    ) -> &Self::Output {
        // TODO: Store Offset in a different form so that we don't need to convert them this way?
        // Perhaps we could also store ShapevecForNik as a HashMap?
        &self.nik[shape_ix_a][a.x_usize() + a.y_usize() * self.width][shape_ix_b]
            [b.x_usize() + b.y_usize() * self.width]
    }
}
impl Nonintersectionkey {
    fn abuse_this_datastructure_for_in_bounds_check(
        &self,
        shape_ix: usize,
        offset: Offset,
    ) -> bool {
        !self.nik[shape_ix][offset.x_usize() + offset.y_usize() * self.width].is_empty()
    }
}

/// An error that can be returned by the solver.
#[derive(Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum SolvePuzzleError {
    MismatchedBounds,
    MismatchedGoalShapes(char),
    GoalblockWithoutStartingblock(char),
    WidthTooLarge,
    HeightTooLarge,
    EmptyPuzzle,
}
impl std::fmt::Display for SolvePuzzleError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MismatchedBounds => write!(f, "Start-bounds don't match goal-bounds."),
            Self::MismatchedGoalShapes(c) => write!(f, "The shape of the goal-block '{c}' in the start-configuration doesn't match its shape in the goal-configuration."),
            Self::GoalblockWithoutStartingblock(c) => write!(f, "The block '{c}' is in the goal-configuration, but not in the start-configuration."),
            Self::WidthTooLarge => write!(f, "The width of the puzzle is too large to fit into the data-type the solver uses. If you encounter this while trying to solve an actual puzzle, please file an issue."),
            Self::HeightTooLarge => write!(f, "The height of the puzzle is too large to fit into the data-type the solver uses. If you encounter this while trying to solve an actual puzzle, please file an issue."),
            Self::EmptyPuzzle => write!(f, "Please add some blocks."),
        }
    }
}
impl std::error::Error for SolvePuzzleError {}
impl From<SolvePuzzleError> for JsValue {
    #[inline]
    fn from(val: SolvePuzzleError) -> Self {
        Self::from(val.to_string())
    }
}

const BOUNDS_CHAR: char = '.';

fn build_nonintersectionkey(
    bounds: &Bounds,
    shapekey: &Shapekey,
    width: Width,
    height: Height,
) -> Nonintersectionkey {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    struct RelativeOffset(isize, isize); // (isize, isize) because Coor is too small and I want negatives
    let mut relative_nik: HashMap<(usize, usize, RelativeOffset), bool> = HashMap::new();
    let mut inbounds: HashMap<(usize, Offset), bool> = HashMap::new();

    // Brace yourselves
    let mut nikvec: ShapevecForNik<Vec<ShapevecForNik<BitVec>>> = Vec::new();
    for (shape_ix_a, shape_a) in shapekey.iter().enumerate() {
        let mut nik_a: Vec<ShapevecForNik<BitVec>> = Vec::new();
        for ya in 0..=(height + 1) {
            for xa in 0..=(width + 1) {
                let mut nik_a_xy: ShapevecForNik<BitVec> = ShapevecForNik::new();

                let shift_a = Offset::new(xa, ya);
                let shifted_a: Points = shift_points(shape_a, shift_a);
                if shifted_a.is_subset(bounds) {
                    // Let the fun begin
                    for (shape_ix_b, shape_b) in shapekey.iter().enumerate() {
                        let mut nik_a_xy_b: BitVec = BitVec::new();
                        for yb in 0..=(height + 1) {
                            for xb in 0..=(width + 1) {
                                let relative_offset = RelativeOffset(
                                    xb as isize + width as isize - xa as isize,
                                    yb as isize + height as isize - ya as isize,
                                );
                                let shift_b = Offset::new(xb, yb);
                                nik_a_xy_b.push(
                                    *relative_nik
                                        .entry((shape_ix_a, shape_ix_b, relative_offset))
                                        .or_insert({
                                            let shifted_b: Points = shift_points(shape_b, shift_b);
                                            shifted_b.is_disjoint(&shifted_a)
                                        })
                                        && *inbounds.entry((shape_ix_b, shift_b)).or_insert({
                                            let shifted_b: Points = shift_points(shape_b, shift_b);
                                            shifted_b.is_subset(bounds)
                                        }),
                                );
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

fn dfs_general(
    initial_offset: &Offset,
    is_legal: &dyn Fn(&Offset) -> bool,
) -> impl Iterator<Item = Offset> {
    let mut legal_offsets: Vec<Offset> = Vec::new();
    let mut seen_offsets: VecSet<[Offset; 8]> = VecSet::empty();
    seen_offsets.insert(*initial_offset);

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
            seen_offsets.insert(new_offset);
            legal_offsets.push(new_offset);
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
            if is_legal(&new_offset) && seen_offsets.insert(new_offset) {
                legal_offsets.push(new_offset);
                stack.push((new_offset, new_dir));
            }
        }
    }

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
    // TODO: This is quite ugly, taking in a reference to the other offsets
    // for _every_ offset. The only reason we do this is to filter out
    // in case of BlockstateJustmoved::nongoal.
    fn dfs_nongoal<'a>(
        moving: impl Iterator<Item = (usize, Offset, &'a Offsets)> + 'a,
        blockstate: &'a Blockstate,
        nonintersectionkey: &'a Nonintersectionkey,
        goal_shapekey_key: &'a GoalShapekeyKey,
    ) -> impl Iterator<Item = BlockstateJustmoved> + 'a {
        moving.flat_map(
            move |(moving_shape_ix, moving_offset, offsets_with_same_shape_ix)| {
                let mut trimmed_movingshape_offsets = offsets_with_same_shape_ix.clone();
                trimmed_movingshape_offsets.remove(&moving_offset);
                let is_legal = |offsety: &Offset| -> bool {
                    for (shape_ix, shape_offsets) in blockstate.nongoal_offsets.iter().enumerate() {
                        if shape_ix == moving_shape_ix {
                            continue;
                        }
                        for offset in shape_offsets.iter() {
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
                    for offset in trimmed_movingshape_offsets.iter() {
                        if !nonintersectionkey[(moving_shape_ix, offset, moving_shape_ix, offsety)]
                        {
                            return false;
                        }
                    }
                    true
                };
                dfs_general(&moving_offset, &is_legal).map(move |mutated_offset| {
                    let mut mutated_trimmed_offsets = trimmed_movingshape_offsets.clone();
                    mutated_trimmed_offsets.insert(mutated_offset);
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
                for (shape_ix, shape_offsets) in blockstate.nongoal_offsets.iter().enumerate() {
                    for offset in shape_offsets.iter() {
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
                            .iter()
                            .filter(move |offset| {
                                // TODO: The .filter checks shape_ix != *justmoved_shape_ix for every single offset.
                                // Can we do better?
                                // And do we even need to? Branch-predictor should do some solid work here,
                                // and these checks really are cheap.
                                shape_ix != *justmoved_shape_ix || **offset != *justmoved_offset
                            })
                            .map(move |offset| (shape_ix, *offset, offsets))
                    }),
                blockstate,
                nonintersectionkey,
                goal_shapekey_key,
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
                            .iter()
                            .map(move |offset| (shape_ix, *offset, offsets))
                    }),
                blockstate,
                nonintersectionkey,
                goal_shapekey_key,
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
                            .iter()
                            .map(move |offset| (shape_ix, *offset, offsets))
                    }),
                blockstate,
                nonintersectionkey,
                goal_shapekey_key,
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
        for offset in offsets.iter() {
            let block: Points = shift_points(shape, *offset);
            blocks.push(block);
        }
    }
    for (shape_ix, offset) in goal_shapekey_key.iter().zip(blockstate.goal_offsets.iter()) {
        let block: Points = shapekey[*shape_ix]
            .iter()
            .map(|p| *p + (*offset).into())
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

        let (shape_min, _) = get_points_dimensions(start_points);
        let shape: Shape = start_points
            .iter()
            .map(|point| *point - shape_min)
            .collect();

        // This is only used to map goal-shapes to their indices in shapekey later
        char_to_shape.insert(*c, shape.clone());

        shape_to_chars_and_offsets
            .entry(shape.clone())
            .or_default()
            .push((*c, (&shape_min).into()));
        reconstruction_map.insert((shape, (&shape_min).into()), *c);

        if let Some(goal_points) = goal_chartopoints.get(c) {
            let (target_min, _) = get_points_dimensions(goal_points);
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
        let (shape_min, _) = get_points_dimensions(goal_points);
        let shape: Shape = goal_points.iter().map(|point| *point - shape_min).collect();
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
                .map(|(_, offset)| *offset)
                .collect()
        })
        .filter(|offsets: &Offsets| !offsets.is_empty())
        .collect_vec();
    nongoal_offsets.shrink_to_fit();

    let mut goal_offsets = goal_chars_startoffset_targetoffset
        .iter()
        .map(|(_, start, _)| *start)
        .collect_vec();
    goal_offsets.shrink_to_fit();

    let blockstate: Blockstate = Blockstate {
        nongoal_offsets,
        goal_offsets,
    };
    let mut goal_target_offsets = goal_chars_startoffset_targetoffset
        .iter()
        .map(|(_, _, target)| *target)
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
                let beginning_offset = start_blockstate.goal_offsets[0];
                let is_legal = |offset: &Offset| -> bool {
                    nonintersectionkey.abuse_this_datastructure_for_in_bounds_check(0, *offset)
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
                get_neighboring_blockstates(blockstate, nonintersectionkey, goal_shapekey_key)
                    .into_iter()
                    .map(|blockstate| (blockstate, 1))
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
                    for point in &bounds {
                        board[(point.1 - 1) as usize][(point.0 - 1) as usize] = BOUNDS_CHAR;
                    }

                    for ((shape, offset), c) in rm.iter() {
                        for Point(x, y) in shape.iter() {
                            board[(y + offset.y() - 1) as usize][(x + offset.x() - 1) as usize] =
                                *c;
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
                                let old_offset = previous_offsets
                                    .iter()
                                    .find(|o| !offsets.contains(o))
                                    .unwrap();
                                let new_offset = offsets
                                    .iter()
                                    .find(|o| !previous_offsets.contains(o))
                                    .unwrap();
                                let shape = &shapekey[shape_ix];

                                let c = reconstruction_map
                                    .remove(&(shape.clone(), *old_offset))
                                    .unwrap();

                                reconstruction_map.insert((shape.clone(), *new_offset), c);

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
                                    .remove(&(shape.clone(), *old_offset))
                                    .unwrap();

                                reconstruction_map.insert((shape.clone(), *new_offset), c);

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

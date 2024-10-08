//! Solving sliding-block puzzles like the [ones from the](https://layton.fandom.com/wiki/Category:Slide) [Professor Layton games](https://layton.fandom.com/wiki/Category:Sliding).
//! To solve puzzles, use [`solve_puzzle`] or [`solve_puzzle_minmoves`].
//! For example-puzzles, see the [`examples`] crate.

#![warn(
    clippy::all,
    clippy::nursery,
    clippy::cargo,
    clippy::exhaustive_enums,
    clippy::exhaustive_structs,
    clippy::redundant_type_annotations
)]

pub mod examples;

use bitvec::prelude::*;
use core::cmp::{max, min, Ordering};
use core::ops::{Add, Sub};
use itertools::Itertools;
use rustc_hash::FxHashMap as HashMap;
use smallvec::SmallVec;
use std::collections::{BTreeMap, BTreeSet};
use vec_collections::{AbstractVecSet, VecSet};
use wasm_bindgen::prelude::*;

/// Type of a single coordinate. If changing this type, also change the type of `Offset`.
type Coor = u8;

/// Alias for `Coor`, for readability.
type Width = Coor;

/// Alias for `Coor`, for readability.
type Height = Coor;

/// A "global" coordinate, in contrast to the `Offset` of a `Shape`.
type Point = Offset;

/// A collection of Points. It's a `BTreeSet` rather than any
/// other set, to implement `Ord`.
type Points = BTreeSet<Point>;

/// Shift a collection of points by an `Offset`.
fn shift_points(points: &Points, offset: Offset) -> Points {
    points.iter().map(|point| *point + offset).collect()
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
        max_x = max(max_x, point.x());
        max_y = max(max_y, point.y());
        min_x = min(min_x, point.x());
        min_y = min(min_y, point.y());
    }
    (Offset::new(min_x, min_y), Offset::new(max_x, max_y))
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
        Self(self.0.wrapping_add(1))
    }
    /// Return new offset shifted down by one cell.
    #[inline]
    const fn down(self) -> Self {
        Self(self.0.wrapping_sub(1))
    }
    /// Return new offset shifted left by one cell.
    #[inline]
    const fn left(self) -> Self {
        Self(self.0.wrapping_sub(0x100))
    }
    /// Return new offset shifted right by one cell.
    #[inline]
    const fn right(self) -> Self {
        Self(self.0.wrapping_add(0x100))
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
impl Add for Offset {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0.wrapping_add(rhs.0))
    }
}
impl Sub for Offset {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0.wrapping_sub(rhs.0))
    }
}

/// A collection of `Offset`s, used for storing
/// offsets of the same shape in one blockstate.
type Offsets = VecSet<[Offset; 4]>;

/// The blockstates, represented as a tuple of:
/// - A `Vec`, each entry being a set of non-goal-`Offsets`.
/// - A `Vec`, each entry being a goal-`Offset`.
///
/// By "non-goal-`Offset`" we mean that the `Offset` refers
/// to a block that does not have a specified target-position.
/// A "goal-`Offset`" corresponds to a block with a specified
/// target-position.
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
    /// A `Vec`, each entry being a set of non-goal-`Offsets`.
    nongoal_offsets: NongoalOffsets,
    /// A `Vec`, each entry being a goal-`Offset`.
    goal_offsets: GoalOffsets,
}
/// A `Vec`, each entry being a set of non-goal-`Offsets`.
type NongoalOffsets = SmallVec<[Offsets; 8]>;
/// A `Vec`, each entry being a goal-`Offset`.
type GoalOffsets = SmallVec<[Offset; 8]>;

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
type GoalTargetOffsets = GoalOffsets;

/// Type-alias to make the definition of `Nonintersectionkey`
/// more readable.
type ShapevecForNik<T> = Vec<T>;

/// So that we not always have to check whether two blocks of
/// two given shape and two given offsets intersect each other,
/// we calculate all that in pre-processing and store it in
/// what I called a Nonintersectionkey, a horrible name.
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
    /// The width of the board
    width: usize,
    /// The lookup-table
    nik: ShapevecForNik<Vec<ShapevecForNik<BitVec>>>,
}
impl core::ops::Index<(usize, Offset, usize, Offset)> for Nonintersectionkey {
    type Output = bool;

    #[inline]
    fn index(
        &self,
        (shape_ix_a, offset_a, shape_ix_b, offset_b): (usize, Offset, usize, Offset),
    ) -> &Self::Output {
        &self.nik[shape_ix_a][offset_a.x_usize() + offset_a.y_usize() * self.width][shape_ix_b]
            [offset_b.x_usize() + offset_b.y_usize() * self.width]
    }
}
impl Nonintersectionkey {
    /// For a given shape and offset, `nik[shape][offset]` is nonempty
    /// iff the block is in-bounds, which we exploit in this function
    /// to check whether a block is in-bounds.
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
impl core::fmt::Display for SolvePuzzleError {
    #[inline]
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
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

/// The character delimiting unoccupied bounds in the puzzle
const BOUNDS_CHAR: char = '.';

/// Run a DFS to determine all legal movements of a block.
/// `initial_offset` is the offset of the block in the start-configuration.
/// `is_legal` determines whether an offset is legal.
fn dfs_general<FN>(initial_offset: Offset, is_legal: FN) -> impl Iterator<Item = Offset>
where
    FN: Fn(Offset) -> bool,
{
    /// To eliminate some backtracking, store
    /// which direction the block just moved
    enum CameFrom {
        Up,
        Down,
        Left,
        Right,
    }

    // We'll return legal_offsets in the end
    let mut legal_offsets: SmallVec<[Offset; 8]> = Default::default();
    let mut seen_offsets: VecSet<[Offset; 16]> = VecSet::empty();
    seen_offsets.insert(initial_offset);

    let mut stack: SmallVec<[(Offset, CameFrom); 8]> = Default::default();

    // Initial setup:
    for (new_offset, new_dir) in [
        (initial_offset.up(), CameFrom::Up),
        (initial_offset.down(), CameFrom::Down),
        (initial_offset.left(), CameFrom::Left),
        (initial_offset.right(), CameFrom::Right),
    ] {
        // No need to check if seen_offsets contains them, as we know it won't
        if is_legal(new_offset) {
            seen_offsets.insert(new_offset);
            legal_offsets.push(new_offset);
            stack.push((new_offset, new_dir));
        }
    }

    // Run DFS
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
            if is_legal(new_offset) && seen_offsets.insert(new_offset) {
                legal_offsets.push(new_offset);
                stack.push((new_offset, new_dir));
            }
        }
    }

    legal_offsets.into_iter()
}

/// Keep track of which block we just moved,
/// so that we don't try to move it again on the next turn.
#[derive(Clone)]
enum Justmoved {
    /// A non-goal-block was just moved, identified by
    /// its `(shape-ix, offset)`.
    Nongoal(usize, Offset),
    /// A goal-block was just moved, identified by
    /// its index in the goalvec
    Goal(usize),
    /// No block was just moved (starting-position)
    Nothing,
}

/// A `Blockstate` together with a `Justmoved`.
/// This is the node-type we'll traverse in the end.
#[derive(Clone)]
struct BlockstateJustmoved {
    /// The node's `Blockstate`
    blockstate: Blockstate,
    /// The node's `Justmoved`
    justmoved: Justmoved,
}
impl PartialEq for BlockstateJustmoved {
    fn eq(&self, other: &Self) -> bool {
        self.blockstate == other.blockstate
    }
}
impl Eq for BlockstateJustmoved {}
impl core::hash::Hash for BlockstateJustmoved {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.blockstate.hash(state);
    }
}

/// Get the neighbors of a `BlockstateJustmoved`, by trying
/// to move each block, which uses [`dfs_general`].
/// TODO: If you're willing to have a lot of code-duplication, it might
/// be worth the effort to write a separate function for puzzles that
/// have exactly one goal-block.
fn get_neighboring_blockstates(
    BlockstateJustmoved {
        blockstate,
        justmoved,
    }: &BlockstateJustmoved,
    nonintersectionkey: &Nonintersectionkey,
    goal_shapekey_key: &GoalShapekeyKey,
) -> Vec<BlockstateJustmoved> {
    // It might be possible to speed up queries using dynamic-programming
    // on the nonintersectionkey: In the end,
    //  `left_nonintersection[shape][i] ∩ right_nonintersection[shape][?-i]
    //   ∩ [(nik[shape][x1][y1][shape]∩…∩(nik[shape][x(i-1)][y(i-1)][shape]))
    //     ∩ (nik[shape][x(i+1)][y(i+1)][shape]∩…∩(nik[shape][x?][y?][shape]))
    //   ]`
    // will describe exactly the positions that block `i` of shape `shape` is allowed to move to.
    // This requires linear (in the number of total blocks) amount of computations for
    // left_nonintersection and right_nonintersection.
    // I did try to implement this though, and the benchmarks showed it was slower.
    // Perhaps that's because the current implementation (which loops against each block until
    // it finds a violation) terminates so early, or it's because my code sucks.

    /// Run DFS by moving several nongoal-shapes.
    /// Here, `moving` is an `Iterator` over `(shape-ix, offset, offsets-with-same-shape-ix)`.
    fn dfs_nongoal<'a>(
        moving: impl Iterator<Item = (usize, Offset, &'a Offsets)> + 'a,
        blockstate: &'a Blockstate,
        nonintersectionkey: &'a Nonintersectionkey,
        goal_shapekey_key: &'a GoalShapekeyKey,
    ) -> impl Iterator<Item = BlockstateJustmoved> + 'a {
        // TODO: This is quite ugly, taking in a reference to the other offsets
        // for _every_ offset. The only reason we do this is to filter out
        // in case of BlockstateJustmoved::nongoal.

        // Try moving each block in `moving`
        moving.flat_map(
            move |(moving_shape_ix, moving_offset, offsets_with_same_shape_ix)| {
                // Remove `moving_offset` from `offsets_with_same_shape_ix`,
                // we'll reinsert its new positions later
                let mut trimmed_movingshape_offsets = offsets_with_same_shape_ix.clone();
                trimmed_movingshape_offsets.remove(&moving_offset);

                let tmo = trimmed_movingshape_offsets.clone();
                let is_legal = move |offsety: Offset| -> bool {
                    blockstate.nongoal_offsets.iter().enumerate().all(
                        // Do we intersect any nongoal-blocks of different shape?
                        |(shape_ix, shape_offsets)| {
                            shape_ix == moving_shape_ix
                                || shape_offsets.iter().all(|offset| {
                                    nonintersectionkey
                                        [(shape_ix, *offset, moving_shape_ix, offsety)]
                                })
                        },
                    ) && blockstate.goal_offsets.iter().enumerate().all(
                        // Do we intersect any goal-blocks?
                        |(goalvec_ix, offset)| {
                            nonintersectionkey[(
                                goal_shapekey_key[goalvec_ix],
                                *offset,
                                moving_shape_ix,
                                offsety,
                            )]
                        },
                    ) && tmo.iter().all(|offset| {
                        // Do we intersect any nongoal-blocks of the same shape?
                        nonintersectionkey[(moving_shape_ix, *offset, moving_shape_ix, offsety)]
                    })
                };
                // Run dfs_general using this `is_legal` function
                dfs_general(moving_offset, is_legal).map(move |mutated_offset| {
                    // Reinsert the new offset
                    let mut mutated_trimmed_offsets = trimmed_movingshape_offsets.clone();
                    mutated_trimmed_offsets.insert(mutated_offset);

                    // Put the new set of offsets into a new blockstate
                    let mut new_blockstate = blockstate.clone();
                    new_blockstate.nongoal_offsets[moving_shape_ix] = mutated_trimmed_offsets;

                    // Keep track of which block we just moved
                    let justmoved = Justmoved::Nongoal(moving_shape_ix, mutated_offset);
                    BlockstateJustmoved {
                        blockstate: new_blockstate,
                        justmoved,
                    }
                })
            },
        )
    }

    /// Run DFS by moving several goal-shapes.
    /// Here, `moving` is an `Iterator` over `(goalvec_ix, offset)`.
    fn dfs_goal<'a>(
        moving: impl Iterator<Item = (usize, &'a Offset)> + 'a,
        blockstate: &'a Blockstate,
        nonintersectionkey: &'a Nonintersectionkey,
        goal_shapekey_key: &'a GoalShapekeyKey,
    ) -> impl Iterator<Item = BlockstateJustmoved> + 'a {
        // Try moving each block in `moving`
        moving.flat_map(move |(moving_goalvec_ix, moving_offset)| {
            let moving_shape_ix = goal_shapekey_key[moving_goalvec_ix];
            let is_legal = move |offsety: Offset| -> bool {
                blockstate.nongoal_offsets.iter().enumerate().all(
                    // Check whether we intersect any nongoal-blocks
                    |(shape_ix, shape_offsets)| {
                        shape_offsets.iter().all(|offset| {
                            nonintersectionkey[(shape_ix, *offset, moving_shape_ix, offsety)]
                        })
                    },
                ) && blockstate
                    .goal_offsets
                    .iter()
                    .enumerate()
                    .all(|(goalvec_ix, offset)| {
                        // Check whether we intersect any other goal-blocks
                        goalvec_ix == moving_goalvec_ix
                            || nonintersectionkey[(
                                goal_shapekey_key[goalvec_ix],
                                *offset,
                                moving_shape_ix,
                                offsety,
                            )]
                    })
            };
            // Record which block we just moved
            let justmoved = Justmoved::Goal(moving_goalvec_ix);

            // Run DFS with that `is_legal`-function
            dfs_general(*moving_offset, is_legal).map(move |mutated_offset| {
                // Build a new `Blockstate` using the new goalblock-offset
                let mut new_blockstate = blockstate.clone();
                new_blockstate.goal_offsets[moving_goalvec_ix] = mutated_offset;

                BlockstateJustmoved {
                    blockstate: new_blockstate,
                    justmoved: justmoved.clone(),
                }
            })
        })
    }

    // We'll return `neighboring_blockstates` in the end.
    // It's okay to gather all these into a vector rather than a set,
    // because all neighbors will be unique.
    let mut neighboring_blockstates: Vec<BlockstateJustmoved> = Vec::new();

    // This is quite ugly code-duplication, but each branch of the following match first
    // runs `dfs_goal` and then `dfs_nongoal`. We run `dfs_nongoal` first, to hopefully
    // find the target a little sooner.
    // Depending on the branch, we throw out the block we just moved.
    match *justmoved {
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
                                shape_ix != justmoved_shape_ix || **offset != justmoved_offset
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
                    .filter(|(goalvec_ix, _)| *goalvec_ix != moved_goalvec_ix),
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

    // Return all the neighbors
    neighboring_blockstates
}

/// Pretty-print a puzzle represented by auxiliary-data.
fn _print_puzzle(
    bounds: &Bounds,
    shapekey: &Shapekey,
    blockstate: &Blockstate,
    goal_shapekey_key: &GoalShapekeyKey,
    width: Width,
    height: Height,
) {
    // Create vec of actual blocks.
    // We do not care about keeping shapes together,
    // nor whether a block is non-goal or goal.
    let blocks: Vec<Points> = {
        let mut blocks: Vec<Points> = Vec::new();

        // Add non-goal-blocks
        for (shape, offsets) in shapekey.iter().zip(blockstate.nongoal_offsets.iter()) {
            for offset in offsets.iter() {
                let block: Points = shift_points(shape, *offset);
                blocks.push(block);
            }
        }

        // Add goal-blocks
        for (shape_ix, offset) in goal_shapekey_key.iter().zip(blockstate.goal_offsets.iter()) {
            let block: Points = shapekey[*shape_ix].iter().map(|p| *p + (*offset)).collect();
            blocks.push(block);
        }
        blocks
    };

    // Iterate over all cells of the puzzle
    for y in 0..=(height + 1) {
        // We'll be printing two lines per cell-row
        let mut line0: Vec<&str> = Vec::new();
        let mut line1: Vec<&str> = Vec::new();
        for x in 0..=(width + 1) {
            // Find block that contains (x, y)
            let c = Offset::new(x, y);
            let option_block_ix: Option<usize> = blocks.iter().position(|block| block.contains(&c));

            match option_block_ix {
                // The coordinate is in a block
                Some(block_ix) => {
                    // Depending on whether the block at `block_ix` also occupies
                    // surrounding cells, print different box-drawing chars.
                    let t = y > 0 && blocks[block_ix].contains(&c.up());
                    let b = blocks[block_ix].contains(&c.down());
                    let l = x > 0 && blocks[block_ix].contains(&c.left());
                    let r = blocks[block_ix].contains(&c.right());
                    let tl = y > 0 && x > 0 && blocks[block_ix].contains(&c.left().down());
                    let tr = y > 0 && blocks[block_ix].contains(&c.right().up());
                    let bl = x > 0 && blocks[block_ix].contains(&c.down().left());
                    let br = blocks[block_ix].contains(&c.right().down());
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
                // The coordinate is not in a block
                None => {
                    // Is the coordinate in-bounds?
                    if bounds.contains(&c) {
                        line0.push("    ");
                        line1.push("    ");
                    } else {
                        line0.push("████");
                        line1.push("████");
                    }
                }
            }
        }
        // Output the two lines
        println!("{}", line0.concat());
        println!("{}", line1.concat());
    }
}

/// Auxiliary data structures.
/// In these, we store pre-processing data, speeding
/// up actual computation.
#[cfg_attr(test, derive(PartialEq, Debug))]
struct Auxiliaries {
    /// The `Bounds` of the board
    bounds: Bounds,
    /// The puzzle's `Shapekey`
    shapekey: Shapekey,
    /// The puzzle's initial `Blockstate`
    start_blockstate: Blockstate,
    /// The puzzle's `Nonintersectionkey`
    nonintersectionkey: Nonintersectionkey,
    /// A goalblock `Blockstate.goal_offsets[ix]` has
    /// shape `shapekey[goal_shapekey_key[ix]]`.
    goal_shapekey_key: GoalShapekeyKey,
    /// A goalblock `Blockstate.goal_offsets[ix]` has
    /// target-offset `goal_target_offsets[ix]`.
    goal_target_offsets: GoalTargetOffsets,
    /// The puzzle's `ReconstructionMap`, used for printing
    reconstruction_map: ReconstructionMap,
    /// The puzzle's `Width`
    width: Width,
    /// The puzzle's `Height`
    height: Height,
}

/// Possible non-error results of pre-processing
#[cfg_attr(test, derive(PartialEq, Debug))]
enum PreprocessingOutput {
    /// The puzzle is empty
    EmptyPuzzle,
    /// The puzzle is not empty, here have its `Auxiliaries`
    ProperPuzzle(Auxiliaries),
}

/// Preprocess a puzzle, given by string-representation of
/// `start` and `goal`.
fn preprocessing(start: &str, goal: &str) -> Result<PreprocessingOutput, SolvePuzzleError> {
    /// Build a `Nonintersectionkey` from a `Bounds` and `Shapekey`.
    fn build_nonintersectionkey(
        bounds: &Bounds,
        shapekey: &Shapekey,
        width: Width,
        height: Height,
    ) -> Nonintersectionkey {
        /// Helper-struct to represent a relative offset,
        /// because I'd like to have negative values here.
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        struct RelativeOffset(isize, isize); // (isize, isize) because Coor is too small and I want negatives

        let mut relative_nik: HashMap<(usize, usize, RelativeOffset), bool> = Default::default();
        let mut inbounds: HashMap<(usize, Offset), bool> = Default::default();

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
                                                let shifted_b: Points =
                                                    shift_points(shape_b, shift_b);
                                                shifted_b.is_disjoint(&shifted_a)
                                            })
                                            && *inbounds.entry((shape_ix_b, shift_b)).or_insert({
                                                let shifted_b: Points =
                                                    shift_points(shape_b, shift_b);
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

    /// Preprocess a `PreprocessingOutput::ProperPuzzle`.
    /// Assumes that `bounds` exist in `CharToPoints`, i.e. that
    /// the puzzle is not empty
    fn preprocess_proper_puzzle(
        start_chartopoints: &CharToPoints,
        goal_chartopoints: &CharToPoints,
        width: Width,
        height: Height,
    ) -> Result<PreprocessingOutput, SolvePuzzleError> {
        // Check that bounds match.
        let start_bounds = &start_chartopoints[&BOUNDS_CHAR];
        let goal_bounds = &goal_chartopoints[&BOUNDS_CHAR];
        if start_bounds != goal_bounds {
            return Err(SolvePuzzleError::MismatchedBounds);
        }

        let bounds: Shape = start_chartopoints[&BOUNDS_CHAR].clone();

        // Keep track of which char is which shape.
        // This is only used to map goal-shapes to their indices in shapekey later
        let mut char_to_shape: HashMap<char, Shape> = Default::default();
        // Keep track of which shape has which chars and offsets
        let mut shape_to_chars_and_offsets: HashMap<Shape, Vec<(char, Offset)>> =
            Default::default();
        // Keep track of which shape-offset-pair has which char
        let mut reconstruction_map: ReconstructionMap = Default::default();
        // Keep track of which goal-char has which start-offset and target-offset
        let mut goal_chars_startoffset_targetoffset: Vec<(char, Offset, Offset)> = Vec::new();
        for (c, start_points) in start_chartopoints {
            if c == &BOUNDS_CHAR {
                continue;
            }

            // Extract the shape
            let (shape_min, _) = get_points_dimensions(start_points);
            let shape: Shape = start_points
                .iter()
                .map(|point| *point - shape_min)
                .collect();

            // Insert into `char_to_shape`
            char_to_shape.insert(*c, shape.clone());
            // Insert into `shape_to_chars_and_offsets`
            shape_to_chars_and_offsets
                .entry(shape.clone())
                .or_default()
                .push((*c, shape_min));
            // Insert into `reconstruction_map`
            reconstruction_map.insert((shape, shape_min), *c);

            // Is this also a goal-block?
            if let Some(goal_points) = goal_chartopoints.get(c) {
                let (target_min, _) = get_points_dimensions(goal_points);
                // If so, insert into `goal_chars_startoffset_targetoffset`
                goal_chars_startoffset_targetoffset.push((*c, shape_min, target_min));
            }
        }
        // Keep a separate `char_to_goalshape` map
        let mut char_to_goalshape: HashMap<char, Shape> = Default::default();
        for (c, goal_points) in goal_chartopoints {
            // We already know the bounds match, so we don't need to care about those
            if c == &BOUNDS_CHAR {
                continue;
            }
            let (shape_min, _) = get_points_dimensions(goal_points);
            let shape: Shape = goal_points.iter().map(|point| *point - shape_min).collect();
            char_to_goalshape.insert(*c, shape.clone());
        }

        // Check that the start and goal shapes are the same
        for (c, goalshape) in &char_to_goalshape {
            let startshape = char_to_shape
                .get(c)
                .ok_or(SolvePuzzleError::GoalblockWithoutStartingblock(*c))?;
            if startshape != goalshape {
                return Err(SolvePuzzleError::MismatchedGoalShapes(*c));
            }
        }

        // Extract the key-value-pairs from `shape_to_chars_and_offsets` into
        // a vector. We'll be sorting them in a second
        let mut raw_shapekey: Vec<(Shape, Vec<(char, Offset)>)> = shape_to_chars_and_offsets
            .iter()
            .map(|(shape, chars_and_offsets)| (shape.clone(), chars_and_offsets.clone()))
            .collect();

        // Sort these values to make intersection-finding terminate earlier, hopefully.
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
                    // If both or neither are only for goals, sort by size first, idea
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

        // Create shapekey from `raw_shapekey`
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

        // Collect the nongoal-offsets
        let mut nongoal_offsets: NongoalOffsets = raw_shapekey
            .iter()
            .map(|(_, chars_and_offsets)| {
                chars_and_offsets
                    .iter()
                    .filter(|(c, _)| !goal_chartopoints.contains_key(c))
                    .map(|(_, offset)| *offset)
                    .collect()
            })
            .filter(|offsets: &Offsets| !offsets.is_empty())
            .collect();
        nongoal_offsets.shrink_to_fit();

        // Collect the goal-blocks' initial-offsets
        let mut goal_offsets: GoalOffsets = goal_chars_startoffset_targetoffset
            .iter()
            .map(|(_, start, _)| *start)
            .collect();
        goal_offsets.shrink_to_fit();

        // Collect the goal-blocks' target-offsets
        let mut goal_target_offsets: GoalTargetOffsets = goal_chars_startoffset_targetoffset
            .iter()
            .map(|(_, _, target)| *target)
            .collect();
        goal_target_offsets.shrink_to_fit();

        // Build the initial `Blockstate`
        let blockstate: Blockstate = Blockstate {
            nongoal_offsets,
            goal_offsets,
        };

        // Build the nonintersection-key
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

    /// Possible results of `string_to_chartopoints`
    enum StringToCharToPointsResult {
        /// The puzzle is empty
        EmptyPuzzle,
        /// The puzzle is nonempty, here have its `CharToPoints`
        /// and dymensions
        ProperPuzzle(CharToPoints, Width, Height),
    }
    /// Convert a string to a `CharToPoints`-map, which
    /// maps a character to the set of points it occupies.
    fn string_to_chartopoints(s: &str) -> Result<StringToCharToPointsResult, SolvePuzzleError> {
        // Extract bounds-extremes
        let mut min_x = usize::MAX;
        let mut min_y = usize::MAX;
        let mut max_x = usize::MIN;
        let mut max_y = usize::MIN;

        // Annotate each coordinate with the char it represents
        let mut char_annotated_coordinates: Vec<(char, (usize, usize))> = Vec::new();

        for (y, l) in s.lines().enumerate() {
            for (x, c) in l.chars().enumerate() {
                // Is this out-of-bounds?
                if !c.is_whitespace() {
                    min_x = min(min_x, x);
                    min_y = min(min_y, y);
                    max_x = max(max_x, x);
                    max_y = max(max_y, y);
                    char_annotated_coordinates.push((c, (x, y)));
                    // We want this to be in-bounds, as well, but we
                    // don't want to insert it twice if it's already bounds.
                    if c != BOUNDS_CHAR {
                        char_annotated_coordinates.push((BOUNDS_CHAR, (x, y)));
                    }
                }
            }
        }
        // Is this puzzle empty?
        if char_annotated_coordinates.is_empty() {
            return Ok(StringToCharToPointsResult::EmptyPuzzle);
        }

        // Get the dimensions
        let width_usize: usize = max_x + 1 - min_x;
        let height_usize: usize = max_y + 1 - min_y;

        // Try converting to Coor, and otherwise throw SolvePuzzleError::WidthTooLarge
        let width: Width = width_usize
            .try_into()
            .map_err(|_| SolvePuzzleError::WidthTooLarge)?;
        let height: Height = height_usize
            .try_into()
            .map_err(|_| SolvePuzzleError::HeightTooLarge)?;

        // Build the map from char to points by iterating
        // over the annotated coordinates
        let mut char_to_points: CharToPoints = CharToPoints::new();
        for (c, pair) in char_annotated_coordinates {
            let (x_usize, y_usize) = (pair.0 + 1 - min_x, pair.1 + 1 - min_y);
            let x: Coor = x_usize
                .try_into()
                .map_err(|_| SolvePuzzleError::WidthTooLarge)?;
            let y: Coor = y_usize
                .try_into()
                .map_err(|_| SolvePuzzleError::HeightTooLarge)?;
            let point = Offset::new(x, y);
            char_to_points.entry(c).or_default().insert(point);
        }

        // Return the result
        Ok(StringToCharToPointsResult::ProperPuzzle(
            char_to_points,
            width,
            height,
        ))
    }

    // Convert start- and goal-strings to `CharToPoints` maps.
    let start_stc_result = string_to_chartopoints(start)?;
    let goal_stc_result = string_to_chartopoints(goal)?;

    match (start_stc_result, goal_stc_result) {
        // Are both start- and goal-config empty?
        (StringToCharToPointsResult::EmptyPuzzle, StringToCharToPointsResult::EmptyPuzzle) => {
            Ok(PreprocessingOutput::EmptyPuzzle)
        }
        // Are both start- and goal-config nonempty?
        (
            StringToCharToPointsResult::ProperPuzzle(start_chartopoints, width, height),
            StringToCharToPointsResult::ProperPuzzle(goal_chartopoints, _goal_width, _goal_height),
        ) => preprocess_proper_puzzle(&start_chartopoints, &goal_chartopoints, width, height),
        // Is exactly one of start- and goal-config nonempty?
        _ => Err(SolvePuzzleError::MismatchedBounds),
    }
}

/// A heuristic for how many moves it'll take to reach the target-config.
/// This simply counts how many goal-blocks haven't reached their
/// target-offsets yet.
fn misplaced_goalblocks_heuristic(
    blockstate: &Blockstate,
    goal_target_offsets: &GoalTargetOffsets,
) -> usize {
    // You might think that, rather than just counting the misplaced blocks,
    // we could also check their distance from their goal_target_offsets.
    // I implemented and benchmarked that, and it was worse.

    blockstate
        .goal_offsets
        .iter()
        .zip(goal_target_offsets)
        .filter(|(goal_offset, target_offset)| goal_offset != target_offset)
        .count()
}

/// Given `Auxiliaries`, compute the puzzle's solution.
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
        // If there are no goal-blocks, the solution is trivial.
        0 => Some(vec![start_blockstate.clone()]),
        1 => {
            if start_blockstate.nongoal_offsets.is_empty() {
                // If this goalblock is also the only block in the puzzle, we can't just call
                // the usual BFS function, as that assumes we have at least two distinct
                // blocks in the puzzle to do its in-bounds-check.
                // However, now we can just move this one goalblock around using
                // BFS, and check if it can reach its target-position.

                if start_blockstate.goal_offsets == *goal_target_offsets {
                    return Some(vec![start_blockstate.clone()]);
                }
                let beginning_offset = start_blockstate.goal_offsets[0];
                let is_legal = |offset: Offset| -> bool {
                    nonintersectionkey.abuse_this_datastructure_for_in_bounds_check(0, offset)
                };
                let neighbors = dfs_general(beginning_offset, &is_legal).collect_vec();
                // Does it contain the target-position?
                neighbors.contains(&goal_target_offsets[0]).then(|| {
                    // Return path of length 1 (having two nodes)
                    let mut goal_blockstate = start_blockstate.clone();
                    goal_blockstate.goal_offsets.clone_from(goal_target_offsets);
                    vec![start_blockstate.clone(), goal_blockstate]
                })
            } else {
                // There is at least one other block, so we can use the
                // BFS-function that exploits other blocks for in-bounds-checking.
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
                        // We're successful once all goal-blocks reach their target-offset
                        blockstate.goal_offsets == *goal_target_offsets
                    },
                )
                .map(|path| {
                    // Extract the blockstates from the path
                    path.into_iter()
                        .map(|BlockstateJustmoved { blockstate, .. }| blockstate)
                        .collect_vec()
                })
            }
        }
        // If there is more than one goalblock, we can use A*
        _ => pathfinding::directed::astar::astar(
            &BlockstateJustmoved {
                blockstate: start_blockstate.clone(),
                justmoved: Justmoved::Nothing,
            },
            |blockstate| {
                get_neighboring_blockstates(blockstate, nonintersectionkey, goal_shapekey_key)
                    .into_iter()
                    .map(|bst| (bst, 1))
            },
            // Use the [`misplaced_goalblocks_heuristic`]
            |BlockstateJustmoved { blockstate, .. }| {
                misplaced_goalblocks_heuristic(blockstate, goal_target_offsets)
            },
            // We're successful once all goal-blocks reach their target-offset
            |BlockstateJustmoved { blockstate, .. }| {
                blockstate.goal_offsets == *goal_target_offsets
            },
        )
        .map(|path| {
            // Extract the blockstates from the path
            path.0
                .into_iter()
                .map(|BlockstateJustmoved { blockstate, .. }| blockstate)
                .collect_vec()
        }),
    }
}

/// Solve a puzzle, returning its path
fn solve_puzzle_path(start: &str, goal: &str) -> Result<Option<Vec<Blockstate>>, SolvePuzzleError> {
    match preprocessing(start, goal)? {
        PreprocessingOutput::ProperPuzzle(auxiliaries) => {
            Ok(solution_from_auxiliaries(&auxiliaries))
        }
        PreprocessingOutput::EmptyPuzzle => Err(SolvePuzzleError::EmptyPuzzle),
    }
}

/// Solve a puzzle, returning its minimum number of moves.
/// See [`solve_puzzle`] for example usage.
#[inline]
pub fn solve_puzzle_minmoves(start: &str, goal: &str) -> Result<Option<usize>, SolvePuzzleError> {
    let maybe_path = solve_puzzle_path(start, goal)?;
    Ok(maybe_path.map(|path| path.len() - 1))
}

/// Solve a sliding-blocks-puzzle from its string-representation.
///
/// Returns:
/// - [`SolvePuzzleError`] if the puzzle is malformed
/// - `Ok(None)` if the puzzle is well-formed, but not solvable
/// - `Ok(Some(Path))` if the puzzle is well-formed and solvable,
///     where `Path` is a `Vec` of string-representations of boardstates,
///     similar to the input-format.
///
/// The representation is a visually intuitive newline-separated string, which makes up
/// a grid. For instance:
/// ```txt
/// abc
/// c .
/// ```
/// would be a grid of width `3` and height `2`.
/// - Whitespace are cells that are out-of-bounds.
/// - Periods `.` are cells that are in-bounds, but not occupied by any block.
/// - Any other character represents a block that is occupying that cell. Such
///     a cell is automatically in-bounds.
///
/// A block is identified by its character, and occupies all the cells that
/// are occupied by that character. Blocks can consist of multiple grid-cells
/// (by notation, they always consist of at least one cell). No other topological
/// assumptions are made -- blocks need not be convex, they need not be
/// simply connected, and they need not even be connected. For example, the above
/// puzzle has three distinct blocks (`a` and `b` occupy one cell each, `c` occupies
/// two cells) and five in-bounds cells, only one of which is not occupied.
///
/// The goal-configuration must have the same bounds as the start-configuration. A block
/// is considered a goal-block if it is present in the goal-configuration. Not every
/// block needs to be a goal-block. A goal-block must have the same shape in the
/// goal-configuration as in the start-configuration. See the [`examples`]-crate for
/// example-puzzles. There may be an arbitrary number of goal-blocks and non-goal-blocks.
///
/// # Examples
/// ```
/// use sliding_blocks::solve_puzzle;
///
/// let start = "
///   A.b
///    .
/// ";
/// let end = "
///   ..A
///    .
/// ";
/// let computed = solve_puzzle(start, end);
///
/// // A.b  ->  A..  ->  ..A
/// //  .   ->   b   ->   b  
/// let path = vec![
///     "A.b\n . ".to_owned(),
///     "A..\n b ".to_owned(),
///     "..A\n b ".to_owned(),
/// ];
/// assert_eq!(computed, Ok(Some(path)));
/// ```
/// The shortest path is unique in this example, but this need not always be the case.
/// Goal-blocks need not be upper-case, non-goal-blocks need not be lower-case.
///
/// See the [`examples`] crate for more example-puzzles.
#[inline]
#[wasm_bindgen]
pub fn solve_puzzle(start: &str, goal: &str) -> Result<Option<Vec<String>>, SolvePuzzleError> {
    match preprocessing(start, goal)? {
        PreprocessingOutput::EmptyPuzzle => Err(SolvePuzzleError::EmptyPuzzle),
        PreprocessingOutput::ProperPuzzle(auxiliaries) => {
            Ok(solution_from_auxiliaries(&auxiliaries).map(|path| {
                // Now that we have a path, we need to reconstruct it to a string.

                // Get shorter auxiliaries-names
                let width = auxiliaries.width;
                let height = auxiliaries.height;
                let bounds = auxiliaries.bounds;
                let shapekey = auxiliaries.shapekey;
                let goal_shapekey_key = auxiliaries.goal_shapekey_key;

                // Output a reconstruction-map to a string:
                let reconstruction_map_to_string = |rm: &ReconstructionMap| -> String {
                    // Initialise the board as a grid of spaces
                    let mut board: Vec<Vec<char>> =
                        vec![vec![' '; width as usize]; height as usize];

                    // Insert the bounds
                    for point in &bounds {
                        board[point.y_usize() - 1][point.x_usize() - 1] = BOUNDS_CHAR;
                    }

                    // Insert all the other blocks
                    for ((shape, offset), c) in rm {
                        for twoffset in shape {
                            board[twoffset.y_usize() + offset.y_usize() - 1]
                                [twoffset.x_usize() + offset.x_usize() - 1] = *c;
                        }
                    }

                    // Convert it into a string using linebreaks
                    board
                        .iter()
                        .map(|row| row.iter().collect::<String>())
                        .collect::<Vec<_>>()
                        .join("\n")
                };
                // We'll mutate `reconstruction_map` along the way
                let mut reconstruction_map = auxiliaries.reconstruction_map;
                // We'll return `string_path` in the end
                let mut string_path = vec![reconstruction_map_to_string(&reconstruction_map)];
                // Remembmer the previous_blockstate, to determine what block was moved
                let mut previous_blockstate = &path[0];

                // Now walk along the path:
                for blockstate in path.iter().skip(1) {
                    // Figure out which block was moved.
                    'check_single_block: {
                        // Maybe it was a non-goal-block?
                        for (shape_ix, (offsets, previous_offsets)) in blockstate
                            .nongoal_offsets
                            .iter()
                            .zip(previous_blockstate.nongoal_offsets.iter())
                            .enumerate()
                        {
                            if *offsets != *previous_offsets {
                                // It was a goal-block of this shape!

                                // But what offset did it have?
                                let old_offset = previous_offsets
                                    .iter()
                                    .find(|o| !offsets.contains(o))
                                    .unwrap();
                                let new_offset = offsets
                                    .iter()
                                    .find(|o| !previous_offsets.contains(o))
                                    .unwrap();
                                let shape = &shapekey[shape_ix];

                                // Remove old char from `reconstruction_map`
                                let c = reconstruction_map
                                    .remove(&(shape.clone(), *old_offset))
                                    .unwrap();
                                // Insert new char into `reconstruction_map`
                                reconstruction_map.insert((shape.clone(), *new_offset), c);

                                break 'check_single_block;
                            }
                        }

                        // Or maybe it was a goal-block?
                        for (goalvec_ix, (offset, previous_offset)) in blockstate
                            .goal_offsets
                            .iter()
                            .zip(previous_blockstate.goal_offsets.iter())
                            .enumerate()
                        {
                            if *offset != *previous_offset {
                                // Yep, it was this one:
                                let old_offset = previous_offset;
                                let new_offset = offset;
                                let shape = &shapekey[goal_shapekey_key[goalvec_ix]];

                                // Remove char from `reconstruction_map`
                                let c = reconstruction_map
                                    .remove(&(shape.clone(), *old_offset))
                                    .unwrap();
                                // Insert it again with new shape
                                reconstruction_map.insert((shape.clone(), *new_offset), c);

                                break 'check_single_block;
                            }
                        }
                    }
                    // Push the stringified `reconstruction_map` to `string_path`
                    string_path.push(reconstruction_map_to_string(&reconstruction_map));
                    // Remember this `blockstate` for the next move
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

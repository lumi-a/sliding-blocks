pub mod examples;

use colored::{self, Colorize};
use itertools::{iproduct, Itertools};
use ordered_float::OrderedFloat as FloatOrd;
use std::cmp::{max, min, Ordering};
use std::collections::hash_map::Entry::Occupied;
use std::collections::hash_map::Entry::Vacant;
use std::collections::{BTreeSet, BinaryHeap, VecDeque};
use std::collections::{HashMap, HashSet};

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

type Floaty = FloatOrd<f32>; // I *need* a total order on floats, *please*

type MinkowskiDiams = Vec<Vec<Floaty>>; // Given an index in the Blockstate.goal_blocks vec, and an index of its shape in the shapekey vec, what is the diameter of the minkowski sum of the two shapes? Used for heuristic.

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
            .or_insert_with(Vec::new)
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
                    a_size.cmp(&b_size)
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
                    .filter(|(c, _)| !goal_chartopoints.get(c).is_some())
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

fn get_minkowski_diams(goal_shapekey_key: &GoalShapekeyKey, shapekey: &Shapekey) -> MinkowskiDiams {
    let mut minkowski_diams = MinkowskiDiams::new();
    for goalvec_ix in goal_shapekey_key.iter() {
        let goal_shape: BTreeSet<(isize, isize)> = shapekey[*goalvec_ix]
            .iter()
            .map(|Point(x, y)| (*x as isize, *y as isize))
            .collect();
        let goal_max_x = *goal_shape.iter().map(|(x, _)| x).max().unwrap(); // TODO: Two iterations, kinda inefficient
        let goal_max_y = *goal_shape.iter().map(|(_, y)| y).max().unwrap(); // We could finally extract this into a function?
        let mut goal_diams: Vec<Floaty> = Vec::new();
        for shape in shapekey {
            let shape: BTreeSet<(isize, isize)> = shape
                .iter()
                .map(|Point(x, y)| (*x as isize, *y as isize))
                .collect();
            let shape_max_x = *shape.iter().map(|(x, _)| x).max().unwrap(); // TODO: Two iterations, kinda inefficient
            let shape_max_y = *shape.iter().map(|(_, y)| y).max().unwrap(); // We could finally extract this into a function?
            
            // TODO: NONE of this actually is a minkowski-sum. Rename all variables
            let mut minkowski_sum: BTreeSet<(isize, isize)> = BTreeSet::new();
            // The loop-variables dx, dy denote the offset of goal_shape.
            // We'll add (goal_shape + offset) to minkowski_sum iff
            // (goal_shape + offset) intersects shape. For that, looping over [(-goal_max_x)..=shape_max_x]
            // suffices (as oppposed to (-width..=width) or something), which I hope you can
            // visualise in your head.
            // TODO: Could we actually loop over an even smaller set?
            for dx in (-goal_max_x)..=shape_max_x {
                for dy in (-goal_max_y)..=shape_max_y {
                    let moved_goal_shape: BTreeSet<(isize, isize)> =
                        goal_shape.iter().map(|(x, y)| (x + dx, y + dy)).collect();
                    if moved_goal_shape.intersection(&shape).count() > 0 {
                        // TODO: How inefficient in this?
                        minkowski_sum = minkowski_sum.union(&moved_goal_shape).cloned().collect();
                    }
                }
            }

            // Note that the maximum taxicab-distance is NOT admissible: Consider a U-shape, and a single-cell-shape.
            // The actual diameter should be 5 (if the U-shape has volume 5), but the taxicab-maximum is 3.
            // Instead, we'll calculate diameter of each of the connected components of the graph of the sum of the two shapes
            // (where two cells are adjacent iff they're adjacent on the grid) and add 1.
            let nodes: Vec<(isize, isize)> = minkowski_sum.iter().cloned().collect();
            let connected_components =
                pathfinding::undirected::connected_components::connected_components(
                    &nodes,
                    |(x, y)| -> Vec<(isize, isize)> {
                        [(*x - 1, *y), (*x + 1, *y), (*x, *y - 1), (*x, *y + 1)]
                            .iter()
                            .filter(|a| minkowski_sum.contains(a))
                            .cloned()
                            .collect_vec()
                    },
                );

            let diam: usize = connected_components
                .iter()
                .map(|component| {
                    // find diameter of component
                    iproduct!(component, component)
                        .filter(|(v, w)| v > w) // Halves processing
                        .map(|(v, w)| {
                            pathfinding::directed::bfs::bfs(
                                v,
                                |(x, y)| -> Vec<(isize, isize)> {
                                    [(*x - 1, *y), (*x + 1, *y), (*x, *y - 1), (*x, *y + 1)]
                                        .iter()
                                        .filter(|a| component.contains(a))
                                        .cloned()
                                        .collect_vec()
                                },
                                |u| u == w,
                            )
                            .unwrap()
                            .len()
                                - 2 // Subtraction will never overflow, because v!=w
                        })
                        .max()
                        .unwrap_or(1)
                })
                .sum();

            goal_diams.push(FloatOrd(1.0 / diam as f32));
        }
        minkowski_diams.push(goal_diams)
    }
    minkowski_diams
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

fn heuristic(
    blockstate: &Blockstate,
    nonintersectionkey: &Nonintersectionkey,
    goal_shapekey_key: &GoalShapekeyKey,
    goal_target_offsets: &GoalTargetOffsets,
    minkowski_diams: &MinkowskiDiams,
    width: Coor,
    height: Coor,
) -> Floaty {
    // TODO: Prove this is admissible
    let nongoal_offsets = &blockstate.nongoal_offsets;
    let goal_offsets = &blockstate.goal_offsets;
    goal_offsets
        .iter()
        .enumerate()
        .map(|(goalvec_ix, current_goal_offset)| {
            let goal_shapeix = goal_shapekey_key[goalvec_ix];
            let nik_goal = &nonintersectionkey[goal_shapeix];
            let minkowskis_goal = &minkowski_diams[goalvec_ix];
            let goal_target_offset = &goal_target_offsets[goalvec_ix];

            // TODO: This dijkstra-search can be improved: Since only the vertices carry costs, we actually
            // calculate them too often (since we calculate them again for each _edge_).
            // I'm afraid you'll have to use another implementation of dijkstra for that, sorry.
            let (_, cost): (_, Floaty) = pathfinding::directed::dijkstra::dijkstra(
                current_goal_offset,
                |goal_offset| {
                    [goal_offset.up(), goal_offset.down(), goal_offset.left(), goal_offset.right()]
                    .iter()
                    .filter(|offset| {
                        // TODO: This is a TERRIBLE hack to check whether offset is in-bounds
                        let x = offset.0;
                        let y = offset.1;
                        nik_goal[x as usize - 1][y as usize - 1].is_empty()
                    })
                    .map(|offset| -> (Offset, Floaty) {
                        (
                            offset.clone(), // new offset
                            {
                                let nongoal_sum: Floaty =
                                    nongoal_offsets // cost
                                        .iter()
                                        .enumerate()
                                        .map(|(shape_ix, offsets)| -> Floaty {
                                            FloatOrd(
                                                offsets
                                                    .iter()
                                                    .filter(|Offset(x, y)| {
                                                        !nik_goal[offset.0 as usize - 1]
                                                            [offset.1 as usize - 1][shape_ix]
                                                            [*x as usize]
                                                            [*y as usize]
                                                    })
                                                    .count()
                                                    as f32,
                                            ) * minkowskis_goal[shape_ix]
                                        })
                                        .sum();

                                let goal_sum: Floaty = goal_offsets
                                    .iter()
                                    .enumerate()
                                    .filter(|(this_goal_ix, Offset(x, y))| {
                                        !nik_goal[offset.0 as usize - 1][offset.1 as usize - 1]
                                            [goal_shapekey_key[*this_goal_ix]]
                                            [*x as usize][*y as usize]
                                            && *this_goal_ix != goalvec_ix
                                    })
                                    .map(|(this_goal_ix, _)| -> Floaty {
                                        minkowskis_goal[goal_shapekey_key[this_goal_ix]]
                                    })
                                    .sum();

                                nongoal_sum + goal_sum
                            },
                        )
                    }) // TODO: implement actual sum
                    .collect_vec()
                },
                |goal_offset| *goal_offset == *goal_target_offset,
            )
            .unwrap(); // TODO: This unwrap might very well fail in practice
            cost
        })
        .max()
        .unwrap() // TODO: Bad if we have no goal blocks
}

fn tighter_heuristic(
    blockstate: &Blockstate,
    nonintersectionkey: &Nonintersectionkey,
    goal_shapekey_key: &GoalShapekeyKey,
    goal_target_offsets: &GoalTargetOffsets,
    width: Coor,
    height: Coor,
) -> Floaty {
    let nongoal_offset_vec: Vec<Offset> = blockstate
        .nongoal_offsets
        .iter()
        .flatten()
        .cloned()
        .collect();
    let goal_offset_vec = blockstate.goal_offsets;

    #[derive(PartialEq, Eq)]
    struct Hypervertex {
        pulverized_blocks: BTreeSet<usize>, // i ∈ pulverized_blocks <=> block (nongoal_offset_vec:goal_offset_vec)[i] is translucent
        offsets: Offsets, // Offsets that are possible to reach while pulverized_blocks are translucent
        fringe: Vec<Offset>, // Subset of offsets. Contains an offset iff it has in-bounds neighbour-offset that can't be accessed without pulverizing more blocks
    }

    let pulverized_bfs = |goalvec_ix: usize| -> usize {
        let sub_bfs = |&Hypervertex{ pulverized_blocks, offsets, fringe}, newly_pulverized: usize| -> Vec<Hypervertex> {
            let mut new_pulverized_blocks = pulverized_blocks.clone();
            new_pulverized_blocks.insert(newly_pulverized);
            let mut seen_offsets: Offsets = offsets.clone();
            let mut stack: Vec<Offset> = fringe.clone();
            let mut new_fringe: Vec<Offset> = Vec::new();
        }
        let target_goal_offset = goal_target_offsets[goalvec_ix];
        let mut seen_offsets: BTreeSet<Offset> = BTreeSet::new();
        let mut stack: Vec<OffsetAndPulverizedBlocks> = vec![OffsetAndPulverizedBlocks {
            offset: goal_offset_vec[goalvec_ix],
            pulverized_blocks: BTreeSet::new(),
        }];
        while let Some(OffsetAndPulverizedBlocks {
            offset,
            pulverized_blocks,
        }) = stack.pop()
        {
            let pulverized_count = pulverized_blocks.len();
            if offset == target_goal_offset {
                return pulverized_count;
            }

            // TODO: This is AWFUL code.
            // And no, we can't just extract the array into match-arms and for-loop
            // over the match-result, because arrays of different length are
            // different types
            // Should we extend the intersection-vector to reach one more cell? Hmm
            let inserty = |new_offset: Offset| {
                if new_offset.0 < width && new_offset.1 < height {
                    panic!("TODO: Iterate over new blocks you may hit now");

                    //&& seen_offsets.insert(new_offset)
                    //&& is_legal(new_offset)
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
            } else if offset.1 > 0 {
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
        legal_offsets
    };
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
    // - To bfs a block, remove it from U, and test its legality by checking
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
    */ // for now, just this bad implementation of bfs:
    // TODO: I have no idea if `&dyn Fn(Offset) -> bool` is the right signature as I didn't learn about `&dyn` yet
    let bfs_general = |initial_offset: Offset, is_legal: &dyn Fn(&Offset) -> bool| {
        let mut legal_offsets: Vec<Offset> = Vec::new();
        let mut seen_offsets: BTreeSet<Offset> = BTreeSet::new(); // TODO: Different data structures?
        seen_offsets.insert(initial_offset.clone());
        let mut stack: Vec<Offset> = vec![initial_offset];
        while let Some(offset) = stack.pop() {
            // TODO: Maybe this can be made faster by not going back in
            // the direction we just came from
            for new_offset in [offset.up(), offset.down(), offset.left(), offset.right()] {
                // TODO: Which of these checks is faster?
                // TODO: Can we check seen_offsets membership without cloning?
                if is_legal(&new_offset) && seen_offsets.insert(new_offset.clone()) {
                    legal_offsets.push(new_offset.clone());
                    stack.push(new_offset);
                }
            }
        }

        legal_offsets
    };
    let bfs_nongoal = |movingshape_ix: usize,
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
        bfs_general(moving_offset, &is_legal)
    };

    let bfs_goal = |moving_goalvec_ix: usize, moving_offset: Offset| -> Vec<Offset> {
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
        bfs_general(moving_offset, &is_legal)
    };

    // It's okay to gather all these into a vector rather than a set,
    // because all neighbors WILL be unique.
    let mut neighboring_blockstates: Vec<Blockstate> = Vec::new();
    // Start with blockstate.goal_offsets first, to hopefully find the goal a little sooner
    for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
        // TODO: avoid .clone() here?
        for mutated_offset in bfs_goal(goalvec_ix, offset.clone()) {
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
            for mutated_offset in bfs_nongoal(shapekey_ix, offset.clone(), &trimmed_shape_offsets) {
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
    MinkowskiDiams,
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
            char_to_points
                .entry(c)
                .or_insert_with(Points::new)
                .insert(shifted_point);
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
    let minkowski_diams = get_minkowski_diams(&goal_shapekey_key, &shapekey);
    let nonintersectionkey = build_nonintersectionkey(&bounds, &shapekey, width, height);

    (
        bounds,
        shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        minkowski_diams,
        width,
        height,
    )
}

// TODO: Implement A*
// TODO: Add tests

pub fn solve_puzzle_bfs(start: &str, goal: &str) {
    let (
        bounds,
        shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        _minkowski_diams,
        width,
        height,
    ) = puzzle_preprocessing(start, goal);
    /*
    print_puzzle(
        &bounds,
        &shapekey,
        &start_blockstate,
        &goal_shapekey_key,
        width,
        height,
    );
    */

    let path = pathfinding::directed::bfs::bfs(
        &start_blockstate,
        |blockstate| {
            get_neighboring_blockstates(&blockstate, &nonintersectionkey, &goal_shapekey_key)
        },
        |blockstate| blockstate.goal_offsets == goal_target_offsets,
    )
    .unwrap();

    assert!(path.len() < 1000); // TODO: Remove this
}

pub fn solve_puzzle_astar(start: &str, goal: &str) {
    let (
        _bounds,
        _shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        minkowski_diams,
        width,
        height,
    ) = puzzle_preprocessing(start, goal);

    let (path, _) = pathfinding::directed::astar::astar(
        &start_blockstate,
        |blockstate| -> Vec<(Blockstate, Floaty)> {
            get_neighboring_blockstates(
                blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
            )
            .iter()
            .map(|blockstate| (blockstate.clone(), FloatOrd(1.0))) // TODO: This is horrible?
            .collect()
        },
        |blockstate| {
            heuristic(
                blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
                &goal_target_offsets,
                &minkowski_diams,
                width,
                height,
            )
        },
        |blockstate| blockstate.goal_offsets == goal_target_offsets,
    )
    .unwrap();

    assert!(path.len() < 1000); // TODO: Remove this
}

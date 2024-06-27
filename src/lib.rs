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
type Coordinates = (Coor, Coor);

type CoordinatesSet = BTreeSet<Coordinates>;
type CharToCoors = HashMap<char, CoordinatesSet>;

type Shape = BTreeSet<Coordinates>;
type Bounds = Shape;
type Offset = (Coor, Coor);

type Offsets = BTreeSet<Offset>;
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Blockstate {
    nongoal_offsets: Vec<Offsets>, // TODO: Perhaps this is better done on the stack, e.g. with https://crates.io/crates/arrayvec
    goal_offsets: Vec<Offset>,
}
type Shapekey = Vec<Shape>;
type GoalShapekeyKey = Vec<usize>; // Given an index in the Blockstate.goal_blocks vec, what is the index of its shape in the shapekey vec?
type GoalTargetOffsets = Vec<Offset>; // At what offset is a block in a goal-position?
type Coortable<T> = Vec<Vec<T>>;

type Floaty = FloatOrd<f32>; // I *need* a total order on floats, *please*

type MinkowskiDiams = Vec<Vec<Floaty>>; // Given an index in the Blockstate.goal_blocks vec, and an index of its shape in the shapekey vec, what is the diameter of the minkowski sum of the two shapes? Used for heuristic.

fn _intersect_coortables(a: &Coortable<bool>, b: &Coortable<bool>) -> Coortable<bool> {
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
        chartocoors.entry(c).or_default().insert(shifted_coor);
    }

    (chartocoors, width, height)
}

// TODO: This function has a bad name by now. I think it can also be moved into
// the pre-processing function?
fn extract_shapekey(
    start_chartocoors: &CharToCoors,
    goal_chartocoors: &CharToCoors,
) -> (
    Bounds,
    Shapekey,
    Blockstate,
    GoalShapekeyKey,
    GoalTargetOffsets,
) {
    let mut chartoshape: HashMap<char, Shape> = HashMap::new();
    let mut shape_to_chars_and_offsets: HashMap<Shape, Vec<(char, Offset)>> = HashMap::new();
    let mut goal_chars_startoff_targetoff: Vec<(char, Offset, Offset)> = Vec::new(); // (char, start-offset, target-offset)

    // TODO: Handle empty strings gracefully
    // TODO: Also maybe don't use clone, but I'm not into ownership enough to think through how to handle this
    let bounds: Shape = start_chartocoors.get(&BOUNDS_CHAR).unwrap().clone();

    fn get_mins(coordinates_set: &CoordinatesSet) -> (Coor, Coor) {
        // Extract min-x and min-y.
        // Assumes that coordinatesSet is nonempty.
        // TODO: Is that a misassumption?
        let mut min_x = Coor::MAX;
        let mut min_y = Coor::MAX;
        for coor in coordinates_set {
            min_x = min(min_x, coor.0);
            min_y = min(min_y, coor.1);
        }
        (min_x, min_y)
    }

    for (c, start_coords) in start_chartocoors.iter() {
        if c == &BOUNDS_CHAR {
            continue;
        }

        let (min_x, min_y) = get_mins(start_coords);
        let shape: Shape = start_coords
            .iter()
            .map(|coor| (coor.0 - min_x, coor.1 - min_y))
            .collect();

        // This is only used to map goal-shapes to their indices in shapekey later
        chartoshape.insert(*c, shape.clone());

        shape_to_chars_and_offsets
            .entry(shape)
            .or_default()
            .push((*c, (min_x, min_y)));
        if let Some(goal_coords) = goal_chartocoors.get(c) {
            // TODO: Because we're currently in the start-loop, and won't separately
            // do a goal-loop, we'll silently ignore any characters in the goalstring
            // that aren't present in the startstring. If there are chars in the goalstring
            // that are not in the startstring, then that puzzle is not correctly posed,
            // which we should communicate instead.

            let (target_min_x, target_min_y) = get_mins(goal_coords);
            goal_chars_startoff_targetoff.push((*c, (min_x, min_y), (target_min_x, target_min_y)));
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
                .all(|(c, _)| goal_chartocoors.get(c).is_some());
            let b_shape_only_for_goals = b_chars_and_offsets
                .iter()
                .all(|(c, _)| goal_chartocoors.get(c).is_some());

            // Sort shapes that are only for goals last
            if a_shape_only_for_goals && !b_shape_only_for_goals {
                Ordering::Greater
            } else if !a_shape_only_for_goals && b_shape_only_for_goals {
                Ordering::Less
            } else {
                // If both or neither are only for goals, sort by size first
                // (Idea being: If shapes are larger, we find intersections earlier)
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
    let goal_shapekey_key: GoalShapekeyKey = goal_chars_startoff_targetoff
        .iter()
        .map(|(c, _, _)| {
            raw_shapekey
                .iter()
                .position(|(shape, _)| shape == chartoshape.get(c).unwrap())
                .unwrap()
        })
        .collect();

    let blockstate: Blockstate = Blockstate {
        nongoal_offsets: raw_shapekey
            .iter()
            .map(|(_, chars_and_offsets)| {
                chars_and_offsets
                    .iter()
                    .filter(|(c, _)| goal_chartocoors.get(c).is_none())
                    .map(|(_, offset)| *offset)
                    .collect()
            })
            .filter(|offsets: &Offsets| !offsets.is_empty())
            .collect(),
        goal_offsets: goal_chars_startoff_targetoff
            .iter()
            .map(|(_, start, _)| *start)
            .collect(),
    };
    let goal_target_offsets = goal_chars_startoff_targetoff
        .iter()
        .map(|(_, _, target)| *target)
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
            .map(|(x, y)| (*x as isize, *y as isize))
            .collect();
        let goal_max_x = *goal_shape.iter().map(|(x, _)| x).max().unwrap(); // TODO: Two iterations, kinda inefficient
        let goal_max_y = *goal_shape.iter().map(|(_, y)| y).max().unwrap(); // We could finally extract this into a function?
        let mut goal_diams: Vec<Floaty> = Vec::new();
        for shape in shapekey {
            let shape: BTreeSet<(isize, isize)> = shape
                .iter()
                .map(|(x, y)| (*x as isize, *y as isize))
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
            let goal_target_offset = goal_target_offsets[goalvec_ix];
            let goal_target_offset_hack = (goal_target_offset.0 + 1, goal_target_offset.1 + 1);

            // TODO: This dijkstra-search can be improved: Since only the vertices carry costs, we actually
            // calculate them too often (since we calculate them again for each _edge_).
            // I'm afraid you'll have to use another implementation of dijkstra for that, sorry.
            let (_, cost): (_, Floaty) = pathfinding::directed::dijkstra::dijkstra(
                // TODO: This is a TERRIBLE hack that's only necessary because
                // I didn't implement the "hey bounds and all coordinates should have a buffer of 1 to
                // avoid checking against over- or underflows" fix yet:
                // We shift everything by 1. I'm sorry.
                // If you undo this, also undo goal_target_offset_hack to goal_target_offset
                &(current_goal_offset.0 + 1, current_goal_offset.1 + 1),
                |goal_offset| {
                    [
                        (goal_offset.0, goal_offset.1 + 1),
                        (goal_offset.0 + 1, goal_offset.1),
                        (goal_offset.0, goal_offset.1 - 1),
                        (goal_offset.0 - 1, goal_offset.1),
                    ]
                    .iter()
                    .filter(|offset| {
                        // TODO: This is a TERRIBLE hack to check whether offset is in-bounds
                        let x = offset.0;
                        let y = offset.1;
                        x > 0
                            && y > 0
                            && x <= width
                            && y <= height
                            && !nik_goal[x as usize - 1][y as usize - 1].is_empty()
                    })
                    .map(|offset| -> (Offset, Floaty) {
                        (
                            *offset, // new offset
                            {
                                let nongoal_sum: Floaty =
                                    nongoal_offsets // cost
                                        .iter()
                                        .enumerate()
                                        .map(|(shape_ix, offsets)| -> Floaty {
                                            FloatOrd(
                                                offsets
                                                    .iter()
                                                    .filter(|(x, y)| {
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
                                    .filter(|(this_goal_ix, (x, y))| {
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
                |goal_offset| *goal_offset == goal_target_offset_hack,
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

    #[derive(Clone, PartialEq, Eq)]
    struct OffsetAndPulverizedBlocks {
        offset: Offset,
        pulverized_blocks: BTreeSet<usize>, // i ∈ pulverized_blocks <=> we hit block (nongoal_offset_vec:goal_offset_vec)[i]
    }
    impl Ord for OffsetAndPulverizedBlocks {
        // a < b iff:
        // - a.pulverized_blocks ⊊ b.pulverized_blocks, or otherwise:
        // - a.pulverized_blocks < b..pulverized_blocks in the set-universe of coordinates, or otherwise:
        // - a.offset ≤ b.offset.
        fn cmp(&self, other: &Self) -> Ordering {
            std::cmp::Ordering::reverse({
                // We want a min-heap, after all
                let selfpb = &self.pulverized_blocks;
                let otherpb = &other.pulverized_blocks;
                if selfpb.is_subset(&otherpb) && selfpb != otherpb {
                    std::cmp::Ordering::Less
                } else if otherpb.is_subset(&selfpb) && selfpb != otherpb {
                    std::cmp::Ordering::Greater
                } else {
                    let selfiter = selfpb.iter();
                    let otheriter = otherpb.iter();
                    // comparison-sweep:
                    selfiter
                        .zip(otheriter)
                        .map(|(i, j)| i.cmp(j))
                        .find(|&cmp| cmp != std::cmp::Ordering::Equal)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(self.offset.cmp(&other.offset))
                }
            })
        }
    }
    impl PartialOrd for OffsetAndPulverizedBlocks {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    let pulverized_bfs = |goalvec_ix: usize| -> usize {
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
    // TODO: I have no idea if `&dyn Fn(Offset) -> bool` is the right signature as I didn't learn about `&dyn` yet
    let bfs_general = |initial_offset: Offset, is_legal: &dyn Fn(Offset) -> bool| {
        let mut legal_offsets: Vec<Offset> = Vec::new();
        let mut seen_offsets: BTreeSet<Offset> = BTreeSet::new();
        seen_offsets.insert(initial_offset);
        let mut stack: Vec<Offset> = vec![initial_offset];
        while let Some(offset) = stack.pop() {
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
    let bfs_nongoal = |movingshape_ix: usize,
                       moving_offset: Offset,
                       trimmed_movingshape_offsets: &Offsets|
     -> Vec<Offset> {
        let is_legal = |offsetty: Offset| -> bool {
            for (shape_ix, shape_offsets) in blockstate.nongoal_offsets.iter().enumerate() {
                if shape_ix == movingshape_ix {
                    continue;
                }
                for offset in shape_offsets {
                    // TODO: How bad are these "as usize" conversions?
                    // If they're really bad, I might just end up using usize as the type for Coor
                    if !nonintersectionkey[shape_ix][offset.0 as usize][offset.1 as usize]
                        [movingshape_ix][offsetty.0 as usize][offsetty.1 as usize]
                    {
                        return false;
                    }
                }
            }
            for (goalvec_ix, offset) in blockstate.goal_offsets.iter().enumerate() {
                let shapekey_ix = goal_shapekey_key[goalvec_ix];
                if !nonintersectionkey[shapekey_ix][offset.0 as usize][offset.1 as usize]
                    [movingshape_ix][offsetty.0 as usize][offsetty.1 as usize]
                {
                    return false;
                }
            }
            for offset in trimmed_movingshape_offsets {
                if !nonintersectionkey[movingshape_ix][offset.0 as usize][offset.1 as usize]
                    [movingshape_ix][offsetty.0 as usize][offsetty.1 as usize]
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
        let is_legal = |offsetty: Offset| -> bool {
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
        for mutated_offset in bfs_goal(goalvec_ix, *offset) {
            let mut new_blockstate = blockstate.clone();
            new_blockstate.goal_offsets[goalvec_ix] = mutated_offset;
            neighboring_blockstates.push(new_blockstate);
        }
    }
    for (shapekey_ix, offsets) in blockstate.nongoal_offsets.iter().enumerate() {
        for offset in offsets.iter() {
            let mut trimmed_shape_offsets = offsets.clone();
            trimmed_shape_offsets.remove(offset);
            for mutated_offset in bfs_nongoal(shapekey_ix, *offset, &trimmed_shape_offsets) {
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
    width: Coor,
    height: Coor,
) {
    // TODO: Colors aren't the best choice here, because colors WILL change after
    // blocks move, giving the illusion of some blocks having changed shapes

    // Create vec of blocks:
    let mut blocks: Vec<CoordinatesSet> = Vec::new();
    for (shape, offsets) in shapekey.iter().zip(blockstate.nongoal_offsets.iter()) {
        for offset in offsets {
            let block: CoordinatesSet = shape
                .iter()
                .map(|coor| (coor.0 + offset.0, coor.1 + offset.1))
                .collect();
            blocks.push(block);
        }
    }
    for (shapekey_ix, offset) in goal_shapekey_key.iter().zip(blockstate.goal_offsets.iter()) {
        let block: CoordinatesSet = shapekey[*shapekey_ix]
            .iter()
            .map(|coor| (coor.0 + offset.0, coor.1 + offset.1))
            .collect();
        blocks.push(block);
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

fn solve_puzzle_preprocessing(
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

pub fn solve_puzzle_own_bfs(start: &str, goal: &str) {
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
    ) = solve_puzzle_preprocessing(start, goal);

    print_puzzle(
        &bounds,
        &shapekey,
        &start_blockstate,
        &goal_shapekey_key,
        width,
        height,
    );
    println!();

    let path = puzzle_bfs_with_path_reconstruction(
        &start_blockstate,
        |blockstate| {
            get_neighboring_blockstates(
                blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
                width,
                height,
            )
        },
        |blockstate| blockstate.goal_offsets == goal_target_offsets,
    )
    .unwrap();
    assert!(path.len() < 1000); // TODO: Remove this
}

pub fn solve_puzzle_lib_bfs(start: &str, goal: &str) {
    let (
        _bounds,
        _shapekey,
        start_blockstate,
        nonintersectionkey,
        goal_shapekey_key,
        goal_target_offsets,
        _minkowski_diams,
        width,
        height,
    ) = solve_puzzle_preprocessing(start, goal);

    let path = pathfinding::directed::bfs::bfs(
        &start_blockstate,
        |blockstate| {
            get_neighboring_blockstates(
                blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
                width,
                height,
            )
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
    ) = solve_puzzle_preprocessing(start, goal);

    let (path, _) = pathfinding::directed::astar::astar(
        &start_blockstate,
        |blockstate| -> Vec<(Blockstate, Floaty)> {
            get_neighboring_blockstates(
                blockstate,
                &nonintersectionkey,
                &goal_shapekey_key,
                width,
                height,
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

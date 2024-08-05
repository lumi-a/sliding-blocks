use core::convert::From;
use wasm_bindgen::prelude::*;

// TODO: Perhaps conflate JsPuzzle and Puzzle by using Cow

#[wasm_bindgen(getter_with_clone)]
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct JsPuzzle {
    pub name: String,
    pub start: String,
    pub goal: String,
    pub min_moves: usize,
}

#[wasm_bindgen]
#[must_use]
#[inline]
pub fn js_get_all() -> Vec<JsPuzzle> {
    MULTI_GOAL_EXAMPLES
        .iter()
        .map(|puzzle| puzzle.clone().into())
        .chain(
            SINGLE_GOAL_EXAMPLES
                .iter()
                .map(|puzzle| puzzle.clone().into()),
        )
        .collect()
}

#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct Puzzle {
    pub name: &'static str,
    pub start: &'static str,
    pub goal: &'static str,
    pub min_moves: usize,
}
impl From<Puzzle> for JsPuzzle {
    #[inline]
    fn from(puzzle: Puzzle) -> Self {
        Self {
            name: puzzle.name.to_owned(),
            start: puzzle.start.to_owned(),
            goal: puzzle.goal.to_owned(),
            min_moves: puzzle.min_moves,
        }
    }
}

pub const MULTI_GOAL_EXAMPLES: [Puzzle; 8] = [
    Puzzle {
        name: "🎩 Ten-Step Solution",
        start: "  
      A  
    11.  
    .22. 
    .33. 
     .44 
     B   
    ",
        goal: "  
      B  
    ...  
    .... 
    .... 
     ... 
     A 
    ",
        min_moves: 10,
    },
    Puzzle {
        name: "🎩 The Elevator Switch",
        start: "
       AA   
       AA   
    ...a....
    ##aabb..
    ## cb XX
    ..ccddXX
    ....d...
       BB
       BB
    ",
        goal: "
       BB   
       BB   
    ........
    ........
    .. .. ..
    ........
    ........
       AA
       AA
    ",
        min_moves: 12,
    },
    Puzzle {
        name: "🎩 The Dragon Bridge",
        start: "
      ;;
     .;;t
    z...tt
    zz....
    .zv.rr
    .vv.r.
      v.
    ",
        goal: "
      rr
     .rv.
    ..vv..
    ..zv..
    ..zz..
    ..tz..
      tt
    ",
        min_moves: 19,
    },
    Puzzle {
        name: "🎩 Four Balls",
        start: "
        A
    B|.--
     |##.
     .##I
     ==.ID
     C
    ",
        goal: "
        C
    D....
     ....
     ....
     ....B
     A
    ",
        min_moves: 28,
    },
    Puzzle {
        name: "🎩 The Time Machine",
        start: "
      .yy
    ....y_
    .z.___.
    zz.#.rr
    zb.v.rq
     bbvvqq
      .v.
    ",
        goal: "
      ...  
    .rrzyy
    .rzzvy.
    .bz#vv.
    .bb_vq.
     ___qq.
      ...
    ",
        min_moves: 39,
    },
    Puzzle {
        name: "🎩 Jewel Swap",
        start: "
     .0011RR
     .0..1RR
    BB2..3.
    BB2233.
    ",
        goal: "
     .0011BB
     .0..1BB
    RR2..3.
    RR2233.
    ",
        min_moves: 61,
    },
    Puzzle {
        name: "🎩 The Greatest Escape 1",
        start: "
    L.MM|
    LlMM|
    KKPPI
    KKPPI
    --;==
    ..___
    ",
        goal: "
    .....
    .....
    .....
    .....
    ....L
    l...L
    ",
        min_moves: 71,
    },
    Puzzle {
        name: "🎩 The Greatest Escape 2",
        start: "
    l--+;L
    MM|.tL
    MM|ttt
    .==__.
     rrqq 
     $  ,
    ",
        goal: "
    ......
    ......
    ......
    ......
     L... 
     L  l
    ",
        min_moves: 107,
    },
];

pub const SINGLE_GOAL_EXAMPLES: [Puzzle; 13] = [
    Puzzle {
        name: "🎩 Cluttered Bag",
        start: "
    ....---
    rr55.n-
    r55.nn.
    cc..vv.
    ctttuvu
    cct#uuu
    ",
        goal: "
    ...#...
    .......
    .......
    .......
    .......
    .......
    ",
        min_moves: 20,
    },
    Puzzle {
        name: "🎩 Grab the Key",
        start: "
        |
        |
      ##.L.
      ##.LL
      l. QQ
      ll!QQ
      ..!..
        .
        .
      ",
        goal: "
        .
        .
      .....
      .....
      .. ..
      .....
      .....
        |
        |
      ",
        min_moves: 27,
    },
    Puzzle {
        name: "🎩 Taking Out the Trash",
        start: "
    ##
    ##
    PPKK..
    PPKK
    ==--
    __$$
      ..
      ..
    ",
        goal: "
    ..
    ..
    ......
    ....
    ....
    ....
      ##
      ##
    ",
        min_moves: 29,
    },
    Puzzle {
        name: "🎩 Garbage Disposal",
        start: "
      tt
      tt
    ......
    .ppoo.
     ypog
     yygg
      bb
      ..
    ",
        goal: "
      ..
      ..
    ......
    ......
     ....
     ....
      tt
      tt
    ",
        min_moves: 31,
    },
    Puzzle {
        name: "🎩 UFO SOS",
        start: "
     ccPPrr
     cbPPyr
     bbppyy
    ...pp...
     ##  ..
     ##  ..
    ",
        goal: "
     ......
     ......
     ......
    ........
     ..  ##
     ..  ##
    ",
        min_moves: 31,
    },
    Puzzle {
        name: "🎩 Impassable Gate",
        start: "
      ..  
      ..
    ---.rr
     .==r
    .qqyy.
    ggqypp
    g.##.p
    ..##..
    ",
        goal: "
      ##
      ##
    ......
     ....
    ......
    ......
    ......
    ......
    ",
        min_moves: 32,
    },
    Puzzle {
        name: "🎩 The Diabolical Box Reopened",
        start: "
     nnn..
    .n|n.cc
    ..|---c
    rr|#Icc
    r===Il.
    __.iIl.
     ..ill
    ",
        goal: "
     ..#..
    .......
    .......
    .......
    .......
    .......
     .....
    ",
        min_moves: 38,
    },
    Puzzle {
        name: "🎩 A Troublesome Table",
        start: "
    ....g.g....
    h.h.gggp.p.
    hhhbgkrppp.
    hbbbkkrrrp.
    .b.bkcrwr..
    ###ccc.www..
    ###c.c.w.w..
    ",
        goal: "
    ...........
    ...........
    ...........
    ...........
    ...........
    .........###
    .........###
    ",
        min_moves: 56,
    },
    Puzzle {
        name: "🎩 Sliding Pipes",
        start: "
    ####vv00
    |mm tv11
    |mmttt22
    .... ...
    33445566
    ",
        goal: "
    ........
    ... ....
    ........
    .... ...
    ....####
    ",
        min_moves: 68,
    },
    Puzzle {
        name: "🎩 The Diabolical Box",
        start: "
      ##    
      ##    
     +..L.  
    +++.LL  
    .+. --. 
    .^^|rr  
     .^|r.  
      ..    
      ..
    ",
        goal: "
      ..
      ..   
     .....  
    ......  
    ... ... 
    ......   
     .....  
      ##    
      ##
    ",
        min_moves: 78,
    },
    Puzzle {
        name: "🎩 Royal Escape",
        start: "
    --==+
    ##|/.
    ##|*.
    __mm$
    ",
        goal: "
    .....
    ...##
    ...##
    .....
    ",
        min_moves: 81,
    },
    Puzzle {
        name: "🎩 From Right to Left",
        start: "
      ippPPI
      ippPPI
    ...--PPI##
    ...==LKK##
      lkkLKK
      lkkLKK
    ",
        goal: "
      ......
      ......
    ##........
    ##........
      ......  
      ......
    ",
        min_moves: 101,
    },
    Puzzle {
        name: "🎩 Trapped Treasure",
        start: "
    0##1
    0##1
    2--3
    2==3
    4567
    4..8
    ",
        goal: "
    ....
    ....
    ....
    ....
    .##.
    .##.
    ",
        min_moves: 128,
    },
];

mod tests {

    #[allow(unused_imports)]
    use super::*;
    #[allow(unused_imports)]
    use itertools::Itertools;

    #[test]
    fn correct_multi_single_division() {
        let count_chars = |s: &str| {
            s.lines()
                .flat_map(|l| l.chars().filter(|c| !c.is_whitespace()))
                .counts()
                .len()
        };
        for p in MULTI_GOAL_EXAMPLES {
            assert!(count_chars(p.goal) > 2, "{}", p.name);
        }
        for p in SINGLE_GOAL_EXAMPLES {
            assert!(
                count_chars(p.goal) == 2,
                "{} {}",
                p.name,
                count_chars(p.goal)
            );
        }
    }

    #[test]
    fn sorted_examples() {
        fn is_sorted(data: impl IntoIterator<Item = Puzzle>) -> bool {
            data.into_iter()
                .tuple_windows()
                .all(|(a, b)| a.min_moves <= b.min_moves)
        }
        assert!(is_sorted(MULTI_GOAL_EXAMPLES));
        assert!(is_sorted(SINGLE_GOAL_EXAMPLES));
    }
}

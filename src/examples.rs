use wasm_bindgen::prelude::*;

// TODO: Perhaps conflate JsPuzzle and Puzzle by using Cow

#[wasm_bindgen(getter_with_clone)]
#[derive(Debug, Clone)]
pub struct JsPuzzle {
    pub name: String,
    pub start: String,
    pub goal: String,
    pub min_moves: usize,
}

#[wasm_bindgen]
pub fn get_all_js_examples() -> Vec<JsPuzzle> {
    ALL_EXAMPLES
        .iter()
        .map(|puzzle| puzzle.clone().into())
        .collect()
}

#[derive(Debug, Clone)]
pub struct Puzzle {
    pub name: &'static str,
    pub start: &'static str,
    pub goal: &'static str,
    pub min_moves: usize,
}
impl std::convert::From<Puzzle> for JsPuzzle {
    fn from(puzzle: Puzzle) -> Self {
        Self {
            name: puzzle.name.to_string(),
            start: puzzle.start.to_string(),
            goal: puzzle.goal.to_string(),
            min_moves: puzzle.min_moves,
        }
    }
}

pub const ALL_EXAMPLES: [Puzzle; 6] = [
    Puzzle {
        name: "Ten-Step Solution",
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
        name: "Four Balls",
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
        name: "The Elevator Switch",
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
        name: "Garbage Disposal",
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
        name: "The Diabolical Box",
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
        name: "Royal Escape",
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
];

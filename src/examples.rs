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
impl std::convert::Into<JsPuzzle> for Puzzle {
    fn into(self) -> JsPuzzle {
        JsPuzzle {
            name: self.name.to_string(),
            start: self.start.to_string(),
            goal: self.goal.to_string(),
            min_moves: self.min_moves,
        }
    }
}

pub const ALL_EXAMPLES: [Puzzle; 6] = [
    TEN_STEP_SOLUTION,
    FOUR_BALLS,
    THE_ELEVATOR_SWITCH,
    GARBAGE_DISPOSAL,
    THE_DIABOLICAL_BOX,
    ROYAL_ESCAPE,
];

pub const TEN_STEP_SOLUTION: Puzzle = Puzzle {
    name: "Ten-Step Solution",
    start: r#"  
      A  
    11.  
    .22. 
    .33. 
     .44 
     B   
    "#,
    goal: r#"  
      B  
    ...  
    .... 
    .... 
     ... 
     A 
    "#,
    min_moves: 10,
};

pub const FOUR_BALLS: Puzzle = Puzzle {
    name: "Four Balls",
    start: r#"
        A
    B|.--
     |##.
     .##I
     ==.ID
     C
    "#,
    goal: r#"
        C
    D....
     ....
     ....
     ....B
     A
    "#,
    min_moves: 28,
};

pub const THE_ELEVATOR_SWITCH: Puzzle = Puzzle {
    name: "The Elevator Switch",
    start: r#"
       AA   
       AA   
    ...a....
    ##aabb..
    ## cb XX
    ..ccddXX
    ....d...
       BB
       BB
    "#,
    goal: r#"
       BB   
       BB   
    ........
    ........
    .. .. ..
    ........
    ........
       AA
       AA
    "#,
    min_moves: 12,
};

pub const GARBAGE_DISPOSAL: Puzzle = Puzzle {
    name: "Garbage Disposal",
    start: r#"
      tt
      tt
    ......
    .ppoo.
     ypog
     yygg
      bb
      ..
    "#,
    goal: r#"
      ..
      ..
    ......
    ......
     ....
     ....
      tt
      tt
    "#,
    min_moves: 31,
};

pub const THE_DIABOLICAL_BOX: Puzzle = Puzzle {
    name: "The Diabolical Box",
    start: r#"
      ##    
      ##    
     +..L.  
    +++.LL  
    .+. --. 
    .^^|rr  
     .^|r.  
      ..    
      ..
    "#,
    goal: r#"
      ..
      ..   
     .....  
    ......  
    ... ... 
    ......   
     .....  
      ##    
      ##
    "#,
    min_moves: 78,
};

pub const ROYAL_ESCAPE: Puzzle = Puzzle {
    name: "Royal Escape",
    start: r#"
    --==+
    ##|/.
    ##|*.
    __mm$
    "#,
    goal: r#"
    .....
    ...##
    ...##
    .....
    "#,
    min_moves: 81,
};

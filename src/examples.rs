pub struct Puzzle {
    pub name: &'static str,
    pub start: &'static str,
    pub goal: &'static str,
    pub moves: usize,
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
    moves: 10,
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
    moves: 28,
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
    moves: 12,
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
    moves: 31,
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
    moves: 78,
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
    moves: 81,
};

pub struct Puzzle {
    pub start: &'static str,
    pub goal: &'static str,
}

pub const GARBAGE_DISPOSAL: Puzzle = Puzzle {
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
};

pub const DIABOLICAL_BOX: Puzzle = Puzzle {
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
};

pub const ROYAL_ESCAPE: Puzzle = Puzzle {
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
};

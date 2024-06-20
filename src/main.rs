pub fn solve_puzzle(start: &str, end: &str) {
    fn clean_string(s: &str) -> Vec<&str> {
        // Removes initial lines that are all whitespace
        // Trim maximum amount of leading whitespace in each line, such that
        // each line has the same amount of leading whitespace trimmed.
        // Same for trailing whitespace.
        let mut leading = 255; // ???????
        let mut trailing = 255; // Smells like recipe for disaster

        // remove final whitespace-only lines:
        let lines: Vec<&str> = s
            .lines()
            .rev()
            .skip_while(|l| l.trim().is_empty())
            .collect();
        // remove initial whitespace-only lines:
        let mut lines: Vec<&str> = lines
            .into_iter()
            .rev()
            .skip_while(|l| l.chars().all(|c| c.is_whitespace()))
            .collect();

        for line in &lines {
            let count_leading_whitespace = line.chars().take_while(|c| c.is_whitespace()).count();
            let count_trailing_whitespace =
                line.chars().rev().take_while(|c| c.is_whitespace()).count();
            leading = leading.min(count_leading_whitespace);
            trailing = trailing.min(count_trailing_whitespace);
        }
        println!("Leading: {}, trailing: {}", leading, trailing);

        // Remove leading and trailing whitespace from each line:
        for line in &mut lines {
            *line = &line[leading..(line.len() - trailing)];
        }
        lines
    }
    println!(
        "{}\n\n{}",
        clean_string(start).join("\r\n"),
        clean_string(end).join("\n")
    );
}

fn main() {
    let puzzle = (
        "


    #t#CCvv
    A.#DFFv

     #.#DDrr
    #.BEEuu
    #s#E
    ",
        "
    #.#D
    B.#DD
    #.#CC
    #.AEE
    #.#EFF
    ",
    );
    solve_puzzle(puzzle.0, puzzle.1);
}

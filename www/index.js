import { solve_puzzle } from "sliding-blocks"

let solution = solve_puzzle(`
    a.b
     .
`, `
    b.a
     .
`)
for (let step of solution) {
    console.log(step)
}
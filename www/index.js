import { solve_puzzle } from "sliding-blocks"

class Block {
    constructor(shape, offset) {
        this.shape = shape
    }
}

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
import { solve_puzzle } from "sliding-blocks"

// All I want for Christmas is proper JS data structures
type Point = number
const COMPONENT_BITS = 16
const COMPONENT_MASK = (1 << COMPONENT_BITS) - 1
const COMPONENT_SIGN_BIT = 1 << (COMPONENT_BITS - 1)
function p(x: number, y: number): Point {
    // x,y ∈ ℤ * 1/8, please
    const packedX = (Math.round(x * 8) & COMPONENT_MASK) >>> 0
    const packedY = (Math.round(y * 8) & COMPONENT_MASK) >>> 0
    return (packedX << COMPONENT_BITS) | packedY
}
function x(v: Point) {
    const x = (v >>> COMPONENT_BITS) & COMPONENT_MASK
    return ((x & COMPONENT_SIGN_BIT) ? (x | ~COMPONENT_MASK) : x) / 8
}
function y(v: Point) {
    const y = v & COMPONENT_MASK
    return ((y & COMPONENT_SIGN_BIT) ? (y | ~COMPONENT_MASK) : y) / 8
}
const shiftTuple = (v: Point, a: number, b: number) => p(a + x(v), b + y(v))
const shift = (v: Point, w: Point) => p(x(w) + x(v), y(w) + y(v))
const scale = (v: Point, s: number) => p(x(v) * s, y(v) * s)

type Shape = Set<Point>
type Offset = Point

function shape_to_path(shape: Shape) {
    type Edgepoint = Point
    const DIRS = [p(1, 0), p(0, -1), p(-1, 0), p(0, 1)]
    let edgepoint_to_dir: Map<Edgepoint, number> = new Map()
    shape.forEach(point => {
        const left = shiftTuple(point, -1, 0)
        const right = shiftTuple(point, 1, 0)
        const up = shiftTuple(point, 0, -1)
        const down = shiftTuple(point, 0, 1)
        // We maintain the invariant that we touch the shape with
        // our right hand. So if we're moving upwards, the shape
        // is not on the left:
        if (!shape.has(left)) {
            edgepoint_to_dir.set(shiftTuple(point, -0.5, 0), 1)
        }
        if (!shape.has(right)) {
            edgepoint_to_dir.set(shiftTuple(point, 0.5, 0), 3)
        }
        if (!shape.has(up)) {
            edgepoint_to_dir.set(shiftTuple(point, 0, -0.5), 0)
        }
        if (!shape.has(down)) {
            edgepoint_to_dir.set(shiftTuple(point, 0, 0.5), 2)
        }
    })
    let path = ""
    let counter_outer = 0
    let start_edgepoint: Edgepoint | undefined
    while ((start_edgepoint = edgepoint_to_dir.keys().next()?.value) && counter_outer++ < 20) {
        let edgepoint: Edgepoint = start_edgepoint
        let dir: number = edgepoint_to_dir.get(edgepoint)

        path += `M${x(edgepoint)} ${y(edgepoint)}`
        let counter_inner = 0
        do {
            const dirvec: Point = DIRS[dir]
            const forward = scale(dirvec, 0.5)
            const left = scale(DIRS[(dir + 1) % 4], 0.5)
            const right = scale(DIRS[(dir + 3) % 4], 0.5)
            const center = shift(edgepoint, forward)
            if (shape.has(shift(center, shift(left, forward)))) {
                edgepoint = shift(center, left)
                dir = (dir + 1) % 4
            } else if (shape.has(shift(center, shift(right, forward)))) {
                edgepoint = shift(center, forward)
            } else {
                edgepoint = shift(center, right)
                dir = (dir + 3) % 4
            }

            //path += `Q${x(center)} ${y(center)} ${x(edgepoint)} ${y(edgepoint)}`
            path += `C${x(center)} ${y(center)} ${x(center)} ${y(center)} ${x(edgepoint)} ${y(edgepoint)}`
            edgepoint_to_dir.delete(edgepoint)
        } while (edgepoint != start_edgepoint && counter_inner++ < 20)
        path += "Z"
    }
    return path
}

class Block {
    constructor(public shape: Shape, public offset: Offset) {
        this.shape = shape
        this.offset = offset
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
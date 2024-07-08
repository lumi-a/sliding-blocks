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

const PATH_CORNER_RADIUS = 0.1
function shape_to_path(shape: Shape) {
    enum Dir {
        Up,
        Down,
        Left,
        Right
    }
    interface EdgePoint {
        point: Point,
        dir: Dir // Invariant: You're always touching `shape` with your right hand
    }
    let edgepoints = new Collections.Set<EdgePoint>()
    shape.forEach(point => {
        const left = add(point, { x: -1, y: 0 })
        const right = add(point, { x: +1, y: 0 })
        const up = add(point, { x: 0, y: -1 })
        const down = add(point, { x: 0, y: 1 })
        if (!shape.contains(left)) {
            edgepoints.add({ point: add(point, { x: -0.5, y: 0 }), dir: Dir.Up })
        }
        if (!shape.contains(right)) {
            edgepoints.add({ point: add(point, { x: +0.5, y: 0 }), dir: Dir.Down })
        }
        if (!shape.contains(up)) {
            edgepoints.add({ point: add(point, { x: 0, y: -0.5 }), dir: Dir.Right })
        }
        if (!shape.contains(down)) {
            edgepoints.add({ point: add(point, { x: 0, y: +0.5 }), dir: Dir.Left })
        }
        console.log(edgepoints)
        console.log(left)
        console.log(right)
        console.log(up)
        console.log(!shape.contains(down))
        console.log(shape.contains(down))
        console.log(down)
        console.log(shape)
    })
    let path = ""
    let start_edgepoint: EdgePoint | undefined
    let counter_outer = 0
    while (!edgepoints.isEmpty() && counter_outer++ < 100) {
        start_edgepoint = edgepoints.toArray()[0] // TODO: This is dumb
        let edgepoint = Object.assign({}, start_edgepoint)
        path += `M${edgepoint.point.x} ${edgepoint.point.y}`
        let counter_inner = 0
        do {
            let center: Point | undefined
            switch (edgepoint.dir) {
                case Dir.Up:
                    center = add(edgepoint.point, { x: 0, y: -0.5 })
                    if (shape.contains(add(center, { x: -0.5, y: -0.5 }))) {
                        edgepoint = { point: add(center, { x: -0.5, y: 0 }), dir: Dir.Left }
                    } else if (shape.contains(add(center, { x: 0.5, y: -0.5 }))) {
                        edgepoint = { point: add(center, { x: 0, y: -0.5 }), dir: Dir.Up }
                    } else {
                        edgepoint = { point: add(center, { x: 0.5, y: 0 }), dir: Dir.Right }
                    }
                    break
                case Dir.Down:
                    center = add(edgepoint.point, { x: 0, y: 0.5 })
                    if (shape.contains(add(center, { x: 0.5, y: 0.5 }))) {
                        edgepoint = { point: add(center, { x: 0.5, y: 0 }), dir: Dir.Right }
                    } else if (shape.contains(add(center, { x: -0.5, y: 0.5 }))) {
                        edgepoint = { point: add(center, { x: 0, y: 0.5 }), dir: Dir.Down }
                    } else {
                        edgepoint = { point: add(center, { x: -0.5, y: 0 }), dir: Dir.Left }
                    }
                    break
                case Dir.Left:
                    center = add(edgepoint.point, { x: -0.5, y: 0 })
                    if (shape.contains(add(center, { x: -0.5, y: 0.5 }))) {
                        edgepoint = { point: add(center, { x: 0, y: 0.5 }), dir: Dir.Down }
                    } else if (shape.contains(add(center, { x: -0.5, y: -0.5 }))) {
                        edgepoint = { point: add(center, { x: -0.5, y: 0 }), dir: Dir.Left }
                    } else {
                        edgepoint = { point: add(center, { x: 0, y: -0.5 }), dir: Dir.Up }
                    }
                    break
                case Dir.Right:
                    center = add(edgepoint.point, { x: 0.5, y: 0 })
                    if (shape.contains(add(center, { x: 0.5, y: -0.5 }))) {
                        edgepoint = { point: add(center, { x: 0, y: -0.5 }), dir: Dir.Up }
                    } else if (shape.contains(add(center, { x: 0.5, y: 0.5 }))) {
                        edgepoint = { point: add(center, { x: 0.5, y: 0 }), dir: Dir.Right }
                    } else {
                        edgepoint = { point: add(center, { x: 0, y: 0.5 }), dir: Dir.Down }
                    }
            }
            path += `Q${edgepoint.point.x} ${edgepoint.point.y} ${center.x} ${center.y}`
            edgepoints.remove(edgepoint)
        } while ((edgepoint.point.x != start_edgepoint.point.x || edgepoint.point.y != start_edgepoint.point.y) && counter_inner++ < 100)
        path += "Z"
    }
    return path
}

const shape_0 = new Collections.Set<Point>()
for (let p of [{ x: 1, y: 0 }, { x: 0, y: 1 }, { x: 1, y: 2 }, { x: 2, y: 1 }]) {
    shape_0.add(Object.assign({}, p))
    console.log(shape_0)
}
console.log(shape_to_path(shape_0))

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
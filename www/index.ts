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
const shift_tuple = (v: Point, a: number, b: number) => p(a + x(v), b + y(v))
const shift = (v: Point, w: Point) => p(x(w) + x(v), y(w) + y(v))
const scale = (v: Point, s: number) => p(x(v) * s, y(v) * s)
const shift_shape = (shape: Shape, v: Point) => new Set([...shape].map(p => shift(p, v)))
const unshift_shape = (shape: Shape, v: Point) => {
    const w = scale(v, -1)
    return new Set([...shape].map(p => shift(p, w)))
}

type Shape = Set<Point>
type Offset = Point

function isSubset(a: Shape, b: Shape): boolean {
    return [...a].every(x => b.has(x))
}
function isDisjoint(a: Shape, b: Shape): boolean {
    return [...a].every(x => !b.has(x))
}

function get_extremes(coordinates_set: Shape): [Point, Point] {
    // Extract min-x and min-y.
    // Assumes that coordinatesSet is nonempty.
    // TODO: Is that a misassumption?
    let min_x = Number.MAX_SAFE_INTEGER
    let min_y = Number.MAX_SAFE_INTEGER
    let max_x = Number.MIN_SAFE_INTEGER
    let max_y = Number.MIN_SAFE_INTEGER
    for (const point of coordinates_set) {
        max_x = Math.min(max_x, x(point))
        max_y = Math.min(max_y, y(point))
        min_x = Math.min(min_x, x(point))
        min_y = Math.min(min_y, y(point))
    }
    return [p(min_x, min_y), p(max_x, max_y)]
}

function shape_to_path(shape: Shape): string {
    type Edgepoint = Point
    const DIRS = [p(1, 0), p(0, -1), p(-1, 0), p(0, 1)]
    let edgepoint_to_dir: Map<Edgepoint, number> = new Map()
    shape.forEach(point => {
        const left = shift_tuple(point, -1, 0)
        const right = shift_tuple(point, 1, 0)
        const up = shift_tuple(point, 0, -1)
        const down = shift_tuple(point, 0, 1)
        // We maintain the invariant that we touch the shape with
        // our right hand. So if we're moving upwards, the shape
        // is not on the left:
        if (!shape.has(left)) {
            edgepoint_to_dir.set(shift_tuple(point, -0.5, 0), 1)
        }
        if (!shape.has(right)) {
            edgepoint_to_dir.set(shift_tuple(point, 0.5, 0), 3)
        }
        if (!shape.has(up)) {
            edgepoint_to_dir.set(shift_tuple(point, 0, -0.5), 0)
        }
        if (!shape.has(down)) {
            edgepoint_to_dir.set(shift_tuple(point, 0, 0.5), 2)
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

function char_to_color(char: string): string {
    let code = char.charCodeAt(0)
    // TODO: Better color palette
    return `rgb(${code * 54979 % 255},${code * 70769 % 255},${code * 10113 % 255})`
}

const svg_elem = document.getElementById("svg")

class Block {
    public shape: Shape
    public offset: Offset
    public char: string
    public path: SVGPathElement

    constructor(coordinates: Shape, char: string) {
        let [min, max] = get_extremes(coordinates)
        this.shape = unshift_shape(coordinates, min)
        this.offset = min
        this.char = char
    }

    initialise() {
        let path = document.createElementNS("http://www.w3.org/2000/svg", "path")
        path.setAttribute("d", shape_to_path(this.shape))
        path.setAttribute("fill", char_to_color(this.char))
        // TODO: Stroke, shadow, letter-pattern-fill
        this.path = path

        svg_elem.appendChild(this.path)

        this.updateTranslation()
    }
    updateTranslation() {
        this.path.setAttribute("transform", `translate(${x(this.offset)}, ${y(this.offset)})`)
    }
    getCoordinates(): Shape {
        return shift_shape(this.shape, this.offset)
    }
    /*
    addDragging(blockstate, puzzle) {
        // this.containerBlock.addEventListener('pointerdown', event => { startDrag(event, this, blockstate, puzzle) })
        let hammer = new Hammer(this.containerBlock)
        hammer.on('panstart', event => startDrag(event, this, blockstate, puzzle, hammer))
    }
    */
}

const BOUNDS_CHAR: string = '.'
class Blockstate {
    public min: Point
    public max: Point
    public bounds: Block
    public blocks: Array<Block>

    constructor(bounds: Block, blocks: Array<Block>) {
        [this.min, this.max] = get_extremes(bounds.shape)
        this.bounds = bounds
        this.blocks = blocks
    }

    isValid() {
        const bounds = this.bounds.getCoordinates()
        const blocks = this.blocks.map(b => b.getCoordinates())
        for (let block_ix = 0; block_ix < blocks.length; block_ix++) {
            const block = blocks[block_ix]
            if (!isSubset(block, bounds)) return false
            for (let j = block_ix + 1; j < blocks.length; j++) {
                const other_block = blocks[j]
                if (!isDisjoint(block, other_block)) return false
            }
        }
        return true
    }

    static blockstateFromString(str: string) {
        const lines = str.split(/\r?\n/g)

        let global_min_x = Number.MAX_SAFE_INTEGER;
        let global_min_y = Number.MAX_SAFE_INTEGER;

        let char_to_blockcoordinates: { [key: string]: Shape } = {}
        let bounds_coordinates: Shape = new Set()

        let Y = 0
        for (let line of lines) {
            let X = 0
            for (let c of [...line]) {
                if (!(/\s/.test(c))) { // isn't out-of-bounds
                    global_min_x = Math.min(global_min_x, X)
                    global_min_y = Math.min(global_min_y, Y)
                    if (c != BOUNDS_CHAR) {
                        char_to_blockcoordinates[c] = char_to_blockcoordinates[c] || new Set()
                        char_to_blockcoordinates[c].add(p(X, Y))
                    }
                    bounds_coordinates.add(p(X, Y))
                }
                X++
            }
            Y++
        }
        const min = p(global_min_x, global_min_y)

        let bounds: Block = new Block(unshift_shape(bounds_coordinates, min), BOUNDS_CHAR)
        let blocks: Array<Block> = Object.entries(char_to_blockcoordinates).map(([c, coordinates]) => new Block(coordinates, c))

        return new Blockstate(bounds, blocks)
    }

    initialise() {
        svg_elem.innerHTML = ""
        const width = x(this.max) - x(this.min) + 1
        const height = y(this.max) - y(this.min) + 1
        svg_elem.setAttribute("viewBox", `${x(this.min) - 1} ${y(this.min) - 1} ${width + 1} ${height + 1}`)

        this.bounds.initialise()
        for (let block of this.blocks) block.initialise()
    }
}

let bs = Blockstate.blockstateFromString(`
      tt
      tt
    ......
    .ppoo.
     ypog
     yygg
      bb
      ..
`)
bs.initialise()

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
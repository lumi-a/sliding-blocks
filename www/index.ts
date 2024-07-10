import { solve_puzzle } from "sliding-blocks"
import interact from 'interactjs' // TODO: Only install the parts you need: https://interactjs.io/docs/installation#npm-streamlined

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
const shift = (v: Point, w: Point) => p(x(v) + x(w), y(v) + y(w))
const unshift = (v: Point, w: Point) => shift(v, scale(w, -1))
const scale = (v: Point, s: number) => p(x(v) * s, y(v) * s)
const shift_shape = (shape: Shape, v: Point) => new Set([...shape].map(p => shift(p, v)))
const unshift_shape = (shape: Shape, v: Point) => shift_shape(shape, scale(v, -1))

type Shape = Set<Point>
type Offset = Point

function is_subset(a: Shape, b: Shape): boolean {
    return [...a].every(x => b.has(x))
}
function is_disjoint(a: Shape, b: Shape): boolean {
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
        max_x = Math.max(max_x, x(point))
        max_y = Math.max(max_y, y(point))
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
    let outer_failsafe_counter = 0
    let start_edgepoint: Edgepoint | undefined
    while ((start_edgepoint = edgepoint_to_dir.keys().next()?.value) && outer_failsafe_counter++ < 1000) {
        let edgepoint: Edgepoint = start_edgepoint
        let dir: number = edgepoint_to_dir.get(edgepoint)

        path += `M${x(edgepoint)} ${y(edgepoint)}`
        let inner_failsafe_counter = 0
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
        } while (edgepoint != start_edgepoint && inner_failsafe_counter++ < 1000)
        path += "Z"
    }
    return path
}

const BOUNDS_CHAR: string = '.'

function char_to_color(char: string): string {
    let code = char.charCodeAt(0)
    // TODO: Better color palette
    // TODO: Dark mode stuff?
    if (char == BOUNDS_CHAR) {
        return `#BBB`
    }
    return `rgb(${code * 54979 % 255},${code * 70769 % 255},${code * 10113 % 255})`
}

const SVG_NAMESPACE = "http://www.w3.org/2000/svg"
const svg_puzzle_container = document.getElementById("svg-puzzle-container")
const svg_puzzle: SVGSVGElement = document.createElementNS(SVG_NAMESPACE, "svg")
svg_puzzle.id = "svg-puzzle"
svg_puzzle.setAttribute("xmlns", SVG_NAMESPACE)
svg_puzzle_container.appendChild(svg_puzzle)


class Block {
    public shape: Shape
    public offset: Offset
    public char: string
    public path: SVGPathElement
    public svg_elem: SVGSVGElement
    private svg_sctm_inverse: SVGMatrix
    private svg_point: SVGPoint


    constructor(coordinates: Shape, char: string) {
        let [min, max] = get_extremes(coordinates)
        this.shape = unshift_shape(coordinates, min)
        this.offset = min
        this.char = char
    }

    client_to_point(clientX: number, clientY: number): Point {
        this.svg_point.x = clientX
        this.svg_point.y = clientY
        const point = this.svg_point.matrixTransform(this.svg_sctm_inverse)
        return p(point.x, point.y)
    }

    update_translation() {
        this.path.style.transform = `translate(${x(this.offset)}px, ${y(this.offset)}px)`;
    }
    get_coordinates(): Shape {
        return shift_shape(this.shape, this.offset)
    }
    initialise_elem(svg_elem: SVGSVGElement) {
        let path = document.createElementNS("http://www.w3.org/2000/svg", "path")
        path.setAttribute("d", shape_to_path(this.shape))
        path.setAttribute("fill", char_to_color(this.char))
        // TODO: Stroke, shadow, letter-pattern-fill
        this.path = path

        this.svg_elem = svg_elem
        this.svg_sctm_inverse = svg_elem.getScreenCTM().inverse()
        this.svg_point = svg_elem.createSVGPoint()

        svg_elem.appendChild(this.path)

        this.update_translation()
    }

    make_interactive(puzzle: Puzzle, blockarr_ix: number) {
        if (!this.path || !this.svg_elem) {
            return
        }
        const blockstate = puzzle.blockstate
        const valid_offset = (offset: Point) => {
            const block = shift_shape(this.shape, offset)
            // TODO: Maybe don't re-calculate bounds.get_coordinates every time?
            // Also see is_valid function on blockstate class
            if (!is_subset(block, blockstate.bounds.get_coordinates())) return false
            for (let j = 0; j < blockstate.blocks.length; j++) {
                if (j == blockarr_ix) continue
                const other_block = blockstate.blocks[j].get_coordinates()
                if (!is_disjoint(block, other_block)) return false
            }
            return true
        }
        let start_cursor_position: Point
        let start_offset = this.offset
        const start = (event: any) => {
            start_offset = this.offset
            start_cursor_position = unshift(this.client_to_point(event.client.x, event.client.y), start_offset)
        }
        const move = (event: any) => {
            const old_offset = this.offset
            const new_offset_unrounded = unshift(this.client_to_point(event.client.x, event.client.y), start_cursor_position)

            const new_offset = p(Math.round(x(new_offset_unrounded)), Math.round(y(new_offset_unrounded)))
            if (old_offset != new_offset && valid_offset(new_offset)) {
                // Do BFS from old_offset to new_offset
                // TODO: Perhaps do proper animations rather than CSS transitions
                let queue = [old_offset]
                let visited = new Set([old_offset])
                let success = false
                while (queue.length > 0) {
                    let point: Point = queue.shift()
                    if (point == new_offset) {
                        success = true
                        break
                    }
                    for (let neighbor of [shift(point, p(0, 1)), shift(point, p(0, -1)), shift(point, p(1, 0)), shift(point, p(-1, 0))]) {
                        if (visited.has(neighbor) || !valid_offset(neighbor)) continue
                        visited.add(neighbor)
                        queue.push(neighbor)
                    }
                }
                if (success) {
                    this.offset = new_offset
                    this.update_translation()
                }

            }
        }
        const end = (event: any) => {
            // Is new blockstate different from before?
            if (start_offset !== this.offset) {
                puzzle.move_counter += 1
                console.log(`move ${puzzle.move_counter}`)
                if (puzzle.won()) {
                    // TODO: Confetti
                    console.log("won!")
                }
            }
        }

        // TODO: Highlighting the block that's currently being dragged

        interact(this.path).draggable({
            listeners: {
                start: start,
                move: move,
                end: end
            }
        })
    }
}

class Blockstate {
    public min: Point
    public max: Point
    public bounds: Block
    public blocks: Array<Block>
    public svg_elem: SVGSVGElement

    constructor(bounds: Block, blocks: Array<Block>) {
        [this.min, this.max] = get_extremes(bounds.shape)
        this.bounds = bounds
        this.blocks = blocks
    }

    is_valid() {
        const bounds = this.bounds.get_coordinates()
        const blocks = this.blocks.map(b => b.get_coordinates())
        for (let block_ix = 0; block_ix < blocks.length; block_ix++) {
            const block = blocks[block_ix]
            if (!is_subset(block, bounds)) return false
            for (let j = block_ix + 1; j < blocks.length; j++) {
                const other_block = blocks[j]
                if (!is_disjoint(block, other_block)) return false
            }
        }
        return true
    }

    static blockstate_from_string(str: string) {
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
        let blocks: Array<Block> = Object.entries(char_to_blockcoordinates).map(([c, coordinates]) => new Block(unshift_shape(coordinates, min), c))

        return new Blockstate(bounds, blocks)
    }

    initialise(svg_elem: SVGSVGElement) {
        svg_elem.innerHTML = ""
        const width = x(this.max) - x(this.min) + 1
        const height = y(this.max) - y(this.min) + 1
        svg_elem.setAttribute("viewBox", `${x(this.min) - 1} ${y(this.min) - 1} ${width + 1} ${height + 1}`)
        this.svg_elem = svg_elem

        this.bounds.initialise_elem(svg_elem)
        for (let block of this.blocks) block.initialise_elem(svg_elem)
    }

    make_interactive(puzzle: Puzzle) {
        for (let blockarr_ix = 0; blockarr_ix < this.blocks.length; blockarr_ix++) {
            const block = this.blocks[blockarr_ix]
            block.make_interactive(puzzle, blockarr_ix)
        }
    }
}

class Puzzle {
    public start_string: string
    public goal_string: string
    public min_moves: number | null
    public move_counter: number // TODO: Should this be a global variable instead?
    public blockstate: Blockstate
    public goal_offsets: Array<Offset | null> // TODO: This is awful

    constructor(start_string: string, goal_string: string, min_moves: number | null = null) {
        this.start_string = start_string
        this.goal_string = goal_string
        this.min_moves = min_moves
        this.move_counter = 0

        // TODO: Lots of error-checking on shapes that's also done in Rust
        this.blockstate = Blockstate.blockstate_from_string(start_string)
        const goal_blockstate = Blockstate.blockstate_from_string(goal_string)
        this.goal_offsets = this.blockstate.blocks.map(block => {
            for (let goal_block of goal_blockstate.blocks) {
                if (goal_block.char === block.char) return goal_block.offset
            }
            return null
        })
    }

    won() {
        const goal_offsets = this.goal_offsets
        console.log(goal_offsets, this.blockstate.blocks.map(b => b.offset))
        return this.blockstate.blocks.every((b, ix) => goal_offsets[ix] === null || b.offset === goal_offsets[ix])
    }

    initialise(svg_elem: SVGSVGElement) {
        this.blockstate = Blockstate.blockstate_from_string(this.start_string)
        this.blockstate.initialise(svg_elem)
        this.blockstate.make_interactive(this)
        this.move_counter = 0
    }
}

const change_puzzle_dialog = document.getElementById("change-puzzle-dialog") as HTMLDialogElement
const change_puzzle_btn = document.getElementById("change-puzzle-btn") as HTMLButtonElement
change_puzzle_btn.addEventListener("click", () => {
    change_puzzle_dialog.showModal()
})

let puzzle = new Puzzle(`
      tt
      tt
    ......
    .ppoo.
     ypog
     ....
      ..
      ..
`, `
      ..
      ..
    ......
    ......
     ....
     ....
      tt
      tt
`)
puzzle.initialise(svg_puzzle)
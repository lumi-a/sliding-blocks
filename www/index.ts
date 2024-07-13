import './style.css'
import { solve_puzzle, get_all_js_examples } from "sliding-blocks"
import interact from 'interactjs' // TODO: Only install the parts you need: https://interactjs.io/docs/installation#npm-streamlined
import JSConfetti from 'js-confetti'
const jsConfetti = new JSConfetti()

// TODO: Indicate goal

// All I want for Christmas is proper JS data structures
type Point = number
const COMPONENT_BITS = 16
const COMPONENT_MASK = (1 << COMPONENT_BITS) - 1
const COMPONENT_SIGN_BIT = 1 << (COMPONENT_BITS - 1)
const COMPONENT_PRECISION = 16
function p(x: number, y: number): Point {
    // x,y ∈ ℤ * 1/COMPONENT_PRECISION, please
    const packedX = (Math.round(x * COMPONENT_PRECISION) & COMPONENT_MASK) >>> 0
    const packedY = (Math.round(y * COMPONENT_PRECISION) & COMPONENT_MASK) >>> 0
    return (packedX << COMPONENT_BITS) | packedY
}
function x(v: Point) {
    const x = (v >>> COMPONENT_BITS) & COMPONENT_MASK
    return ((x & COMPONENT_SIGN_BIT) ? (x | ~COMPONENT_MASK) : x) / COMPONENT_PRECISION
}
function y(v: Point) {
    const y = v & COMPONENT_MASK
    return ((y & COMPONENT_SIGN_BIT) ? (y | ~COMPONENT_MASK) : y) / COMPONENT_PRECISION
}
const shift_tuple = (v: Point, a: number, b: number) => p(a + x(v), b + y(v))
const shift = (v: Point, w: Point) => p(x(v) + x(w), y(v) + y(w))
const unshift = (v: Point, w: Point) => shift(v, scale(w, -1))
const scale = (v: Point, s: number) => p(x(v) * s, y(v) * s)
const shift_shape = (shape: Shape, v: Point) => new Set([...shape].map(p => shift(p, v)))
const unshift_shape = (shape: Shape, v: Point) => shift_shape(shape, scale(v, -1))
const dot = (v: Point, w: Point) => x(v) * x(w) + y(v) * y(w)

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

function shape_to_path(shape: Shape, inset: number): string {
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
    let failsafe_counter = 0
    let start_edgepoint: Edgepoint | undefined
    while ((start_edgepoint = edgepoint_to_dir.keys().next()?.value) && failsafe_counter++ < 10000) {
        let edgepoint: Edgepoint = start_edgepoint
        let dir: number = edgepoint_to_dir.get(edgepoint)
        let drawpoint = shift(edgepoint, scale(DIRS[(dir + 3) % 4], inset))

        path += `M${x(drawpoint)} ${y(drawpoint)}`
        do {
            const dirvec: Point = DIRS[dir]
            const forward = scale(dirvec, 0.5)
            const backward = scale(dirvec, -0.5)
            const left = scale(DIRS[(dir + 1) % 4], 0.5)
            const right = scale(DIRS[(dir + 3) % 4], 0.5)
            const center = shift(edgepoint, forward)
            const horizontal = DIRS[(dir + 1) % 4]
            const vertical = dirvec
            if (shape.has(shift(center, shift(left, forward)))) {
                edgepoint = shift(center, left)
                dir = (dir + 1) % 4
            } else if (shape.has(shift(center, shift(right, forward)))) {
                edgepoint = shift(center, forward)
            } else {
                edgepoint = shift(center, right)
                dir = (dir + 3) % 4
            }

            const old_drawpoint = drawpoint
            drawpoint = shift(edgepoint, scale(DIRS[(dir + 3) % 4], inset))
            const center_drawpoint = shift(center, shift(scale(horizontal, dot(horizontal, unshift(old_drawpoint, center))), scale(vertical, dot(vertical, unshift(drawpoint, center)))))

            path += `C${x(center_drawpoint)} ${y(center_drawpoint)} ${x(center_drawpoint)} ${y(center_drawpoint)} ${x(drawpoint)} ${y(drawpoint)}`
            edgepoint_to_dir.delete(edgepoint)
        } while (edgepoint != start_edgepoint && failsafe_counter++ < 10000)
        path += "Z"
    }
    return path
}

const BOUNDS_CHAR: string = '.'

function char_to_color(char: string, lightness: number = 0.75, alpha: number = 1): string {
    // TODO: Better color palette
    // TODO: Dark mode stuff?
    if (char == BOUNDS_CHAR) {
        return `#BBB`
    }
    const code = char.charCodeAt(0)
    const chroma = 0.2
    const hue = (code * 65557) % 360
    return `oklch(${lightness} ${chroma} ${hue} / ${alpha})`
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
    construct_elem_path(inset: number = 1 / 32): SVGPathElement {
        let path = document.createElementNS(SVG_NAMESPACE, "path")
        path.setAttribute("d", shape_to_path(this.shape, this.char === BOUNDS_CHAR ? -inset : inset))
        path.setAttribute("fill-rule", "evenodd")
        path.setAttribute("stroke", char_to_color(this.char, 0.25))
        if (this.char !== BOUNDS_CHAR) {
            const pattern_id = `block-pattern-${this.char}`
            path.setAttribute("fill", `url(#${pattern_id})`)
        } else {
            path.setAttribute("fill", "#CCC")
        }

        return path
    }
    create_elem_pattern(svg_elem: SVGSVGElement) {
        if (this.char !== BOUNDS_CHAR) {
            const defs: SVGDefsElement = svg_elem.querySelector("defs") ?? svg_elem.appendChild(document.createElementNS(SVG_NAMESPACE, "defs"))
            const pattern = document.createElementNS(SVG_NAMESPACE, "pattern")
            const pattern_id = `block-pattern-${this.char}`

            const [_, max] = get_extremes(this.shape)
            pattern.setAttribute("id", pattern_id)
            pattern.setAttribute("patternUnits", "userSpaceOnUse")
            pattern.setAttribute("x", "-0.5")
            pattern.setAttribute("y", "-0.5")
            pattern.setAttribute("width", `${x(max) + 1}`)
            pattern.setAttribute("height", `${y(max) + 1}`)
            const rect = document.createElementNS(SVG_NAMESPACE, "rect")
            rect.setAttribute("width", `${x(max) + 1}`)
            rect.setAttribute("height", `${y(max) + 1}`)
            rect.setAttribute("fill", char_to_color(this.char))
            pattern.appendChild(rect)
            const code = this.char.charCodeAt(0)
            for (let p of this.shape) {
                for (let i = 0; i < 10; i++) {
                    const letter = document.createElementNS(SVG_NAMESPACE, "text") as SVGTextElement
                    const X = (52.092819851131311425 * code + 13.087032255978894794 * i) % 1 + x(p)
                    const Y = (28.640673508054986905 * code + 94.824207530838495049 * i) % 1 + y(p)
                    const d = (14.336965871263130613 * code + 40.125163576904165817 * i) % 360
                    const a = (72.313644289540589845 * code + 61.413884855320691933 * i) % 0.5
                    const f = (75.427959404814958242 * code + 85.346753489292519779 * i) % 0.4 + 0.2
                    letter.setAttribute("x", X.toString())
                    letter.setAttribute("y", Y.toString())
                    letter.setAttribute("fill", char_to_color(this.char, 0.5, a))
                    letter.setAttribute("font-size", f.toString())
                    letter.setAttribute("transform", `rotate(${d} ${X} ${Y})`)
                    letter.textContent = this.char
                    pattern.appendChild(letter)
                }
            }
            defs.appendChild(pattern)
        }
    }
    initialise_elem(svg_elem: SVGSVGElement) {
        this.create_elem_pattern(svg_elem)
        const path = this.construct_elem_path()

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
            this.path.classList.add("dragging")
        }
        const move = (event: any) => {
            const old_offset = this.offset
            const new_offset_unrounded = unshift(this.client_to_point(event.client.x, event.client.y), start_cursor_position)

            const new_offset = p(Math.round(x(new_offset_unrounded)), Math.round(y(new_offset_unrounded)))
            if (old_offset != new_offset && valid_offset(new_offset)) {
                // Do BFS from old_offset to new_offset
                // TODO: Perhaps do proper animations rather than CSS transitions
                // TODO: Rather than do BFS every move-call, do bfs once for valid offsets in start-call and check if this is contained in valid offsets
                // TODO: Better yet, do Floyd-Warshall and then "guess" the position the user wanted to move the block to by using the closest legal position (there always exists one because the initial position is legal)
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
                move_counter_elem.textContent = puzzle.move_counter.toString()
                if (puzzle.won()) {
                    if (puzzle.min_moves !== null && puzzle.move_counter < puzzle.min_moves) jsConfetti.addConfetti({ emojis: ["🐞"] })
                    else if (puzzle.min_moves !== null && puzzle.move_counter === puzzle.min_moves) jsConfetti.addConfetti({ emojis: ["🏆"] })
                    else jsConfetti.addConfetti()
                }
            }
            this.path.classList.remove("dragging")
        }

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

    to_string() {
        // TODO: How many of these methods should I import from Rust?
        let str = ""
        const block_coordinates: Array<[Block, Shape]> = this.blocks.map(b => [b, b.get_coordinates()])
        const bounds = this.bounds.get_coordinates()
        for (let Y = y(this.min); Y <= y(this.max); Y++) {
            for (let X = x(this.min); X <= x(this.max); X++) {
                let found = false
                for (let [b, coordinates] of block_coordinates) {
                    if (coordinates.has(p(X, Y))) {
                        str += b.char
                        found = true
                        break
                    }
                }
                if (!found) {
                    if (bounds.has(p(X, Y))) str += BOUNDS_CHAR
                    else str += " "
                }
            }
            str += "\n"
        }
        return str
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
        this.min_moves = min_moves
        this.move_counter = 0

        // TODO: Lots of error-checking on shapes that's also done in Rust
        const start_blockstate = Blockstate.blockstate_from_string(start_string)
        this.blockstate = start_blockstate
        const goal_blockstate = Blockstate.blockstate_from_string(goal_string)
        this.goal_offsets = this.blockstate.blocks.map(block => {
            for (let goal_block of goal_blockstate.blocks) {
                if (goal_block.char === block.char) return goal_block.offset
            }
            return null
        })
        this.start_string = start_blockstate.to_string()
        this.goal_string = goal_blockstate.to_string()
    }

    won() {
        const goal_offsets = this.goal_offsets
        return this.blockstate.blocks.every((b, ix) => goal_offsets[ix] === null || b.offset === goal_offsets[ix])
    }

    initialise(svg_elem: SVGSVGElement) {
        this.blockstate = Blockstate.blockstate_from_string(this.start_string)
        this.blockstate.initialise(svg_elem)
        this.blockstate.make_interactive(this)
        this.move_counter = 0
        move_counter_elem.textContent = "0"
    }
}

const change_puzzle_dialog = document.getElementById("change-puzzle-dialog") as HTMLDialogElement
const change_puzzle_btn = document.getElementById("change-puzzle-btn") as HTMLButtonElement
const puzzle_textarea_start = document.getElementById("puzzle-textarea-start") as HTMLTextAreaElement
const puzzle_textarea_goal = document.getElementById("puzzle-textarea-goal") as HTMLTextAreaElement
const puzzle_submit_btn = document.getElementById("puzzle-submit-btn") as HTMLButtonElement
const move_counter_elem = document.getElementById("move-counter") as HTMLSpanElement
change_puzzle_btn.addEventListener("click", () => {
    puzzle_textarea_goal.value = current_puzzle.goal_string
    puzzle_textarea_start.value = current_puzzle.start_string
    change_puzzle_dialog.showModal()
})

puzzle_submit_btn.addEventListener("click", e => {
    e.preventDefault()
    const start_string = puzzle_textarea_start.value
    const goal_string = puzzle_textarea_goal.value
    change_puzzle_dialog.close()
    current_puzzle = new Puzzle(start_string, goal_string)
    current_puzzle.initialise(svg_puzzle)
})

const puzzle_selection = document.getElementById("puzzle-selection") as HTMLSelectElement
const predefined_puzzles = get_all_js_examples()
for (let js_puzzle of predefined_puzzles) {
    const option = document.createElement("option")
    option.textContent = js_puzzle.name
    option.value = js_puzzle.name
    puzzle_selection.appendChild(option)
}
puzzle_selection.addEventListener("change", () => {
    // TODO: Use indices as names instead
    const name = puzzle_selection.value
    const js_puzzle = predefined_puzzles.find(js_puzzle => js_puzzle.name === name)!
    puzzle_textarea_start.value = js_puzzle.start
    puzzle_textarea_goal.value = js_puzzle.goal
})

let current_puzzle: Puzzle = new Puzzle(`
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
`, 5)
current_puzzle.initialise(svg_puzzle)
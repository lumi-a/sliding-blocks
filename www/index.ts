import './style.css'
import { solve_puzzle, get_all_js_examples } from "sliding-blocks"
import interact from 'interactjs' // TODO: Only install the parts you need: https://interactjs.io/docs/installation#npm-streamlined
import JSConfetti from 'js-confetti'
const jsConfetti = new JSConfetti()

// All I want for Christmas is proper JS data structures
type Point = number
const COMPONENT_BITS = 16
const COMPONENT_MASK = (1 << COMPONENT_BITS) - 1
const COMPONENT_SIGN_BIT = 1 << (COMPONENT_BITS - 1)
const COMPONENT_PRECISION = 32
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
const square = (x: number) => x * x
const dist2 = (v: Point, w: Point) => Math.sqrt(square(x(v) - x(w)) + square(y(v) - y(w)))
const dist1 = (v: Point, w: Point) => Math.abs(x(v) - x(w)) + Math.abs(y(v) - y(w))

type Shape = Set<Point>
type Offset = Point

function is_subset(a: Shape, b: Shape): boolean {
    return [...a].every(x => b.has(x))
}
function is_disjoint(a: Shape, b: Shape): boolean {
    return [...a].every(x => !b.has(x))
}
function union(shapes: Shape[]): Shape {
    let union: Shape = new Set()
    for (const shape of shapes) {
        for (const point of shape) {
            union.add(point)
        }
    }
    return union
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
    // TODO: Dark mode stuff?
    if (char == BOUNDS_CHAR) {
        return `#BBB`
    }
    const code = char.codePointAt(0)
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
    private walking_offsets: Array<Offset>
    private walking_offsets_time: number


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

    update_walking_offsets(walk: Array<Offset>) {
        this.walking_offsets_time += walk.length
        this.walking_offsets = walk.concat(this.walking_offsets)
    }

    update_translation(dt: number) {
        const time_scale = 100
        const walking_offsets = this.walking_offsets
        const walking_offsets_time = this.walking_offsets_time
        if (walking_offsets && walking_offsets.length > 0) {
            const new_offsets_time = Math.max(0, walking_offsets_time - dt * Math.ceil(walking_offsets_time) ** 0.75 / time_scale)
            this.walking_offsets_time = new_offsets_time
            while (walking_offsets.length > new_offsets_time + 2) {
                walking_offsets.pop()
            }
            if (walking_offsets.length === 1) {
                const offy = this.walking_offsets[0]
                this.path.style.transform = `translate(${x(offy)}px, ${y(offy)}px)`
            } else {
                const lambda = new_offsets_time % 1
                const from = this.walking_offsets[Math.ceil(new_offsets_time)]
                const to = this.walking_offsets[Math.floor(new_offsets_time)]
                this.path.style.transform = `translate(${lambda * x(from) + (1 - lambda) * x(to)}px, ${lambda * y(from) + (1 - lambda) * y(to)}px)`
            }
        }
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
            const code = this.char.codePointAt(0)
            for (let p of this.shape) {
                for (let i = 0; i < 10; i++) {
                    const letter = document.createElementNS(SVG_NAMESPACE, "text") as SVGTextElement
                    const X = (52.092819851131311425 * code + 13.087032255978894794 * i) % 1 + x(p)
                    const Y = (28.640673508054986905 * code + 94.824207530838495049 * i) % 1 + y(p)
                    const d = (14.336965871263130613 * code + 40.125163576904165817 * i) % 360
                    const a = (72.313644289540589845 * code + 61.413884855320691933 * i) % 0.25
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

        this.walking_offsets = []
        this.walking_offsets_time = 0
        this.update_walking_offsets([this.offset])
        this.update_translation(0)
    }

    make_interactive(puzzle: Puzzle, blockarr_ix: number) {
        if (!this.path || !this.svg_elem) {
            return
        }
        const blockstate = puzzle.blockstate
        const bounds = blockstate.bounds.get_coordinates()
        let start_cursor_position: Point
        let start_offset = this.offset
        let valid_offsets: Set<Point>
        const start = (event: any) => {
            start_offset = this.offset
            start_cursor_position = unshift(this.client_to_point(event.client.x, event.client.y), start_offset)
            this.path.classList.add("dragging")

            const other_blocks_union = union(blockstate.blocks.filter((b, i) => i != blockarr_ix).map(b => b.get_coordinates()))
            const offset_is_valid = (offset: Point) => {
                const block = shift_shape(this.shape, offset)
                if (!is_subset(block, bounds)) return false
                if (!is_disjoint(block, other_blocks_union)) return false
                return true
            }

            let queue = [start_offset]
            valid_offsets = new Set([start_offset])
            while (queue.length > 0) {
                let point: Point = queue.shift()
                for (let neighbor of [shift(point, p(0, 1)), shift(point, p(0, -1)), shift(point, p(1, 0)), shift(point, p(-1, 0))]) {
                    if (valid_offsets.has(neighbor) || !offset_is_valid(neighbor)) continue
                    valid_offsets.add(neighbor)
                    queue.push(neighbor)
                }
            }
        }
        const move = (event: any) => {
            const old_offset = this.offset
            const new_offset_unrounded = unshift(this.client_to_point(event.client.x, event.client.y), start_cursor_position)

            // find closest offset to current position:
            let new_offset = start_offset
            let closest_distance = Infinity
            for (let offset of valid_offsets) {
                const dist = dist2(new_offset_unrounded, offset)
                if (dist < closest_distance) {
                    closest_distance = dist
                    new_offset = offset
                }
            }
            if (old_offset != new_offset) {
                // Do BFS with path reconstruction from old_offset to new_offset
                let queue = [old_offset]
                let parents: { [key: string]: Point | null } = { old_offset: null }
                while (queue.length > 0) {
                    let point: Point = queue.shift()!
                    if (point == new_offset) break
                    for (let neighbor of [shift(point, p(0, 1)), shift(point, p(0, -1)), shift(point, p(1, 0)), shift(point, p(-1, 0))]) {
                        if (neighbor in parents || !valid_offsets.has(neighbor)) continue
                        parents[neighbor] = point
                        queue.push(neighbor)
                    }
                }
                let path = [new_offset]
                let backtrack = parents[new_offset]
                while (backtrack !== old_offset) {
                    path.push(backtrack)
                    backtrack = parents[backtrack]
                }
                this.offset = new_offset
                this.update_walking_offsets(path)
            }
        }
        const end = (event: any) => {
            // Is new blockstate different from before?
            if (start_offset !== this.offset) {
                puzzle.add_to_history(blockstate)
                if (puzzle.won()) {
                    const confetti_number = 24
                    if (puzzle.min_moves !== null && puzzle.history_ix < puzzle.min_moves) jsConfetti.addConfetti({ emojis: ["🐞"], confettiNumber: confetti_number })
                    else if (puzzle.min_moves !== null && puzzle.history_ix === puzzle.min_moves) jsConfetti.addConfetti({ emojis: ["🏆"], confettiNumber: confetti_number })
                    else jsConfetti.addConfetti({ confettiNumber: confetti_number, confettiColors: puzzle.blockstate.blocks.map(b => char_to_color(b.char)) })
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
        blocks.sort((a, b) => a.char.localeCompare(b.char)) // Any sort order suffices

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

    initialise(svg_elem: SVGSVGElement, goal_offsets: Array<Offset | null>) {
        svg_elem.innerHTML = ""
        const width = x(this.max) - x(this.min) + 1
        const height = y(this.max) - y(this.min) + 1
        svg_elem.setAttribute("viewBox", `${x(this.min) - 1} ${y(this.min) - 1} ${width + 1} ${height + 1}`)
        this.svg_elem = svg_elem

        this.bounds.initialise_elem(svg_elem)
        for (let block_ix = 0; block_ix < this.blocks.length; block_ix++) {
            if (goal_offsets[block_ix] !== null) {
                const block = this.blocks[block_ix]
                const shadow_elem = block.construct_elem_path(0)
                const shadow_offset = goal_offsets[block_ix]

                shadow_elem.classList.add("shadowGoal")
                shadow_elem.setAttribute("transform", `translate(${x(shadow_offset)}, ${y(shadow_offset)})`)
                svg_elem.appendChild(shadow_elem)
            }
        }
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
    public history_ix: number
    public history: Array<Array<Offset>>
    public blockstate: Blockstate
    public goal_offsets: Array<Offset | null> // TODO: This is awful

    constructor(start_string: string, goal_string: string, min_moves: number | null = null) {
        this.min_moves = min_moves

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

    add_to_history(blockstate: Blockstate) {
        while (this.history.length > this.history_ix + 1) this.history.pop()
        this.history.push(blockstate.blocks.map(x => x.offset))
        this.history_ix++
        move_counter_elem.textContent = this.history_ix.toString()
    }
    // TODO: Disable buttons when history at limit
    history_forward() {
        if (this.history_ix < this.history.length - 1) {
            this.history_ix++
            this.update_blockstate_from_history()
        }
    }
    history_backward() {
        if (this.history_ix > 0) {
            this.history_ix--
            this.update_blockstate_from_history()
        }
    }
    update_blockstate_from_history() {
        const blockstate = this.blockstate
        const blocks = blockstate.blocks
        const bounds = blockstate.bounds.get_coordinates()
        for (let blockarr_ix = 0; blockarr_ix < blocks.length; blockarr_ix++) {
            const block = blocks[blockarr_ix]
            const old_offset = block.offset
            const new_offset = this.history[this.history_ix][blockarr_ix]
            if (old_offset !== new_offset) {
                const other_blocks_union = union(blocks.filter((b, i) => i != blockarr_ix).map(b => b.get_coordinates()))
                const offset_is_valid = (offset: Point) => {
                    const offset_block = shift_shape(block.shape, offset)
                    if (!is_subset(offset_block, bounds)) return false
                    if (!is_disjoint(offset_block, other_blocks_union)) return false
                    return true
                }
                // Do BFS with path reconstruction from old_offset to new_offset
                let queue = [old_offset]
                let parents: { [key: Point]: Point | null } = {}
                parents[old_offset] = null
                let success = false
                while (queue.length > 0) {
                    let point: Point = queue.shift()!
                    if (point == new_offset) {
                        success = true
                        break
                    }
                    for (let neighbor of [shift(point, p(0, 1)), shift(point, p(0, -1)), shift(point, p(1, 0)), shift(point, p(-1, 0))]) {
                        if (neighbor in parents || !offset_is_valid(neighbor)) continue
                        parents[neighbor] = point
                        queue.push(neighbor)
                    }
                }
                if (success) {
                    let path = []
                    let backtrack = new_offset
                    while (backtrack !== old_offset) {
                        path.push(backtrack)
                        backtrack = parents[backtrack]
                    }
                    block.update_walking_offsets(path)
                } else {
                    // This should never happen.
                    block.update_walking_offsets([new_offset])
                }
            }
        }
        // Finally, update all the offsets for real
        // (We should only ever one offset in the first place, so we could
        //  get away with doing this in the previous loop already, but if
        //  anything ever breaks, this makes it easier to hunt down)
        for (let blockarr_ix = 0; blockarr_ix < blocks.length; blockarr_ix++) {
            const block = blocks[blockarr_ix]
            block.offset = this.history[this.history_ix][blockarr_ix]
        }
        move_counter_elem.textContent = (this.history_ix === this.history.length - 1) ? this.history_ix.toString() : `${this.history_ix}/${this.history.length - 1}`
    }

    initialise(svg_elem: SVGSVGElement) {
        this.blockstate = Blockstate.blockstate_from_string(this.start_string)
        this.blockstate.initialise(svg_elem, this.goal_offsets)
        this.blockstate.make_interactive(this)

        this.history_ix = -1
        this.history = []
        this.add_to_history(this.blockstate)
    }
}

const change_puzzle_dialog = document.getElementById("change-puzzle-dialog") as HTMLDialogElement
const change_puzzle_btn = document.getElementById("change-puzzle-btn") as HTMLButtonElement
const reset_puzzle_btn = document.getElementById("reset-puzzle-btn") as HTMLButtonElement
const puzzle_textarea_start = document.getElementById("puzzle-textarea-start") as HTMLTextAreaElement
const puzzle_textarea_goal = document.getElementById("puzzle-textarea-goal") as HTMLTextAreaElement
const puzzle_submit_btn = document.getElementById("puzzle-submit-btn") as HTMLButtonElement
const puzzle_solve_btn = document.getElementById("puzzle-solve-btn") as HTMLButtonElement
const move_counter_elem = document.getElementById("move-counter") as HTMLSpanElement
const history_forward_btn = document.getElementById("history-forward") as HTMLButtonElement
const history_backward_btn = document.getElementById("history-backward") as HTMLButtonElement
change_puzzle_btn.addEventListener("click", () => {
    puzzle_textarea_goal.value = current_puzzle.goal_string
    puzzle_textarea_start.value = current_puzzle.start_string
    change_puzzle_dialog.showModal()
})
reset_puzzle_btn.addEventListener("click", () => {
    current_puzzle.initialise(svg_puzzle)
})
history_backward_btn.addEventListener("click", () => {
    current_puzzle.history_backward()
})
history_forward_btn.addEventListener("click", () => {
    current_puzzle.history_forward()
})

puzzle_submit_btn.addEventListener("click", e => {
    e.preventDefault()
    const start_string = puzzle_textarea_start.value
    const goal_string = puzzle_textarea_goal.value
    change_puzzle_dialog.close()
    current_puzzle = new Puzzle(start_string, goal_string)
    current_puzzle.initialise(svg_puzzle)
})
puzzle_solve_btn.addEventListener("click", e => {
    e.preventDefault()
    const start_string = puzzle_textarea_start.value
    const goal_string = puzzle_textarea_goal.value
    const solution = solve_puzzle(start_string, goal_string)
    const solution_blockstates = solution.map(s => Blockstate.blockstate_from_string(s))
    change_puzzle_dialog.close()
    current_puzzle = new Puzzle(start_string, goal_string)
    current_puzzle.initialise(svg_puzzle)
    for (let blockstate of solution_blockstates.slice(1)) {
        current_puzzle.add_to_history(blockstate)
    }
    current_puzzle.history_ix = 0
    current_puzzle.update_blockstate_from_history()
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
    puzzle_textarea_start.value = Blockstate.blockstate_from_string(js_puzzle.start).to_string()
    puzzle_textarea_goal.value = Blockstate.blockstate_from_string(js_puzzle.goal).to_string()
})

let current_puzzle: Puzzle = new Puzzle(predefined_puzzles[0].start, predefined_puzzles[0].goal)
current_puzzle.initialise(svg_puzzle)

let timestamp: number
function animate(t: number) {
    const dt = timestamp ? t - timestamp : 0
    timestamp = t
    for (let block of current_puzzle.blockstate.blocks) block.update_translation(dt)
    requestAnimationFrame(animate)
}
requestAnimationFrame(animate)
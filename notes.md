If a puzzle has exactly one block, the puzzle is either:
- Unsolvable, or
- Solvable in 0 moves, or
- Solvable in 1 move.

If a puzzle has two blocks, it can already take an arbitrary number of moves to solve. For instance, consider the following start- and goal-configuration:
```
++++............
+A.+............
++++............

................
..............A.
................
```
This takes 25 moves to solve, because the blocks `A` and `+` need to be moved alternatingly one position to the right. If you don't demand the blocks to be connected, this can even be simplified to a single line each:
```
+A.+............

..............A.
```
Another example, based on a different idea, looks like this:
```
.............................
bbbbbbbbbbbbbbbbbbbbbbbbbbbbb
b   b   b   b   b   b   b   b
bA..b...b...b...b...b...b...b
b b   b   b   b   b   b   b b
b.b...b...b...b...b...b...b.b
bbbbbbbbbbbbbbbbbbbbbbbbbbbbb

.............................
.............................
.   .   .   .   .   .   .   .
...........................A.
. .   .   .   .   .   .   . .
.............................
.............................
```
Here, the block `b` must be moved up and down alternatingly, while block `A` moves two positions to the right at a time.
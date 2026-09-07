# canonaut (Python)

Graph canonicalization and automorphism groups: Python bindings for
[canonaut](https://github.com/um-univie/canonaut), a safe Rust port of
[nauty](https://pallini.di.uniroma1.it/). The PyPI package and the Rust crate
share the name.

## Install

```
pip install canonaut
```

## Usage

Nodes can be any hashable value — ints, strings, tuples, whatever you already
have. No need to relabel your graph to `0..n-1` first.

```python
import canonaut

# Canonize a 4-cycle with arbitrary node labels
r = canonaut.canonize([("a","b"), ("b","c"), ("c","d"), ("d","a")])
print(r.num_orbits)          # 1  (vertex-transitive)
print(r.group_size1, r.group_size2)  # |Aut| ~= group_size1 * 10**group_size2
print(r.canonical_edges)     # edges of the canonical form (0..n-1 index space)

# Vertex colors (partition): vertices with different colors are never mapped
# onto each other by any automorphism, and the canonical form respects it.
r = canonaut.canonize(
    [(0,1), (1,2), (2,3)],
    colors={0: "red", 1: "red", 2: "blue", 3: "blue"},
)
print(r.num_orbits)          # 4 (colors broke the path's end-to-end symmetry)

# Automorphism generators, as {node: image_node} maps
r = canonaut.canonize([("a","b"),("b","c"),("c","d"),("d","a")], generators=True)
print(r.generators)

# Orbits, grouped by original node label
print(r.orbits)               # [['a','b','c','d']] for the 4-cycle

# Isomorphism test — works with different label sets on each side
canonaut.are_isomorphic(
    [("p","q"),("q","r"),("r","s"),("s","p")],
    [(1,2),(2,3),(3,4),(4,1)],
)  # True

# Hashable canonical key for deduplication — equal iff the graphs are isomorphic
key = r.canonical_bytes()

# networkx graphs work directly (duck-typed, no hard dependency)
import networkx as nx
r = canonaut.canonize(nx.petersen_graph())
```

### Isolated vertices

`nodes=` includes vertices with no edges at all:

```python
canonaut.canonize([("a","b")], nodes=["a", "b", "isolated"])
```

### Directed graphs

```python
canonaut.canonize(edges, directed=True)
```

Detected automatically for graph-like objects exposing `.is_directed()` (e.g.
`networkx.DiGraph`).

### Fixed vertex colors and canonical keys

Color values are fixed labels: an isomorphism must map each color to the same
value in the other graph. Supported labels are integers, floats (except NaN),
strings, bytes, `None`, and tuples/frozensets containing these. Python-equal
numeric values such as `1` and `1.0` denote the same color. `None` and missing
dictionary entries share the implicit uncolored label. Map custom color objects
to stable labels before calling `canonize`.

`canonical_bytes()` includes graph order, directedness, and fixed colors as well
as canonical adjacency. Its versioned format has changed; regenerate previously
stored keys. The high-level and raw APIs have separate color encodings, so compare
keys produced through the same API. `canonical_edges` alone does not encode colors.

## Low-level API

`canonaut.raw` exposes the compiled extension directly: `0..n-1` integer
vertex indices only, no node-label translation. Use it if you already have
integer-labeled graphs and want to skip that overhead:

```python
r = canonaut.raw.canonize(4, [(0,1),(1,2),(2,3)], colors=[0,0,1,1])
```

## `CanonResult`

| Attribute | Type | Meaning |
|---|---|---|
| `labeling` | `list` | Original nodes in canonical order — `labeling[i]` is the node placed at canonical position `i` |
| `orbits` | `list[list]` | Automorphism-group orbits, as groups of original nodes |
| `canonical_edges` | `list[tuple[int,int]]` | Edges of the canonical form, in canonical-position (`0..n-1`) space |
| `group_size1`, `group_size2` | `float`, `int` | `\|Aut(G)\| ~= group_size1 * 10**group_size2` |
| `num_orbits`, `num_generators` | `int` | |
| `generators` | `list[dict]` | Automorphism generators as `{node: image_node}` maps (only when `generators=True`) |
| `.canonical_bytes()` | `bytes` | Hashable canonical key |

"""
canonaut — ergonomic Python bindings for canonaut (a safe Rust port of nauty).

High-level API (this module): works with arbitrary hashable node labels
(ints, strings, tuples, ...), dict-based vertex colors, and duck-types
networkx-like graph objects. See `canonize()` and `are_isomorphic()`.

Low-level API (`canonaut.raw`): the compiled extension's raw functions,
operating on 0..n-1 integer vertex indices only. Use these if you already
have integer-labeled graphs and want to skip the node-translation overhead.
"""

from . import _canonaut as raw
from ._canonaut import CanonResult as _RawCanonResult

__all__ = ["canonize", "are_isomorphic", "CanonResult", "raw"]

_UNCOLORED = object()


def _graph_like(g):
    """Duck-type networkx-like graph objects: anything with .edges()/.nodes()."""
    if isinstance(g, (list, tuple)):
        return None
    if not (hasattr(g, "edges") and hasattr(g, "nodes")):
        return None
    is_multigraph = getattr(g, "is_multigraph", None)
    if callable(is_multigraph) and is_multigraph():
        raise TypeError("multigraphs are not supported; convert explicitly to a simple graph")
    is_directed = getattr(g, "is_directed", None)
    directed = is_directed() if callable(is_directed) else None
    return list(g.edges()), list(g.nodes()), directed


def _resolve_nodes(edges, nodes):
    """Node list in a deterministic order: `nodes` first (if given), then any
    further nodes seen only in `edges`, in first-appearance order."""
    resolved = list(dict.fromkeys(nodes)) if nodes is not None else []
    seen = set(resolved)
    for u, v in edges:
        for node in (u, v):
            if node not in seen:
                resolved.append(node)
                seen.add(node)
    return resolved


def _pack_color_parts(tag, parts):
    return tag + b"".join(len(part).to_bytes(8, "big") + part for part in parts)


def _encode_color(value):
    """Stable, injective tokens for supported color labels, respecting 1 == 1.0."""
    if value is _UNCOLORED or value is None:
        return b"u"
    if isinstance(value, (int, float)):
        if isinstance(value, float):
            if value != value:
                raise ValueError("NaN is not a valid vertex color")
            if value in (float("inf"), float("-inf")):
                return b"n+inf" if value > 0 else b"n-inf"
            numerator, denominator = value.as_integer_ratio()
        else:
            numerator, denominator = int(value), 1
        return f"n{numerator}/{denominator}".encode("ascii")
    if isinstance(value, str):
        return b"s" + value.encode("utf-8", "surrogatepass")
    if isinstance(value, bytes):
        return b"b" + value
    if isinstance(value, tuple):
        return _pack_color_parts(b"t", [_encode_color(item) for item in value])
    if isinstance(value, frozenset):
        return _pack_color_parts(b"f", sorted({_encode_color(item) for item in value}))
    raise TypeError(
        "vertex colors must be numbers, strings, bytes, None, or tuples/frozensets "
        "of these; map custom objects to stable color labels first"
    )


def _resolve_colors(colors, node_list):
    if colors is None:
        return None, None
    if hasattr(colors, "get"):
        raw_values = [colors.get(node, _UNCOLORED) for node in node_list]
    else:
        raw_values = list(colors)
        if len(raw_values) != len(node_list):
            raise ValueError(
                f"colors has length {len(raw_values)} but there are {len(node_list)} nodes"
            )
    if not raw_values:
        return None, None
    raw_values = [_UNCOLORED if value is None else value for value in raw_values]
    if all(value is _UNCOLORED for value in raw_values):
        return None, None
    tokens = {value: _encode_color(value) for value in dict.fromkeys(raw_values)}
    palette = sorted(set(tokens.values()))
    ids = {token: index for index, token in enumerate(palette)}
    return [ids[tokens[value]] for value in raw_values], palette


class CanonResult:
    """Ergonomic wrapper around the raw extension's ``CanonResult``.

    Attributes
    ----------
    labeling : list
        Original nodes in canonical order: ``labeling[i]`` is the node placed
        at canonical position ``i``.
    orbits : list[list]
        Automorphism-group orbits, as groups of original nodes — nodes that
        some automorphism maps onto each other share a group.
    canonical_edges : list[tuple[int, int]]
        Edges of the canonical form, in canonical-position space (``0..n-1``)
        — a structural invariant, independent of your node labels.
    group_size1, group_size2 : float, int
        ``|Aut(G)| ~= group_size1 * 10**group_size2``.
    num_generators, num_orbits : int
    generators : list[dict]
        Automorphism generators as ``{node: image_node}`` maps (populated
        only when ``canonize(..., generators=True)``).
    """

    __slots__ = (
        "labeling",
        "orbits",
        "canonical_edges",
        "group_size1",
        "group_size2",
        "num_generators",
        "num_orbits",
        "generators",
        "_raw",
        "_color_palette",
    )

    def __init__(self, raw_result: _RawCanonResult, node_list: list, color_palette=None):
        n = len(node_list)
        self._raw = raw_result
        self._color_palette = color_palette
        # PyO3 Vec getters copy the list; fetch each once to keep conversion linear.
        raw_labeling = raw_result.labeling
        raw_orbits = raw_result.orbits
        self.labeling = [node_list[vertex] for vertex in raw_labeling]
        groups: dict = {}
        for v in range(n):
            groups.setdefault(raw_orbits[v], []).append(node_list[v])
        self.orbits = list(groups.values())
        self.canonical_edges = raw_result.canonical_edges
        self.group_size1 = raw_result.group_size1
        self.group_size2 = raw_result.group_size2
        self.num_generators = raw_result.num_generators
        self.num_orbits = raw_result.num_orbits
        self.generators = [
            {node_list[v]: node_list[perm[v]] for v in range(n)}
            for perm in raw_result.generators
        ]

    def canonical_bytes(self) -> bytes:
        """Hashable canonical key: equal iff the two graphs are isomorphic."""
        key = self._raw.canonical_bytes()
        if self._color_palette is not None:
            key += _pack_color_parts(b"COL\x01", self._color_palette)
        return key

    def __repr__(self) -> str:
        return (
            f"CanonResult(num_orbits={self.num_orbits}, "
            f"|Aut|~={self.group_size1:.4g}e{self.group_size2}, "
            f"generators={len(self.generators)})"
        )


def canonize(edges, *, nodes=None, colors=None, directed=False, generators=False):
    """
    Compute the canonical form of a graph with arbitrary hashable node labels.

    Parameters
    ----------
    edges : iterable of (node, node) pairs, or a graph-like object
        Edge list, or any object exposing ``.edges()``/``.nodes()`` (e.g. a
        ``networkx.Graph``/``DiGraph``). Nodes can be any hashable value.
        Multigraph objects are rejected; convert them explicitly to simple
        graphs if edge multiplicities should be discarded.
    nodes : iterable, optional
        All nodes, in the order you want them considered — only needed to
        include isolated vertices with no edges, or to fix ordering.
    colors : dict or sequence, optional
        Vertex coloring (partition). As a dict: ``{node: color}`` — nodes not
        present get an implicit shared color distinct from any color you
        supplied. As a sequence: one color per node, aligned with the
        resolved node order (see `nodes`). Colors may be numbers, strings,
        bytes, None, or tuples/frozensets of these;
        color values are fixed labels across graphs. Vertices with different
        colors are never mapped onto each other by any automorphism, and the
        canonical form respects the coloring.
    directed : bool, optional
        Treat edges as one-directional. Overridden if `edges` is a
        graph-like object exposing ``.is_directed()``.
    generators : bool, optional
        Also collect automorphism generators. Slightly slower.

    Returns
    -------
    CanonResult

    Examples
    --------
    >>> import canonaut
    >>> r = canonaut.canonize([("a","b"),("b","c"),("c","d"),("d","a")])
    >>> r.num_orbits
    1
    >>> r2 = canonaut.canonize(
    ...     [(0,1),(1,2),(2,3)],
    ...     colors={0: "red", 1: "red", 2: "blue", 3: "blue"},
    ... )
    >>> r2.num_orbits
    4
    """
    like = _graph_like(edges)
    if like is not None:
        edges, inferred_nodes, inferred_directed = like
        if nodes is None:
            nodes = inferred_nodes
        if inferred_directed is not None:
            directed = inferred_directed

    edges = list(edges)
    node_list = _resolve_nodes(edges, nodes)
    index_of = {node: i for i, node in enumerate(node_list)}
    n = len(node_list)

    edge_idx = [(index_of[u], index_of[v]) for u, v in edges]
    color_idx, color_palette = _resolve_colors(colors, node_list)

    raw_result = raw.canonize(
        n, edge_idx, directed=directed, generators=generators, colors=color_idx
    )
    return CanonResult(raw_result, node_list, color_palette)


def are_isomorphic(edges1, edges2, *, nodes1=None, nodes2=None, directed=False):
    """
    Test whether two graphs (arbitrary hashable node labels) are isomorphic.

    Parameters
    ----------
    edges1, edges2 : iterable of (node, node) pairs, or graph-like objects
        Multigraph objects are rejected.
    nodes1, nodes2 : iterable, optional
        Full node sets — only needed for isolated vertices, or explicit order.
    directed : bool, optional
        Overridden per-graph if that graph is a graph-like object exposing
        ``.is_directed()``.

    Returns
    -------
    bool
    """
    directed1 = directed2 = directed
    like1 = _graph_like(edges1)
    if like1 is not None:
        edges1, inferred_nodes1, inferred_directed1 = like1
        if nodes1 is None:
            nodes1 = inferred_nodes1
        if inferred_directed1 is not None:
            directed1 = inferred_directed1

    like2 = _graph_like(edges2)
    if like2 is not None:
        edges2, inferred_nodes2, inferred_directed2 = like2
        if nodes2 is None:
            nodes2 = inferred_nodes2
        if inferred_directed2 is not None:
            directed2 = inferred_directed2

    if directed1 != directed2:
        return False

    edges1 = list(edges1)
    edges2 = list(edges2)
    node_list1 = _resolve_nodes(edges1, nodes1)
    node_list2 = _resolve_nodes(edges2, nodes2)

    if len(node_list1) != len(node_list2):
        return False

    index_of1 = {node: i for i, node in enumerate(node_list1)}
    index_of2 = {node: i for i, node in enumerate(node_list2)}
    edge_idx1 = [(index_of1[u], index_of1[v]) for u, v in edges1]
    edge_idx2 = [(index_of2[u], index_of2[v]) for u, v in edges2]

    return raw.are_isomorphic(
        len(node_list1), edge_idx1, len(node_list2), edge_idx2, directed=directed1
    )

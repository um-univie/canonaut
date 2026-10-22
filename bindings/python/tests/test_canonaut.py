"""
Tests for the ergonomic canonaut Python API.

Run directly (no pytest required):
    .venv-testdata/bin/python3 -m unittest bindings/python/tests/test_canonaut.py -v
or, if pytest is installed:
    pytest bindings/python/tests/
"""

import unittest
import itertools

import canonaut


class TestCanonize(unittest.TestCase):
    def test_arbitrary_string_labels(self):
        r = canonaut.canonize([("a", "b"), ("b", "c"), ("c", "d"), ("d", "a")])
        self.assertEqual(r.num_orbits, 1)
        self.assertEqual(r.group_size1, 8.0)
        self.assertEqual(set(r.labeling), {"a", "b", "c", "d"})

    def test_uncolored_p4_has_two_orbits(self):
        r = canonaut.canonize([(0, 1), (1, 2), (2, 3)])
        self.assertEqual(r.num_orbits, 2)
        self.assertEqual(r.group_size1, 2.0)
        orbit_sets = {frozenset(g) for g in r.orbits}
        self.assertEqual(orbit_sets, {frozenset({0, 3}), frozenset({1, 2})})

    def test_dict_colors_break_symmetry(self):
        r = canonaut.canonize(
            [(0, 1), (1, 2), (2, 3)],
            colors={0: "red", 1: "red", 2: "blue", 3: "blue"},
        )
        self.assertEqual(r.num_orbits, 4)
        self.assertEqual(r.group_size1, 1.0)

    def test_sequence_colors_aligned_to_nodes(self):
        r = canonaut.canonize(
            [(0, 1), (1, 2), (2, 3)],
            nodes=[0, 1, 2, 3],
            colors=[0, 0, 1, 1],
        )
        self.assertEqual(r.num_orbits, 4)

    def test_colors_missing_nodes_share_implicit_bucket(self):
        # P4's only nontrivial automorphism is the reversal 0<->3, 1<->2, which
        # maps {0,3} onto itself and {1,2} onto itself. Coloring exactly {0,3}
        # and leaving 1,2 as an implicit shared bucket must NOT break that
        # automorphism, since it still respects both classes as sets.
        colored = canonaut.canonize(
            [(0, 1), (1, 2), (2, 3)], colors={0: "end", 3: "end"}
        )
        uncolored = canonaut.canonize([(0, 1), (1, 2), (2, 3)])
        # The automorphism group is unchanged, but colors are part of the key.
        self.assertNotEqual(colored.canonical_bytes(), uncolored.canonical_bytes())
        self.assertEqual(colored.num_orbits, uncolored.num_orbits)

    def test_colors_on_one_vertex_breaks_symmetry(self):
        # Coloring only vertex 0 forces it out of {0,3}'s pairing, forbidding
        # the reversal automorphism entirely.
        r = canonaut.canonize([(0, 1), (1, 2), (2, 3)], colors={0: "special"})
        self.assertEqual(r.num_orbits, 4)

    def test_colors_wrong_sequence_length_raises(self):
        with self.assertRaises(ValueError):
            canonaut.canonize([(0, 1), (1, 2)], colors=[0, 1])

    def test_isolated_node_via_nodes(self):
        r = canonaut.canonize([("a", "b")], nodes=["a", "b", "isolated"])
        orbit_sets = {frozenset(g) for g in r.orbits}
        self.assertIn(frozenset({"isolated"}), orbit_sets)

    def test_generators_are_node_to_node_maps(self):
        r = canonaut.canonize(
            [("a", "b"), ("b", "c"), ("c", "d"), ("d", "a")], generators=True
        )
        self.assertEqual(len(r.generators), r.num_generators)
        for gen in r.generators:
            self.assertEqual(set(gen.keys()), {"a", "b", "c", "d"})
            self.assertEqual(set(gen.values()), {"a", "b", "c", "d"})

    def test_canonical_bytes_is_isomorphism_invariant(self):
        k1 = canonaut.canonize(
            [("a", "b"), ("b", "c"), ("c", "d"), ("d", "a")]
        ).canonical_bytes()
        k2 = canonaut.canonize([(1, 2), (2, 3), (3, 4), (4, 1)]).canonical_bytes()
        self.assertEqual(k1, k2)

    def test_networkx_duck_typing(self):
        nx = _import_networkx_or_skip(self)
        g = nx.cycle_graph(["x", "y", "z", "w"])
        r = canonaut.canonize(g)
        self.assertEqual(r.num_orbits, 1)
        self.assertEqual(r.group_size1, 8.0)


class TestCorrectnessRegressions(unittest.TestCase):
    def test_multigraph_objects_are_rejected_before_losing_multiplicity(self):
        nx = _import_networkx_or_skip(self)
        for cls in [nx.MultiGraph, nx.MultiDiGraph]:
            a, b = cls(), cls()
            a.add_edges_from([(0, 1), (0, 1)])
            b.add_edge(0, 1)
            for graph in [a, b, cls()]:
                with self.assertRaisesRegex(TypeError, "multigraphs are not supported"):
                    canonaut.canonize(graph)
            for left, right in [(a, b), (b, a), (nx.Graph(), a)]:
                with self.assertRaisesRegex(TypeError, "multigraphs are not supported"):
                    canonaut.are_isomorphic(left, right)

    def test_directed_path_has_no_nontrivial_automorphism(self):
        r = canonaut.raw.canonize(4, [(0, 2), (1, 0)], directed=True)
        self.assertEqual(r.group_size1, 1.0)
        self.assertEqual(r.num_orbits, 4)

    def test_all_three_vertex_digraph_groups_against_brute_force(self):
        slots = list(itertools.product(range(3), repeat=2))
        permutations = list(itertools.permutations(range(3)))
        for mask in range(1 << len(slots)):
            edges = {edge for i, edge in enumerate(slots) if mask & (1 << i)}
            automorphisms = [p for p in permutations if {(p[u], p[v]) for u, v in edges} == edges]
            r = canonaut.raw.canonize(3, list(edges), directed=True, generators=True)
            self.assertEqual(r.group_size1, len(automorphisms), (mask, edges))
            for p in r.generators:
                self.assertIn(tuple(p), automorphisms)

    def test_colored_keys_do_not_depend_on_node_order(self):
        edges = [(0, 1), (1, 2)]
        colors = {0: "red", 1: "blue", 2: "blue"}
        keys = {canonaut.canonize(edges, nodes=order, colors=colors).canonical_bytes()
                for order in itertools.permutations(range(3))}
        self.assertEqual(len(keys), 1)

    def test_fixed_color_labels_and_partitions_are_in_the_key(self):
        keys = {canonaut.canonize([], nodes=[0, 1], colors=colors).canonical_bytes()
                for colors in [None, [0, 0], [0, 1], [1, 1], ["red", "red"]]}
        self.assertEqual(len(keys), 5)
        raw_keys = {canonaut.raw.canonize(2, [], colors=colors).canonical_bytes()
                    for colors in [None, [0, 0], [0, 1], [1, 1]]}
        self.assertEqual(len(raw_keys), 4)

    def test_colored_keys_are_invariant_under_relabeling(self):
        a = canonaut.canonize([(0, 1), (1, 2)], colors={0: "red", 1: "blue", 2: "blue"})
        b = canonaut.canonize([("z", "x"), ("y", "z")],
                             colors={"y": "red", "z": "blue", "x": "blue"})
        self.assertEqual(a.canonical_bytes(), b.canonical_bytes())

    def test_stable_composite_and_numeric_colors(self):
        for a, b in [(1, 1.0), (True, 1), (0.0, -0.0),
                     (("C", frozenset([1, 2])), ("C", frozenset([2.0, 1.0])))]:
            self.assertEqual(canonaut.canonize([], nodes=[0], colors=[a]).canonical_bytes(),
                             canonaut.canonize([], nodes=[0], colors=[b]).canonical_bytes())
        with self.assertRaises(TypeError):
            canonaut.canonize([], nodes=[0], colors=[object()])
        with self.assertRaises(ValueError):
            canonaut.canonize([], nodes=[0], colors=[float("nan")])

    def test_directedness_is_in_the_key(self):
        self.assertNotEqual(canonaut.canonize([], nodes=[0]).canonical_bytes(),
                            canonaut.canonize([], nodes=[0], directed=True).canonical_bytes())

    def test_empty_graphs(self):
        r = canonaut.canonize([])
        self.assertEqual(r.labeling, [])
        self.assertEqual(r.group_size1, 1.0)
        self.assertEqual(r.num_orbits, 0)
        self.assertTrue(canonaut.are_isomorphic([], []))
        self.assertEqual(r.canonical_bytes(), canonaut.canonize([], colors=[]).canonical_bytes())

    def test_implicit_uncolored_labels_match_an_uncolored_graph(self):
        expected = canonaut.canonize([(0, 1)]).canonical_bytes()
        for colors in [{}, [None, None]]:
            self.assertEqual(expected, canonaut.canonize([(0, 1)], colors=colors).canonical_bytes())

    def test_raw_invalid_vertices_are_value_errors(self):
        for directed in [False, True]:
            with self.assertRaises(ValueError):
                canonaut.raw.canonize(3, [(0, 3)], directed=directed)
            with self.assertRaises(ValueError):
                canonaut.raw.are_isomorphic(3, [(3, 0)], 3, [(0, 1)], directed=directed)

    def test_duplicate_nodes_and_edges_do_not_change_a_simple_graph(self):
        self.assertEqual(len(canonaut.canonize([(0, 1)], nodes=[0, 0, 1]).labeling), 2)
        self.assertTrue(canonaut.are_isomorphic([(0, 1)], [(0, 1), (1, 0)]))

    def test_mixed_graph_kinds_are_not_isomorphic(self):
        nx = _import_networkx_or_skip(self)
        a, b = nx.Graph([(0, 1)]), nx.DiGraph([(0, 1)])
        self.assertFalse(canonaut.are_isomorphic(a, b))
        self.assertFalse(canonaut.are_isomorphic(b, a))

    def test_undirected_loops_preserve_the_correct_group(self):
        r = canonaut.raw.canonize(3, [(0, 0)])
        self.assertEqual(r.group_size1, 2.0)
        self.assertEqual(r.num_orbits, 2)

    def test_networkx_directed_is_detected(self):
        nx = _import_networkx_or_skip(self)
        g = nx.DiGraph([(0, 1), (1, 2), (2, 0)])
        r = canonaut.canonize(g)
        self.assertEqual(r.num_orbits, 1)
        self.assertEqual(r.group_size1, 3.0)


class TestAreIsomorphic(unittest.TestCase):
    def test_string_labels(self):
        self.assertTrue(
            canonaut.are_isomorphic(
                [("p", "q"), ("q", "r"), ("r", "s"), ("s", "p")],
                [(1, 2), (2, 3), (3, 4), (4, 1)],
            )
        )

    def test_non_isomorphic_graphs(self):
        # Path vs cycle on the same vertex count are never isomorphic.
        self.assertFalse(
            canonaut.are_isomorphic(
                [(0, 1), (1, 2), (2, 3)],
                [(0, 1), (1, 2), (2, 3), (3, 0)],
            )
        )

    def test_mismatched_vertex_count_is_false(self):
        self.assertFalse(
            canonaut.are_isomorphic([(0, 1)], [(0, 1), (1, 2)])
        )


def _import_networkx_or_skip(testcase):
    try:
        import networkx as nx
    except ImportError:
        testcase.skipTest("networkx not installed")
    return nx


if __name__ == "__main__":
    unittest.main()

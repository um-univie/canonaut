#include <nauty.h>

/*
 * Thin wrapper around densenauty for benchmarking against canonaut.
 *
 * g      : input graph  (n*m setwords, MSB-0 adjacency bit-matrix)
 * lab    : pre-allocated int[n]   (filled by nauty: canonical labeling)
 * ptn    : pre-allocated int[n]   (re-initialised by nauty: defaultptn=TRUE)
 * orbits : pre-allocated int[n]   (filled by nauty: orbit indices)
 * canong : pre-allocated graph[n*m] (canonical graph output)
 * m      : setwords per row  (= ceil(n/WORDSIZE))
 * n      : number of vertices
 */
void c_nauty_dense_canonize(
    graph *g,
    int *lab, int *ptn, int *orbits,
    graph *canong,
    int m, int n
) {
    DEFAULTOPTIONS_GRAPH(opts);
    opts.getcanon = TRUE;
    statsblk stats;
    densenauty(g, lab, ptn, orbits, &opts, &stats, m, n, canong);
}

/*
 * Same call, but also reporting nauty's search statistics, so the port's
 * search tree can be compared against the reference's node for node.
 */
void c_nauty_dense_canonize_stats(
    graph *g,
    int *lab, int *ptn, int *orbits,
    graph *canong,
    int m, int n,
    double *group_size1, int *group_size2, int *num_orbits, int *num_generators,
    long *num_nodes, long *num_bad_leaves, int *max_level, long *canon_updates
) {
    DEFAULTOPTIONS_GRAPH(opts);
    opts.getcanon = TRUE;
    statsblk stats;
    densenauty(g, lab, ptn, orbits, &opts, &stats, m, n, canong);

    *group_size1 = stats.grpsize1;
    *group_size2 = stats.grpsize2;
    *num_orbits = stats.numorbits;
    *num_generators = stats.numgenerators;
    *num_nodes = stats.numnodes;
    *num_bad_leaves = stats.numbadleaves;
    *max_level = stats.maxlevel;
    *canon_updates = stats.canupdates;
}

/* Correctness oracle for colored graphs and digraphs. lab/ptn are initialized
 * by the caller when colored is true. Keep this independent of Rust's setup. */
void c_nauty_reference(
    graph *g, int *lab, int *ptn, int *orbits, graph *canong,
    int m, int n, int directed, int colored,
    double *group_size1, int *group_size2
) {
    DEFAULTOPTIONS_GRAPH(opts);
    opts.getcanon = TRUE;
    opts.digraph = directed;
    opts.defaultptn = !colored;
    if (directed) opts.tc_level = 0;
    statsblk stats;
    densenauty(g, lab, ptn, orbits, &opts, &stats, m, n, canong);
    *group_size1 = stats.grpsize1;
    *group_size2 = stats.grpsize2;
}

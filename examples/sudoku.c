// From Laurent Fousse - sudoku graph example, released into the public domain

#include "vn_graph.h"

int main (int argc, char *argv[]) {
    int n, i, j, k, l, si, sj;
    graph_t g;
    if (argc > 1) {
	n = atoi (argv[1]);
    } else
	n = 3;

    g = graph_new (n * n * n * n);
    /* lines and columns conditions */
    for (k = 0; k < n * n; k++) {
	for (i = 0; i < n * n; i++) {
	    for (j = i + 1; j < n * n; j++) {
		graph_add_edge (g, k * n * n + i, k * n * n + j);
		graph_add_edge (g, i * n * n + k, j * n * n + k);
	    }
	}
    }

    /* squares conditions */
    for (si = 0; si < n; si++) {
	for (sj = 0; sj < n; sj++) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
		    for (k = 0; k < n; k++) {
			for (l = 0; l < n; l++) {
			    if ((i != k) || (j != l)) {
				graph_add_edge (g, ((n * si) + i) * n * n + sj + j,
						   ((n * si) + k) * n * n + sj + l);
			    }
			}
		    }
		}
	    }
	}
    }
    printf ("chi = %d\n", graph_chromatic_number (g, 0));
    graph_make_dotfile_colored(g, "sudoku.dot");
    printf ("omega = %d\n", graph_clique_number (g));
    graph_to_theta(g,"sudoku.theta");
    graph_clear (g);
    return 0;
}

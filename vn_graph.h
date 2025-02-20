// KMB 2006 Dec 13

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <sys/times.h>
#define NDEBUG
#include <assert.h>

#define VERY_NAUTY_VERSION 1.1

//! useful macros for the user
#define degree(g,i)      (g->d[i])
#define visited(g,i)     (g->v[i])
#define color(g,i)       (g->c[i])
#define cluster(g,i)     (g->l[i])
#define neighbour(g,i,j) (g->a[i][j])
#define nnodes(g)        (g->nnodes)
#define nedges(g)        (g->nedges)

//! type of a node
typedef unsigned int node_t;

/*
struct _graph {
  node_t **a; // adjacency lists
  int  *b; // allocated length of adjacency list size_t
  int  *d; // degree (d[i]<=b[i]) node_t
  int  *v; // visited char
  int  *c; // color
  int  *l; // cluster number
  unsigned int nnodes;
  unsigned int nedges;
};
*/

struct _graph {
  node_t **a; // adjacency lists
  node_t  *d; // degree (d[i]<=b[i])
  size_t  *b; // allocated length of adjacency list
  char    *v; // visited 
  int     *c; // color
  int     *l; // cluster number
  unsigned int nnodes;
  unsigned int nedges;
};

//! version string
const char* graph_get_version(void);

//! type of a graph
typedef struct _graph* graph_t;
//! make a new empty graph with n nodes
graph_t graph_new(unsigned int n);
//! set initial size of adjacency lists
void set_adj_list_initial_size(size_t sz);
//! destroy a graph
void graph_clear(graph_t g);
//! set to the empty graph
void graph_empty(graph_t g);
//! set to the complete graph
void graph_complete(graph_t g);
//! add undirected edge (i,j) to the graph
void graph_add_edge(graph_t g, node_t i, node_t j);
//! 1 if edge i--j is present in g
int graph_has_edge(graph_t g, node_t i, node_t j);
//! delete edge u--v from the graph (if it is present)
// return 1 if edge was found and deleted, else 0
int graph_del_edge(graph_t g, node_t u, node_t v);
//! in-place local_complementation at node i; i.e. subgraph N_i(g) is complemented, where N_i(g) is the set of neighbours of i
void graph_local_complement(graph_t g, node_t i);

void graph_add_node(graph_t g);
void graph_show(graph_t g);
int  graph_greedy_color(graph_t g, int perm[]);
void graph_gnp(graph_t g, double p);
void graph_gnm(graph_t g, unsigned long m);
void graph_grg(graph_t g, double r);
void graph_grg_torus(graph_t g, double r);
void graph_lognormal_grg_torus(graph_t g, double r, double alpha);
int graph_grg_torus_iterator(unsigned int n, double r, int f(graph_t,int,void*), void* cd);
int graph_gnm_iterator(unsigned int n, unsigned long m, int f(graph_t,int,void*), void* cd);
int graph_gnp_iterator(unsigned int n, double p, int f(graph_t,int,void*), void* cd);
int graph_random_line_graph_iterator(unsigned int nlines, int f(graph_t,int,void*), void* cd);
int graph_random_line_graph_poisson_iterator(double tau, int f(graph_t,int,void*), void* cd);
int  graph_ncolors(graph_t g);
int  graph_chromatic_number(graph_t g, clock_t timeout);
int  graph_edge_chromatic_number(graph_t g, clock_t timeout);
int  graph_clique_number(graph_t g);
int  graph_check_coloring(graph_t g);
int  graph_sequential_color(graph_t g, int perm[], int ub);
int  graph_sequential_color_repeat(graph_t g, int n);
void graph_make_dotfile(graph_t g, char* fn);
void graph_make_dotfile_colored(graph_t g, char* fn);
graph_t graph_line_graph(graph_t g);
void graph_to_dimacs(graph_t g, char* fn);
void graph_to_theta(graph_t g, char* fn);
int  graph_node_degree(graph_t g, node_t i);
int  graph_min_degree(graph_t g);
int  graph_max_degree(graph_t g);
double graph_mean_degree(graph_t g);
int  graph_nclusters(graph_t g);
int* graph_cluster_sizes(graph_t g);
int  graph_max_cluster(graph_t g);
int  graph_connected(graph_t g);

struct _histogram {
  unsigned long n; // total count
  int a; // allocation 0..a-1
  double *h; // bins - double for portability, but actually only unsigned longs used
  double *m; // means
  double *ss; // sums-of-squares
  char *lbl;
};

//! histogram functions
//! type of a histogram
typedef struct _histogram* histogram_t;

histogram_t histogram_new(char*);
void histogram_add(histogram_t h, int i);
void histogram_empty(histogram_t h);
void histogram_clear(histogram_t h);
void histogram_empty(histogram_t h); // TODO
int histogram_max(histogram_t h);
void histogram_show(histogram_t h);
void histogram_write(FILE*,histogram_t h);
int  histogram_biggest_done(const histogram_t h, double eps);
int histogram_biggest2_done(const histogram_t h, double eps);
int histogram_mean_done(const histogram_t h, double eps);
double histogram_quantile(const histogram_t h, double y);
double histogram_mean(const histogram_t h);
double histogram_variance(const histogram_t h);
unsigned int histogram_mode(histogram_t h);
histogram_t graph_showg_reader(FILE* f, int (graph_t), char* lbl);
histogram_t graph_geng_reader(FILE* f, int (graph_t), char* lbl);

#define histogram_bin(x,i) ((unsigned long)x->h[i])
#define histogram_rel_bin(x,i) (x->h[i]/(double)x->n)

//! geng functions
char* geng_getline(FILE *f);
graph_t geng_stringtograph(char *s);
int geng_graphsize(char *s);

//! RNGs
double poisson_deviate(double lambda);
double geometric_deviate(double p);
double binomial_deviate(double pp, int n);

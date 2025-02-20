// KMB 2006 Dec 04
// gcc -Wall -O3  -I.. -fomit-frame-pointer gnm_generator.c -o gnm_generator -L.. -lvn_graph -lm -D TEST_GNM_GENERATOR

#include <stdlib.h>
#include "vn_graph.h"

graph_t graph_gnm(int n, int m) {
/*
  Return the random graph G_{n,m}.
 
  Gives a graph picked randomly out of the set of all graphs
  with n nodes and m edges.
 
  :Parameters:
   - n: the number of nodes
   - m: the number of edges
 
  Algorithm by Keith M. Briggs - python origianl Mar 31, 2006.
  Inspired by Knuth's Algorithm S (Selection sampling technique),
  in section 3.4.2 of The Art of Computer Programming by Donald E. Knuth
  Volume 2 / Seminumerical algorithms
  Third Edition, Addison-Wesley, 1997.
*/
  int u=0,v=1,t=0,k=0,mmax=n*(n-1)/2;
  graph_t g=graph_new(n);
  if (m>=mmax) {
    graph_complete(g);
    return g;
  }
  if (n==1) return g;
  while (1) {
    if ((mmax-t)*(rand()/(1.0+RAND_MAX))<m-k) {
      graph_add_edge(g,u,v);
      if (++k==m) return g;
    }
    t++;
    v++;
    if (v==n) {
      u++;
      v=u+1;
    }
  }
}

#ifdef TEST_GNM_GENERATOR
int main() {
  int n=20, m=4;
  graph_t g=graph_gnm(n,m);
  graph_show(g);
  printf("%d edges\n",nedges(g));
  return 0;
}
#endif

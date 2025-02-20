#include <string.h>
#include "../vn_graph.h"

void graph_gnp_slow(graph_t g, double p);

/*
//! delete edge u--v from the graph (if it is present)
// return 1 if edge was found and deleted, else 0
int graph_del_edge(graph_t g, node_t u, node_t v) {
  int i,j;
  for (i=0; i<g->d[u]; i++) if (g->a[u][i]==v) {
    memmove(&g->a[u][i],&g->a[u][i+1],(g->d[u]-1-i)*sizeof(node_t));
    g->d[u]--;
    for (j=0; j<g->d[v]; j++) if (g->a[v][j]==u) {
      memmove(&g->a[v][j],&g->a[v][j+1],(g->d[v]-1-j)*sizeof(node_t));
      g->d[v]--;
      break;
    }
    g->nedges--;
    return 1;
  }
  return 0;
}

int graph_has_edge(graph_t g, node_t i, node_t j) {
  int k;
  for (k=0; k<g->d[i]; k++) if (g->a[i][k]==j) return 1;
  return 0;
}

void local_complement(graph_t g, node_t i) {
  // in-place local_complementation at node i
  // i.e. subgraph N_i(g) is complemented, where N_i(g) is the set of
  // neighbours of i
  const int verbose=0;
  int a,d,nd=0,na=0,n=g->d[i]*(g->d[i]-1)/2;
  node_t j,k,nj,nk,addj[n],addk[n],delj[n],delk[n];
  for (j=0; j<degree(g,i); j++) {
    nj=g->a[i][j];
    for (k=0; k<degree(g,i); k++) {
      nk=g->a[i][k];
      if (nj<nk) { 
        if (graph_has_edge(g,nj,nk)) { 
          delj[nd]=nj; delk[nd]=nk; nd++;
        } else { 
          addj[na]=nj; addk[na]=nk; na++;
        }
      }
    }
  }
  if (verbose) {
    printf("N_%d(g): ",i);
    for (d=0; d<nd; d++) printf("%d--%d ",delj[d],delk[d]);
    printf("\nna=%d nd=%d\n",na,nd);
  }
  for (a=0; a<na; a++) graph_add_edge(g,addj[a],addk[a]);
  for (d=0; d<nd; d++) graph_del_edge(g,delj[d],delk[d]);
}
*/

int main() {
  int i,c;
  graph_t g=graph_new(100);
  for (i=0; i<10000; i++) {
    graph_gnp(g,0.2);
    //printf("g\n"); graph_show(g);
    c=graph_chromatic_number(g,10);
    //c=graph_clique_number(g);
    //printf("chi(g)=%d\n",c);
    //graph_make_dotfile(g,"foo.dot");
    //graph_del_edge(g,0,1);
    //printf("del 0-1\n");
    //graph_show(g);
    graph_local_complement(g,0);
    //printf("l_0(g)\n"); graph_show(g);
    //local_complement(g,0); printf("l_0^2(g)\n"); graph_show(g);
    printf("%d\n",graph_chromatic_number(g,10)-c);
    //printf("%d\n",graph_clique_number(g)-c);
  }
  return 0;
}

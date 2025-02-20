// KMB 2007 Feb 07
// Article math.PR/0702061 Line-of-sight percolation, Bela Bollobas, Svante Janson, Oliver Riordan Mon 5 Feb 2007
// my model is not exactly the same
// make try_line_of_sight_model
// ./try_line_of_sight_model >| ~/try_planarg/los.his
// ./try_line_of_sight_model | /home/kbriggs/try_planarg/plot_histograms_01.py | p

#include <math.h>
#include "../vn_graph.h"

histogram_t chi,Del,mc;

int graph_los_iterator(unsigned int n, double omega, double p, int f(graph_t,int,void*), void* cd) {
  int i,j,k,l,ig=0;
  double dn=n,x[n*n],y[n*n];
  graph_t g;
  double distance(double xi, double yi, double xj, double yj) {
    if (xi==xj) return fabs(yi-yj);
    if (yi==yj) return fabs(xi-xj);
    return 1e300;
  }
  do { // iterate over graphs
    k=0; // nodes so far
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        if ((double)rand()<p*(double)RAND_MAX) { // pick nodes to be in graph
          x[k]=i/dn;
          y[k]=j/dn;
          k++;
        } 
      }
    }
    g=graph_new(k);
    for (i=0; i<k; i++) {
      for (j=i+1; j<k; j++) {
        if (distance(x[i],y[i],x[j],y[j])<omega) {
          if (1) //degree(g,i)<4 && degree(g,j)<4) // try to keep it planar
            graph_add_edge(g,i,j);
        } 
      }
    }
    l=f(g,ig++,cd);
    graph_clear(g);
  } while (l);
  return l;
}

int f(graph_t g, int ig, void* cd) {
  int c,ng=1000;
  c=graph_max_cluster(g);
  if (c>0 && c<100000)
    histogram_add(mc,c);
  else {
    fprintf(stderr,"mc=%d\n",c);
    exit(1);
  }
  //histogram_add(Del,graph_max_degree(g));
  //c=graph_chromatic_number(g,10);
  //if (c>0) histogram_add(chi,c);
  c=1;
  if (c<=4 && ig<=ng) return 1;
  if (c==5 || ig>ng) graph_make_dotfile_colored(g,"los.dot");
  return 0;
}

int main() {
  int i;
  double p;
  mc=histogram_new("max_cluster(los)");
  chi=histogram_new("chi(los)");
  Del=histogram_new("Delta(los)");
  for (i=1; i<10; i++) {
    histogram_empty(mc);
    histogram_empty(Del);
    histogram_empty(chi);
    p=0.1*i;
    printf("p=%g\n",p);
    graph_los_iterator(32,0.3,p,f,NULL);
    histogram_show(mc);
    //histogram_show(chi);
  }
  //histogram_show(Del);
  return 0;
}

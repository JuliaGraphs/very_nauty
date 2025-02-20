// KMB 2007 Jan 17
// 
// make rlg_poisson_properties && ./rlg_poisson_properties
// /home/kbriggs/python/plot_histograms.py < rlg_poisson_n.his | p -q1 -X"n" -Y"\*x"
// /home/kbriggs/python/plot_histograms.py < rlg_poisson_chi.his | p -q1 -X"\*t" -Y"\*x"
// /home/kbriggs/python/plot_histograms.py < rlg_poisson_eta.his | p -q1 -X"\*t" -Y"\*x'"
// /home/kbriggs/python/plot_histograms.py < rlg_poisson_ome.his | p -q1 -X"\*t" -Y"\*w"

#include <math.h>
#include "vn_graph.h"

typedef struct {
  histogram_t chi;
  histogram_t eta;
  histogram_t ome;
  histogram_t n;
} graph_data_t;

int f(graph_t g, int i, void *cd) {
  histogram_add(((graph_data_t*)cd)->n,nnodes(g)); 
  // timeout will return -1, which histogram_add ignores...
  histogram_add(((graph_data_t*)cd)->eta,graph_edge_chromatic_number(g,10)); 
  histogram_add(((graph_data_t*)cd)->ome,graph_clique_number(g)); 
  histogram_add(((graph_data_t*)cd)->chi,graph_chromatic_number(g,10)); 
  if (i>=999) {
    graph_make_dotfile_colored(g,"rlg_poisson.dot");
    return 0; // done
  }
  //if (histogram_mean_done(h,1.0e-2)) return 0; // done
  return 1; // continue
}

int main(int argc, char* argv[]) {
  double tau,maxtau=7.0;
  graph_data_t cd;
  FILE* fchi=fopen("rlg_poisson_chi.his","w");
  FILE* feta=fopen("rlg_poisson_eta.his","w");
  FILE* fome=fopen("rlg_poisson_ome.his","w");
  FILE* fn=fopen("rlg_poisson_n.his","w");
  cd.chi=histogram_new("chromatic number of random line graph");
  cd.eta=histogram_new("edge chromatic number of random line graph");
  cd.ome=histogram_new("clique number of random line graph");
  cd.n  =histogram_new("size of random line graph");
  if (argc>1) maxtau=atof(argv[1]);
  for (tau=1.0; tau<=maxtau; tau*=1.2) {
    fprintf(fchi,"tau=%g\n",tau);
    fprintf(feta,"tau=%g\n",tau);
    fprintf(fome,"tau=%g\n",tau);
    fprintf(fn,"tau=%g\n",tau);
    histogram_empty(cd.chi);
    histogram_empty(cd.eta);
    histogram_empty(cd.ome);
    histogram_empty(cd.n);
    graph_random_line_graph_poisson_iterator(tau,f,(void*)&cd);
    histogram_write(fchi,cd.chi);
    histogram_write(feta,cd.eta);
    histogram_write(fome,cd.ome);
    histogram_write(fn,cd.n);
  }
  histogram_clear(cd.chi);
  histogram_clear(cd.eta);
  histogram_clear(cd.ome);
  histogram_clear(cd.n);
  fclose(fchi);
  fclose(feta);
  fclose(fome);
  fclose(fn);
  return 0;
}


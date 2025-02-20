// KMB 2007 Jan 12
// 
// make rlg_properties && ./rlg_properties
// /home/kbriggs/python/plot_histograms.py < rlg_chi.his | p -q1 -Xn -Y"\*x"
// /home/kbriggs/python/plot_histograms.py < rlg_eta.his | p -q1 -Xn -Y"\*x'"
// /home/kbriggs/python/plot_histograms.py < rlg_ome.his | p -q1 -Xn -Y"\*w"

#include <math.h>
#include "vn_graph.h"

typedef struct {
  histogram_t chi;
  histogram_t eta;
  histogram_t ome;
} graph_data_t;

int f(graph_t g, int i, void *cd) {
  // timeout will return -1, which histogram_add ignores...
  histogram_add(((graph_data_t*)cd)->chi,graph_chromatic_number(g,10)); 
  histogram_add(((graph_data_t*)cd)->eta,graph_edge_chromatic_number(g,10)); 
  histogram_add(((graph_data_t*)cd)->ome,graph_clique_number(g)); 
  if (i>=999) {
    graph_make_dotfile_colored(g,"rlg.dot");
    return 0; // done
  }
  //if (histogram_mean_done(h,1.0e-2)) return 0; // done
  return 1; // continue
}

int main(int argc, char* argv[]) {
  unsigned int nlines,maxnlines=12;
  graph_data_t cd;
  FILE* fchi=fopen("rlg_chi.his","w");
  FILE* feta=fopen("rlg_eta.his","w");
  FILE* fome=fopen("rlg_ome.his","w");
  cd.chi=histogram_new("chromatic number of random line graph");
  cd.eta=histogram_new("edge chromatic number of random line graph");
  cd.ome=histogram_new("clique number of random line graph");
  if (argc>1) maxnlines=atoi(argv[1]);    // nlines
  for (nlines=2; nlines<=maxnlines; nlines+=1) {
    fprintf(fchi,"nlines=%d\n",nlines);
    fprintf(feta,"nlines=%d\n",nlines);
    fprintf(fome,"nlines=%d\n",nlines);
    histogram_empty(cd.chi);
    histogram_empty(cd.eta);
    histogram_empty(cd.ome);
    graph_random_line_graph_iterator(nlines,f,(void*)&cd);
    histogram_write(fchi,cd.chi);
    histogram_write(feta,cd.eta);
    histogram_write(fome,cd.ome);
  }
  histogram_clear(cd.chi);
  histogram_clear(cd.eta);
  histogram_clear(cd.ome);
  fclose(fchi);
  fclose(feta);
  fclose(fome);
  return 0;
}


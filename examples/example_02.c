// KMB 2006 Jul 17
// test graph_sequential_color against graph_chromatic_number
// with timeout on graph_sequential_color - set it to 0 for no timeout

#include "vn_graph.h"

int main(int argc, char* argv[]) {
  const int verbose=0;
  const int timeout=5;   // timeout in seconds for graph_chromatic_number
  int n=100;              // GRG parameter
  double r=0.2;           // GRG parameter
  int ngraphs=10000;
  int i,chi,printfreq;
  char chi_lbl[200],fn[100];
  histogram_t h;
  if (argc>1) n=atoi(argv[1]);
  if (argc>2) r=atof(argv[2]);
  if (argc>3) ngraphs=atoi(argv[3]);
  printfreq=ngraphs/100;
  graph_t g=graph_new(n);
  printf("example_02.c: GRG(%d,%.3f) chi\n",n,r);
  snprintf(chi_lbl,200,"Chromatic number of GRG(%d,%.3f) chi",n,r);
  h=histogram_new(chi_lbl);
  for (i=0; i<ngraphs; i++) {
    if (verbose || (i+1)%printfreq==0) printf("%dth graph, mean=%g\n",i+1,histogram_mean(h));
    graph_grg_torus(g,r);
    chi=graph_chromatic_number(g,timeout);
    if (chi>0) { // got exact chi
      if (verbose) printf("chi=%2d\n",chi);
      histogram_add(h,chi);
    } else {     // failed to get exact chi within timeout
      printf("graph %d timeout at %d seconds\n",i,timeout);
      if (graph_check_coloring(g)) {
        snprintf(fn,100,"/tmp/grg_n%d_r%.3f_%d.dot",n,r,i);
        fprintf(stderr,"graph_check_coloring=%d\n",graph_check_coloring(g));
        graph_make_dotfile_colored(g,fn);
      }
    }
    if (i>100 && histogram_mean_done(h,1e-2)) goto done;
  }
  fprintf(stderr,"mean not converged after %d graphs!\n",i);
  done:
  graph_clear(g);
  histogram_show(h);
  histogram_clear(h);
  return 0;
}

/* K M Briggs 2007 Jan 17

Pearson correlation coefficient of edge degrees in grg.

make grg_degree_correlation_fast && ./grg_degree_correlation_fast | p -x 0.5 20 -lx -y 0 1

mcopy -o -v grg_degree_correlation_fast.c a:
   
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../vn_graph.h"

#ifndef M_PI
#define M_PI 3.14159263538979
#endif

static const double RAND_MAXp1=((double)RAND_MAX+1.0);
#define UNIFORM_01 (rand()/RAND_MAXp1) /* uniform(0,1) */

double m,s0,s1,s2; // globals for statistics

int f(graph_t g,int k,void* cd) {
  // return after getting data from l edges
  int i,j,di,dj,ne=0,l=*(int*)cd;
  do {
    i=(int)floor(nnodes(g)*UNIFORM_01); // a random node
    if ((di=degree(g,i))) { // node i is not isolated
      j=(int)floor(di*UNIFORM_01); // i--j is a random edge
      dj=degree(g,g->a[i][j]);
      ne++;
      m++;
      s0+=di+dj;
      s1+=di*dj;
      s2+=di*di+dj*dj;
    }
  } while (ne<l);
  if (k>100000) return 0;
  return 1;
}

double edge_prob(double r) {
  if (r<=0.5) return M_PI*r*r;
  if (r*r<=0.5) return 2*sqrt(r*r-0.25)+(M_PI-4*acos(0.5/r))*r*r;
  return 1.0;
}

int main() {
  int n=100; // nnodes
  int l=5;   // edges looked at per graph
  double c,d,e,r=0.25;
  //FILE* fl=fopen("grg_degree_correlation.dat","w");
  d=edge_prob(r)*(n-1);
  do { // r loop
    m=s0=s1=s2=0.0;
    fprintf(stderr,"r=%g mean degree=%g\n",r,d);
    graph_grg_torus_iterator(n,r,f,(void*)&l);
    e=s0/2.0/m;
    e=e*e;
    c=(s1/m-e)/(s2/2.0/m-e); // Pearson correlation coefficient
    printf("%g %g\n",d,c);
    fflush(stdout);
    r*=0.95;
    d=edge_prob(r)*(n-1);
  } while (d>0.5);
  //fclose(fl);
  return 0;
}

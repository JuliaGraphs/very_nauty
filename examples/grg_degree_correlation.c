/* K M Briggs 2006 Dec 14

make grg_degree_correlation
   
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../vn_graph.h"
#include "cov_update.c"

#ifndef M_PI
#define M_PI 3.14159263538979
#endif

static const double RAND_MAXp1=((double)RAND_MAX+1.0);
#define UNIFORM_01 (rand()/RAND_MAXp1) /* uniform[0,1) */

int n=1000;
graph_t g;
FILE* f;

void edge_deg(int l) {
  // return after getting data from at most l edges
  int i,j,di,ne=0;
  do {
    i=(int)floor(nnodes(g)*UNIFORM_01); // a random node
    if ((di=degree(g,i))) {
      j=(int)floor(di*UNIFORM_01); // a random edge
      fprintf(f,"%d %d\n",di,degree(g,j));
      update(di,degree(g,j));
    }
  } while (++ne<l);
}

double edge_prob(double r) {
  if (r<=0.5) return M_PI*r*r;
  if (r*r<=0.5) return 2*sqrt(r*r-0.25)+(M_PI-4*acos(0.5/r))*r*r;
  return 1.0;
}

int main() {
  int j,k;
  double d,z,mc,ssc,c,rd,r=0.4;
  f=fopen("foo","w");
  g=graph_new(n);
  do { // r loop
    r*=0.95;
    d=edge_prob(r)*(n-1);
    fprintf(stderr,"r=%g mean degree=%g ",r,d);
    k=mc=ssc=count=0;
    rd=nan("");
    for (j=0; j<100; j++) { // trials
      graph_grg_torus(g,r);
      edge_deg(n/10);
      if (j>10 && j%10==0) { // check convergence
        c=correlation();
        fprintf(stderr,"c=%g\n",c);
        if (c==c) { // might be NaN
          if (c>0.99) {
            fprintf(stderr,"converged>0.99: reldev(c)=%g\n",rd);
            goto done; // corrrelation done
          }
          k++;
          z=mc;
          mc+=(c-mc)/k;      // mean correlation
          ssc+=(c-mc)*(c-z); // sum of squares
          rd=sqrt(ssc)/k/mc; // relative deviation
          if (k>5 && (rd<1.0e-2 || c>0.999)) {
            fprintf(stderr,"converged: mean corr=%g reldev(c)=%g\n",mc,rd);
            goto done; // corrrelation done
          }
        }
      }
    }
    fprintf(stderr,"not converged: mean correlation=%g reldev(c)=%g\n",mc,rd);
    //printf("standard deviations of x0 = %g\n",stdev(0));
    //printf("standard deviations of x1 = %g\n",stdev(1));
    //printf("covar(x0,x1)= %g\n",covar());
    //printf("correlation coefficient r(x0,x1) = %g\n",correlation());
    //printf("radius=%g meandeg=%g correlation=%g\n",r,3.14159263538979*(n-1)*r*r,correlation());
    done:
    c=correlation();
    if (c==c) printf("%g %g\n",d,c);
    break;
    fflush(stdout);
  } while (d>1e-1);
  graph_clear(g);
  fclose(f);
  return 0;
}

/* K M Briggs 2007 Jan 24

Properties of edge degrees in grg.

make graph_degree_mutual_information && ./graph_degree_mutual_information | p -x 0.5 20 -lx -y 0 1

mcopy -o -v graph_degree_mutual_information.c a:
   
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../vn_graph.h"

static const double RAND_MAXp1=((double)RAND_MAX+1.0);
#define UNIFORM_01 (rand()/RAND_MAXp1) /* uniform(0,1) */

// globals for statistics
double m,s0,s1,s2;
double *px,*py,**pxy;

int f(graph_t g,int k,void* cd) {
  // return after getting data from l edges
  int i,j,di,dj,ne=0,l=*(int*)cd;
  do {
    i=(int)floor(nnodes(g)*UNIFORM_01); // a random node
    di=degree(g,i);
    if (1) {
      do {
        j=(int)floor(nnodes(g)*UNIFORM_01); // i--j is a random edge
        //j=(int)floor(di*UNIFORM_01); // i--j is a random edge
      } while (j==i);
      dj=degree(g,j);
      //dj=degree(g,g->a[i][j]);
      ne++;
      m++;
      px[di]++;
      py[dj]++;
      pxy[di][dj]++;
      //s0+=di+dj;
      //s1+=di*dj;
      //s2+=di*di+dj*dj;
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

double mutual_information(int n) {
  double s=0.0;
  int i,j,k=n;
  do {
    k--;
  } while (k>=1 && (px[k]<4 || py[k]<4 || pxy[k][k]<4));
  fprintf(stderr,"k=%d\n",k);
  for (i=1; i<k; i++) for (j=1; j<k; j++) {
    //fprintf(stderr,"%g %g %g\n",px[i],py[j],pxy[i][j]);
    s+=pxy[i][j]*log(m*pxy[i][j]/px[i]/px[j]);
  }
  return s/m;
}

int main() {
  int n=100; // nnodes
  int l=1;   // edges looked at per graph
  int i,j;
  double c,d,e,mi,r=0.5;
  //FILE* fl=fopen("grg_degree_correlation.dat","w");
  d=edge_prob(r)*(n-1);
  px =malloc(n*sizeof(double));
  py =malloc(n*sizeof(double));
  pxy=malloc(n*sizeof(double*));
  for (i=0; i<n; i++) pxy[i]=malloc(n*sizeof(double));
  do { // r loop
    m=s0=s1=s2=0.0;
    for (i=0; i<n; i++) {
       px[i]=py[i]=0.0;
       for (j=0; j<n; j++) pxy[i][j]=0.0;
    }
    fprintf(stderr,"r=%g mean degree=%g\n",r,d);
    //graph_gnp_iterator(n,r,f,(void*)&l);
    graph_grg_torus_iterator(n,r,f,(void*)&l);
    e=s0/2.0/m;
    e=e*e;
    c=(s1/m-e)/(s2/2.0/m-e); // Pearson correlation coefficient
    mi=mutual_information(n);
    if (mi==mi) printf("%g %g\n",d,mi);
    fflush(stdout);
    r*=0.95;
    d=edge_prob(r)*(n-1);
  } while (d>0.5);
  //fclose(fl);
  return 0;
}

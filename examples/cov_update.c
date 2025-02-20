/* K M Briggs 2006 Oct 24

   Covariance and correlation in 2d

   Taken from pipemath-1.2/covariance.c

   Self-test:
     gcc -std=c99 -Wall -O -DTEXT_COV_UPDATE cov_update.c -lm  && ./a.out
   
   mcopy -o -v cov_update.c a:
   
   FIXME: get rid of globals!
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double nan(const char *tagp);
double *mu,*delta,**zeta,count=0;

int update(double x0, double x1) {
  // add pair (x0,x1) to calculations.
  // returns number of pairs used so far.
  int i,j;
  static int m=0;
  size_t n=2,szdb=sizeof(double);
  double f,g;
  if (count==0) { // reset - new data set
    mu   =malloc(n*szdb);
    delta=malloc(n*szdb);
    zeta =malloc(n*sizeof(double*));
    for (i=0; i<n; i++) {
      zeta[i]=malloc((i+1)*szdb);
      for (j=0; j<=i; j++) zeta[i][j]=0.0;
      mu[i]=0.0;
    }
    m=0;
  } 
  m++;
  count++;
  f=(m-1.0)/m;
  // i=0...
  delta[0]=x0-mu[0];
  mu[0]+=delta[0]/m;
  g=f*delta[0];
  for (j=0; j<=0; j++) zeta[0][j]+=g*delta[j];
  // i=1...
  delta[1]=x1-mu[1];
  mu[1]+=delta[1]/m;
  g=f*delta[1];
  for (j=0; j<=1; j++) zeta[1][j]+=g*delta[j];
  return m;
} 

double stdev(int i) { // stdev of ith component (i=0,1) after m values 
  if (count==0) return nan("");
  return sqrt(zeta[i][i]/(count-1));
}

double covar() { // covariance after m values 
  if (count==0) return nan("");
  return zeta[1][0]/(count-1);
}

double correlation() { // correlation after m values 
  if (count==0) return 0.0;
  return zeta[1][0]/sqrt(zeta[1][1]*zeta[0][0]);
}

#ifdef TEXT_COV_UPDATE
int main() { // self-test program
  int i;
  double x;
  for (i=0; i<1000; i++) { // test data
    x=450+0.001*i;
    update(x,x+0.5*(rand()/(double)RAND_MAX-0.5));
  }
  printf("standard deviations of x0 = %g\n",stdev(0));
  printf("standard deviations of x1 = %g\n",stdev(1));
  printf("covar(x0,x1)= %g\n",covar());
  printf("correlation coefficient r(x0,x1) = %g\n",correlation());
  return 0;
}
#endif

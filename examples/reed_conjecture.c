// KMB 2006 Jul 25
// B. Reed 'omega, Delta, and chi' J Graph Theory 27, 177-212 (1998)
// Can every graph be colored with at most ceil(Delta+1+omega)/2 colors?
// http://www.imada.sdu.dk/Research/Graphcol/new_problems.html
// conj 1: let gap1(g)=ceil((Delta(g)+1+omega(g))/2)-chi(g) and look for gap1<0
// conj 2: let gap2(g)=(2(Delta(g)+1)+omega(g))/3-chi(g) and look for gap2<0

/*
geng -s 11 | ./reed_conjecture
for i in 1 2 3 4 5 6 7 8 9 10 11 12; do geng -t $i | ./reed_conjecture; done
for i in 13; do geng -t -d3 -D5 $i | ./reed_conjecture; done # => 372

Below, tf1 means triangle-free cases of gap1=0 (n=14 467871369 total graphs)
Below, tf2 means triangle-free cases of gap2=0
(All with with Delta<=n-3)
n:     1 2 3 4 5 6 7 8 9 10 11 12  13    14   15 
tf1:   0 0 0 0 0 0 0 0 0  0  1 21 826 39889
tf2:   0 0 0 0 0 0 0 0 0  0  1 21 826 39889

All graphs with Delta<=n-4:
n:     1 2 3 4 5 6 7 8 9 10 11 12  13    14   15 
       0 0 0 0 0 0 0 0 0  1 47
*/

/*
http://epubs.siam.org/SIDMA/volume-20/art_62654.html
SIDMA, Volume 20 Issue 2
SIAM Journal on Discrete Mathematics
pages 529–535, © 2006 Society for Industrial and Applied Mathematics
On Graph Associations
Landon Rabern

Abstract.

We introduce a notion of vertex association and consider sequences of these associations. This allows for slick proofs of a few known theorems as well as showing that for any induced subgraph H of G, chi(G) \leq chi(H) + 1/2 (omega (G) + |G| - |H| - 1). As a special case of this, we have chi(G) \leq \lceil omega(G) + tau(G) / 2 \rceil (here chi(G) denotes the chromatic number, omega(G) the clique number, and tau(G) the vertex cover number), which is a generalization of the Nordhaus--Gaddum upper bound. In addition, this settles a conjecture of Reed that chi(G) \leq \lceil omega(G) + Delta(G) + 1 / 2 \rceil in the case when delta(\overline G) \leq omega(\overline G).
*/



#include <stdlib.h>
#include <math.h>
#include "vn_graph.h"

const int v=0; // verbosity

histogram_t geng_reader(FILE* f, int tf) {
  char *s;
  unsigned long i=0;
  int n=0,Delta,omega=0,chi,gap1;
  double gap2;
  histogram_t gap1s=histogram_new("gap1");
  histogram_t gap2s=histogram_new("gap2");
  graph_t g;
  char fn[50];
  while ((s=geng_getline(f))) {
    i++;
    n=geng_graphsize(s);
    if (i%10000000==0) {
      printf("not finished: n=%d i=%lu s=%s",n,i,s);
      histogram_show(gap1s);
      histogram_show(gap2s);
    }
    g=geng_stringtograph(s);
    //printf("i=%d s=%s",i,s);
    Delta=graph_max_degree(g);
    if (tf) omega=2; // triangle-free
    else    omega=graph_clique_number(g);
    if (omega<Delta-1 && Delta<n-1) // exclude cases known to be true
     if (Delta<n-3) { // required by Rabern result SIDMA 2006
      chi=graph_chromatic_number(g,1000);
      if (chi==0) {
        printf("timeout at g=%s\n",s);
      } else {
        gap1=(int)ceil((Delta+1.0+omega)/2.0)-chi;
        if (gap1==0) {
          if (v) printf("gap1=0 Delta=%d omega=%d chi=%d s=%s",Delta,omega,chi,s);
	  if (v && n<13) {
            graph_show(g);
            snprintf(fn,50,"zero_gap1_%lu.dot",i);
            graph_make_dotfile_colored(g,fn);
	  }
        }
        if (gap1<0) {
          printf("counterexample to conjecture 1: n=%d Delta=%d omega=%d chi=%d g=%s\n",n,Delta,omega,chi,s);
          snprintf(fn,50,"counterexample1_%lu.dot",i);
          graph_make_dotfile_colored(g,fn);
          graph_show(g);
        } else {
          histogram_add(gap1s,gap1);
	}
        gap2=(2.0*(Delta+1.0)+omega)/3.0-chi;
        if (Delta>=3 && gap2==0) {
          if (v) printf("gap2=0 Delta=%d omega=%d chi=%d s=%s",Delta,omega,chi,s);
	  if (v && n<13) {
            graph_show(g);
            snprintf(fn,50,"zero_gap2_%lu.dot",i);
            graph_make_dotfile_colored(g,fn);
	  }
        }
        if (Delta>=3 && gap2<0) {
          printf("counterexample to conjecture 2: n=%d Delta=%d omega=%d chi=%d g=%s\n",n,Delta,omega,chi,s);
          snprintf(fn,50,"counterexample2_%lu.dot",i);
          graph_make_dotfile_colored(g,fn);
          graph_show(g);
        } else {
          histogram_add(gap2s,gap2);
	}
      }
    }
    if (s) free(s);
    graph_clear(g);
  }
  printf("n=%d done\n",n);
  histogram_show(gap1s);
  histogram_show(gap2s);
  return gap1s;
}

int main(int argc, char* argv[]) {
  int tf=0; // triangle-free
  if (argc>1) tf=1;
  if (tf) printf("Assuming triangle-free\n");
  //histogram_t gap1s=geng_reader(stdin);
  geng_reader(stdin,tf);
  return 0;
}

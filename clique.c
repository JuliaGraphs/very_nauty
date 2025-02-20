/* Modified by K M Briggs: dynamic memory allocation and various tidy-ups */
/* 

http://locus.siam.org/SICOMP/volume-06/art_0206036.html

S. Tsukiyama, M. Ide, H. Ariyoshi, and I. Shirawaka. A new algorithm for generating all the maximal independent sets. SIAM J. Computing, 6:505--517, 1977.

Abstract.  The problem of generating all the maximal independent sets (or maximal cliques) of a given graph is fundamental in graph theory and is also one of the most important in terms of the application of graph theory. In this paper, we present a new efficient algorithm for generating all the maximal independent sets, for which processing time and memory space are bounded by $O(nm\mu)$ and $O(n+m)$, respectively, where $n$, $m$, and $\mu$ are the numbers of vertices, edges, and maximal independent sets of a graph.

Keywords. algorithm, backtracking, graph, maximal clique, maximal independent set

 */

struct {
  int maxsz;
  int verbose;
  int matrix_size;
  int mode;
  int *deg;
  int *IS;
  int **adj_list;
  int **bucket;
} clqdata;

/* MaxClique - for finding maximal cliques and independent sets of graphs
 * Written by Kevin O'Neill
 * Based on algorithm of Tsukiyama, Ide, Ariyoshi, and Shirakawa
 * Latest version: December 18, 2003
 *
 * Copyright (C) 2003  Kevin O'Neill
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 * The author can be contacted at oneill@cs.cornell.edu, or by mail at
 *   Computer Science Department
 *   4130 Upson Hall
 *   Cornell University
 *   Ithaca, NY 14853-7501
 */


static void backtrack(int i);
static void output_IS(int* is);

static void backtrack(int i) {
  int c,x,f,j,y,k,z,*ax,*ay;
  if (i>=(clqdata.matrix_size-1)) {
    output_IS(clqdata.IS);
  } else {
    x=i+1;
    ax=clqdata.adj_list[x];
    c=0;
    for (j=0; j<clqdata.deg[x] && ((y=ax[j])<=i); j++) if (clqdata.IS[y]==0) c++;
    if (c==0) {
      for (j=0; j<clqdata.deg[x] && ((y=ax[j])<=i); j++) clqdata.IS[y]++;
      backtrack(x);
      for (j=0; j<clqdata.deg[x] && ((y=ax[j])<=i); j++) clqdata.IS[y]--;
    } else {
      clqdata.IS[x]=c;
      backtrack(x);
      clqdata.IS[x]=0;
      f=1;
      for (j=0; j<clqdata.deg[x] && ((y=ax[j])<=i); j++) {
        if (clqdata.IS[y]==0) {
          clqdata.bucket[x][j]=1;
          ay=clqdata.adj_list[y];
          for (k=0; k<clqdata.deg[y] && ((z=ay[k])<=i); k++)
            if (--clqdata.IS[z]==0) f=0;
        }
        clqdata.IS[y]++;
      }
      if (f) backtrack(x);
      for (j=0; j<clqdata.deg[x] && ((y=ax[j])<=i); j++) clqdata.IS[y]--;
      for (j=0; j<clqdata.deg[x]; j++) if (clqdata.bucket[x][j]) {
        y=ax[j];
        ay=clqdata.adj_list[y];
        for (k=0; k<clqdata.deg[y] && ((z=ay[k])<=i); k++) clqdata.IS[z]++;
        clqdata.bucket[x][j]=0;
      }
    }
  }
}

static void output_IS(int* is) {
  int ISIndex,sz=0;
  if (clqdata.verbose) {
    if (clqdata.mode==0) fprintf(stdout, "Maximal clique: ");
    else fprintf(stdout, "Maximal independent set: ");
  }
  for (ISIndex=0; ISIndex<clqdata.matrix_size; ISIndex++) {
    if (is[ISIndex]==0) { // IS[v]=0 iff v is in MIS
      sz++;
      if (clqdata.verbose) fprintf(stdout, "%d ", ISIndex);
    }
  }
  if (clqdata.verbose) { fprintf(stdout, "\n"); fflush(stdout); }
  clqdata.maxsz=sz>clqdata.maxsz?sz:clqdata.maxsz;
}

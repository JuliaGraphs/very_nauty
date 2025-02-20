/* Modified by Keith M Briggs:
     dynamic memory allocation
     global struct cg
     timeout 
     various tidy-ups
     Code is probably sufficiently debugged, and may be useful to some people.
*/

struct { // coloring globals
  graph_t g;
  struct tms buffer;		/* structure for timing  */
  int tps; 
  clock_t current_time, start_time, timeout;
  double utime;
  int verbose;
  int BestColoring;
  int num_node;
  int prob_count;
  int lb;
  int num_prob, max_prob;
  int best_clique;
  int *ColorClass;
  int *Order;
  int *Handled;
  int *ColorCount;
  int **ColorAdj;
  int **adj;
} cg;

/*
 * The author of this software is Michael Trick.  Copyright (c) 1994 by
 * Michael Trick.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR DOES NOT MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/*
   COLOR.C: Easy code for graph coloring
   Author: Michael A. Trick, Carnegie Mellon University, trick+@cmu.edu
   Last Modified: November 2, 1994


Code is probably insufficiently debugged, but may be useful to some people.

For more information on this code, see Anuj Mehrotra and Michael A. Trick,
"A column generation approach to graph coloring", GSIA Technical report series.

*/


static int greedy_clique(int *valid, int *clique) {
  int i, j, k;
  int max;
  int place, done;
  int order[cg.num_node+1]; // allocation on stack
  int weight[cg.num_node];
  bzero(order,(cg.num_node+1)*sizeof(int));
  bzero(weight,cg.num_node*sizeof(int));
  bzero(clique,cg.num_node*sizeof(int));
  place=0;
  for (i=0; i<cg.num_node; i++) if (valid[i]) order[place++]=i;
  for (i=0; i<cg.num_node; i++)
    if (valid[i])
      for (j=0; j<cg.num_node; j++) if (valid[j] && cg.adj[i][j]) weight[i]++;
  done = 0;
  while (!done) {
    done = 1;
    for (i = 0; i < place - 1; i++) {
      j = order[i];
      k = order[i + 1];
      if (weight[j] < weight[k]) {
	order[i] = k;
	order[i + 1] = j;
	done = 0;
      }
    }
  }
  clique[order[0]] = 1;
  for (i=1; i<place; i++) {
    j=order[i];
    for (k=0; k<i; k++) if (clique[order[k]] && !cg.adj[j][order[k]]) break;
    clique[j]=(k==i);
  }
  max=0;
  for (i=0; i<place; i++) if (clique[order[i]]) max++;
  return max;
}

/* Target is a goal value:  once a clique is found with value target
   it is possible to return

   Lower is a bound representing an already found clique:  once it is
   determined that no clique exists with value better than lower, it
   is permitted to return with a suboptimal clique.

   Note, to find a clique of value 1, it is not permitted to just set
   the lower to 1:  the recursion will not work.  Lower represents a
   value that is the goal for the recursion.
*/

static int max_w_clique(int valid[], int clique[], int lower, int target) {
  int start,i,j,k,incumb,new_weight,place,finish,done,place1,total_left;
  int *valid1,*clique1,*order,*value;
  cg.num_prob++;
  if (cg.num_prob>cg.max_prob) return -1;
  for (j=0; j<cg.num_node; j++) clique[j] = 0;
  total_left=0;
  for (i=0; i<cg.num_node; i++) if (valid[i]) total_left++;
  if (total_left<lower) return 0;
  value=(int*)Calloc(cg.num_node,sizeof(int));
  incumb=greedy_clique(valid, clique);
  if (incumb>=target) {
    Free(value,cg.num_node*sizeof(int));
    return incumb;
  }
  if (incumb>cg.best_clique) cg.best_clique=incumb;
  place=0;
  order=Calloc(cg.num_node+1,sizeof(int));
  if (order==NULL) {
    fprintf(stderr,"Calloc failure at trick.c line 131, num_node=%d\n",cg.num_node);
    exit(1);
  }
  for (i = 0; i < cg.num_node; i++) {
    if (clique[i]) {
      order[place] = i;
      total_left--;
      place++;
    }
  }
  start = place;
  for (i = 0; i < cg.num_node; i++) {
    if (!clique[i] && valid[i]) {
      order[place] = i;
      place++;
    }
  }
  finish = place;
  for (place = start; place < finish; place++) {
    i = order[place];
    value[i] = 0;
    for (j = 0; j < cg.num_node; j++) if (valid[j] && cg.adj[i][j]) value[i]++;
  }
  done = 0;
  while (!done) {
    done = 1;
    for (place = start; place < finish - 1; place++) {
      i = order[place];
      j = order[place + 1];
      if (value[i] < value[j]) {
	order[place] = j;
	order[place + 1] = i;
	done = 0;
      }
    }
  }
  Free(value,cg.num_node*sizeof(int));
  valid1=(int*)Calloc(cg.num_node,sizeof(int));
  clique1=(int*)Calloc(cg.num_node,sizeof(int));
  for (place = start; place < finish; place++) {
    if (incumb + total_left < lower) { 
      incumb=0;
      break;
    }
    j = order[place];
    total_left--;
    if (clique[j]) continue;
    for (place1 = 0; place1 < place; place1++) {
      k = order[place1];
      valid1[k]=valid[k] && cg.adj[j][k];
    }
    new_weight = max_w_clique(valid1, clique1, incumb - 1, target - 1);
    if (new_weight+1>incumb) {
      incumb=new_weight+1;
      for (k=0; k<cg.num_node; k++) clique[k]=clique1[k];
      clique[j]=1;
      if (incumb>cg.best_clique) cg.best_clique=incumb;
    }
    if (incumb >= target) break;
  }
  Free(valid1,cg.num_node*sizeof(int));
  Free(clique1,cg.num_node*sizeof(int));
  Free(order,(cg.num_node+1)*sizeof(int));
  return incumb;
}

static void AssignColor(int node, int color) {
  int node1;
  cg.ColorClass[node]=color;
  for (node1=0; node1<node; node1++) if (cg.adj[node][node1]) {
    if (!cg.ColorAdj[node1][color]) cg.ColorCount[node1]++;
    cg.ColorAdj[node1][color]++;
    cg.ColorAdj[node1][0]--;
  }
  for (node1=node+1; node1<cg.num_node; node1++) if (cg.adj[node][node1]) {
    if (!cg.ColorAdj[node1][color]) cg.ColorCount[node1]++;
    cg.ColorAdj[node1][color]++;
    cg.ColorAdj[node1][0]--;
  }
}

static void RemoveColor(int node, int color) {
  int node1;
  cg.ColorClass[node]=0;
  for (node1=0; node1<node; node1++) {
    if (cg.adj[node][node1]) {
      cg.ColorAdj[node1][color]--;
      if (cg.ColorAdj[node1][color]==0) cg.ColorCount[node1]--;
      cg.ColorAdj[node1][0]++;
    }
  }
  for (node1=node+1; node1<cg.num_node; node1++) {
    if (cg.adj[node][node1]) {
      cg.ColorAdj[node1][color]--;
      if (cg.ColorAdj[node1][color]==0) cg.ColorCount[node1]--;
      cg.ColorAdj[node1][0]++;
    }
  }
}

static int Color(int i, int current_color) {
  int l,j,new_val;
  int k, max, place;
  cg.prob_count++;
  if (current_color>=cg.BestColoring) return current_color;
  if (cg.BestColoring<=cg.lb) return cg.BestColoring;
  if (i>=cg.num_node) return current_color;
  if (cg.timeout) {  // a timeout is specified
    times(&cg.buffer);
    cg.current_time=cg.buffer.tms_utime;
    if ((cg.current_time-cg.start_time)>cg.tps*cg.timeout) return -1;
  }
  /* Find node with maximum color_adj */
  max=-1;
  place=-1;
  for (k=0; k<cg.num_node; k++) if (!cg.Handled[k]) {
    if ((cg.ColorCount[k]>max) || ((cg.ColorCount[k]==max) && (cg.ColorAdj[k][0]>cg.ColorAdj[place][0]))) {
      max=cg.ColorCount[k];
      place=k;
    }
  }
  if (place==-1) {  // FIXME doesn't work? - did he mean "Graph is empty"?
    printf("Graph is disconnected.  This code needs to be updated for that case.\n");
    exit(1);
  }
  cg.Order[i]=place;
  cg.Handled[place]=1;
  for (j=1; j<=current_color; j++) if (!cg.ColorAdj[place][j]) {
    cg.ColorClass[place]=j;
    AssignColor(place,j);
    new_val=Color(i+1,current_color);
    if (new_val<cg.BestColoring) {
      cg.BestColoring=new_val;
      for (l=0; l<cg.num_node; l++) cg.g->c[l]=cg.ColorClass[l];
    }
    RemoveColor(place,j);
    if (cg.BestColoring<=current_color) {
      cg.Handled[place]=0;
      return cg.BestColoring;
    }
  }
  if (current_color+1<cg.BestColoring) {
    cg.ColorClass[place]=current_color+1;
    AssignColor(place,current_color+1);
    new_val=Color(i+1,current_color+1);
    if (new_val<cg.BestColoring) {
      cg.BestColoring=new_val;
      for (l=0; l<cg.num_node; l++) cg.g->c[l]=cg.ColorClass[l];
    }
    RemoveColor(place,current_color+1);
  }
  cg.Handled[place]=0;
  return cg.BestColoring;
}

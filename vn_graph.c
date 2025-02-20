// KMB 2006 Dec 05, 13
//   - better sizes for graph_t components
//   - non-recursive dfs
//   - new internals with g->b and g->d
//   - neighbour lists are no longer -1-terminated
//   - neighbour lists grow as required - can now do bigger graphs
//   - g->d[i] gives degree of node i
// KMB 2006 Jul 19; Nov 22
// TODO tidy up debug prints

#include <limits.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>

#include "vn_graph.h"

static const double RAND_MAXp1=((double)RAND_MAX+1.0);
#define UNIFORM_01 (rand()/RAND_MAXp1) /* uniform[0,1) */

#undef VN_COUNT_RNG_CALLS
#ifdef VN_COUNT_RNG_CALLS
long long vn_rng_calls=0; 
#define VN_GRAPH_UNIFORM_01 (vn_rng_calls++,rand()/((double)RAND_MAX)) /* uniform(0,1) */
#else
#define VN_GRAPH_UNIFORM_01 (rand()/((double)RAND_MAX)) /* uniform(0,1) */
#endif

//#define MALLOC_FREE_CHECK
#ifdef MALLOC_FREE_CHECK
  static size_t alloc=0;
  #define Malloc(sz) malloc((alloc+=(sz),sz))
  #define Calloc(n,sz) calloc(n,(alloc+=((n)*(sz)),sz))
  #define Realloc(p,oldsz,newsz) realloc(p,(alloc+=((newsz)-(oldsz)),(newsz)))
  #define Free(p,sz) free((alloc-=(sz),p))
  #define print_alloc(msg) fprintf(stderr,"%s: alloc=%ld\n",msg,(long)alloc);
  #define check_alloc(msg) if (alloc) { fprintf(stderr,"%s: alloc=%ld\n",msg,(long)alloc); exit(1); }
#else
  #define Malloc(sz) malloc(sz)
  #define Calloc(n,sz) calloc(n,sz)
  #define Realloc(p,oldsz,newsz) realloc(p,(newsz))
  #define Free(p,sz) free(p)
#endif

inline static void graph_add_edge1(graph_t g, double k);

#include "trick.c"
#include "clique.c"
#include "rainbow.c"
#include "geng_reader.c"

const char* graph_get_version(void) {
  // FIXME use VERY_NAUTY_VERSION
  const char* graph_version="very_nauty-1.1 compiled" " "__DATE__" "__TIME__;
  return graph_version;
}

static size_t adj_list_initial_size=8; // initial size of adjacency lists

void set_adj_list_initial_size(size_t sz) {
  adj_list_initial_size=sz;
}

static void shuffle(int p[], int n, int m);

//! make a new empty graph with n nodes
graph_t graph_new(unsigned int n) {
  size_t sz=adj_list_initial_size;
  int i;
  graph_t g=Malloc(sizeof(struct _graph));
  g->v=Malloc(n*sizeof(char));
  g->c=Malloc(n*sizeof(int));
  g->l=Malloc(n*sizeof(int));
  g->b=Malloc(n*sizeof(size_t));
  g->d=Malloc(n*sizeof(node_t));
  g->a=Malloc(n*sizeof(node_t*));
  if (n<=sz) sz=n; // initial allocation of neighbour lists
  for (i=0; i<n; i++) {
    g->c[i]=-1;
    g->l[i]=-1;
    g->v[i]=0;
    g->d[i]=0;
    g->b[i]=sz;
    g->a[i]=Malloc(sz*sizeof(node_t));
  }
  g->nnodes=n;
  g->nedges=0;
  return g;
}

//! destroy a graph
void graph_clear(graph_t g) {
  int i,n=g->nnodes;
  for (i=0; i<n; i++) Free(g->a[i],n*sizeof(node_t*));
  Free(g->v,n*sizeof(char));
  Free(g->c,n*sizeof(int));
  Free(g->l,n*sizeof(int));
  Free(g->b,n*sizeof(size_t));
  Free(g->d,n*sizeof(node_t));
  Free(g->a,n*sizeof(node_t*));
  Free(g,sizeof(struct _graph));
}

//! add a node to a graph
void graph_add_node(graph_t g) {
  size_t sz=adj_list_initial_size;
  unsigned int i,n=g->nnodes,np1=n+1;
  if (np1<=sz) sz=np1; // initial allocation of neighbour lists
  g->a=Realloc(g->a,n*sizeof(node_t*),np1*sizeof(node_t*));
  g->b=Realloc(g->b,n*sizeof(size_t),np1*sizeof(size_t));
  g->b[n]=sz;
  g->a[n]=Malloc(sz*sizeof(node_t));
  g->v=Realloc(g->v,n*sizeof(char),np1*sizeof(char));
  g->c=Realloc(g->c,n*sizeof(int),np1*sizeof(int));
  g->l=Realloc(g->l,n*sizeof(int),np1*sizeof(int));
  g->d=Realloc(g->d,n*sizeof(node_t),np1*sizeof(node_t));
  g->d[n]=0;
  for (i=0; i<n; i++) {
    g->v[i]=0;
    g->c[i]=-1;
    g->l[i]=-1;
  }
  g->nnodes=np1;
}

//! delete edge u--v from the graph (if it is present)
// return 1 if edge was found and deleted, else 0
int graph_del_edge(graph_t g, node_t u, node_t v) {
  int i,j;
  for (i=0; i<g->d[u]; i++) if (g->a[u][i]==v) {
    memmove(&g->a[u][i],&g->a[u][i+1],(g->d[u]-1-i)*sizeof(node_t));
    g->d[u]--;
    for (j=0; j<g->d[v]; j++) if (g->a[v][j]==u) {
      memmove(&g->a[v][j],&g->a[v][j+1],(g->d[v]-1-j)*sizeof(node_t));
      g->d[v]--;
      break;
    }
    g->nedges--;
    return 1;
  }
  return 0;
}

//! 1 if edge i--j is present in g
int graph_has_edge(graph_t g, node_t i, node_t j) {
  int k;
  for (k=0; k<g->d[i]; k++) if (g->a[i][k]==j) return 1;
  return 0;
}

//! empty a graph
void graph_empty(graph_t g) {
  int i;
  bzero(g->d,g->nnodes*sizeof(node_t)); // degree 0
  bzero(g->v,g->nnodes*sizeof(char)); // not visited
  for (i=0; i<g->nnodes; i++) {
    g->c[i]=-1; // not colored
    g->l[i]=-1; // no cluster
  }
  g->nedges=0;
}

//! a complete graph
void graph_complete(graph_t g) {
  int i,j,n=g->nnodes;
  for (i=0; i<n; i++) {
    g->a[i]=Realloc(g->a[i],g->b[i]*sizeof(node_t),n*sizeof(node_t));
    for (j=0; j<n; j++) if (j!=i) g->a[i][j]=j;
    g->d[i]=n-1;
    g->b[i]=n;
    g->c[i]=-1;
    g->l[i]=-1;
    g->v[i]=0;
  }
  g->nedges=n*(n-1)/2;
}

//! add undirected edge u--v to the graph, if it is not already present
void graph_add_edge(graph_t g, node_t u, node_t v) {
  int i;
  for (i=0; i<g->d[u]; i++) if (g->a[u][i]==v) return;
  assert(g->d[u]<=g->b[u]);
  if (g->d[u]==g->b[u]) {
    g->b[u]+=16;
    g->a[u]=Realloc(g->a[u],(g->b[u]-16)*sizeof(node_t),g->b[u]*sizeof(node_t));
  }
  g->a[u][g->d[u]++]=v;
  assert(g->d[v]<=g->b[v]);
  if (g->d[v]==g->b[v]) {
    g->b[v]+=16;
    g->a[v]=Realloc(g->a[v],(g->b[v]-16)*sizeof(node_t),g->b[v]*sizeof(node_t));
  }
  g->a[v][g->d[v]++]=u;
  g->nedges++;
}

//! add undirected edge u--v to the graph, whether present already or not
// warning - can make multigraphs if abused
void graph_append_edge(graph_t g, node_t u, node_t v) {
  assert(g->d[u]<=g->b[u]);
  if (g->d[u]==g->b[u]) {
    g->b[u]+=16;
    g->a[u]=Realloc(g->a[u],(g->b[u]-16)*sizeof(node_t),g->b[u]*sizeof(node_t));
  }
  g->a[u][g->d[u]++]=v;
  assert(g->d[v]<=g->b[v]);
  if (g->d[v]==g->b[v]) {
    g->b[v]+=16;
    g->a[v]=Realloc(g->a[v],(g->b[v]-16)*sizeof(node_t),g->b[v]*sizeof(node_t));
  }
  g->a[v][g->d[v]++]=u;
  g->nedges++;
}

//! degree of a node
int graph_node_degree(graph_t g, node_t i) {
  return g->d[i];
}

//! minimum degree of a node in g
int graph_min_degree(graph_t g) {
  int i,d,k=INT_MAX;
  for (i=0; i<g->nnodes; i++) if ((d=graph_node_degree(g,i))<k) k=d;
  return k;
}

//! maximum degree of a node in g
int graph_max_degree(graph_t g) {
  int i,d,k=-1;
  for (i=0; i<g->nnodes; i++) if ((d=graph_node_degree(g,i))>k) k=d;
  return k;
}

//! mean degree
double graph_mean_degree(graph_t g) {
  return (2.0*g->nedges)/g->nnodes;
}

//! simple listing of graph
void graph_show(graph_t g) {
  int i,j;
  for (i=0; i<g->nnodes; i++) if (g->d[i]>0) { // i has neighbours
    printf("%d: ",i);
    for (j=0; j<g->d[i]; j++) printf("%d ",g->a[i][j]);
    printf("\n");
  }
}

void graph_show_edge_list(graph_t g, int offset) {
  int i,j;
  printf("%d %d\n",g->nnodes,g->nedges);
  for (i=0; i<g->nnodes; i++) if (g->d[i]>0) { // i has neighbours
    for (j=0; j<g->nnodes; j++) if (g->d[i]>0 && g->a[i][j]>i) printf("%d %d\n",i+offset,g->a[i][j]+offset);
  }
}

void graph_dump_to_file(graph_t g, char* fn) {
  // my edge list format
  int i,j;
  FILE* f=fopen(fn,"w");
  fprintf(f,"n=%d m=%d.\n",g->nnodes,g->nedges);
  for (i=0; i<g->nnodes; i++) if (g->d[i]>0) { // i has neighbours
    for (j=0; j<g->d[i]; j++) if (j>i) fprintf(f,"%d %d\n",i,j);
  }
  fclose(f);
  printf("%s written\n",fn);
}

void graph_to_showg(graph_t g, int n, char* fn) {
  // "showg -e -l0" format
  int i,j;
  FILE* f=fopen(fn,"w");
  fprintf(f,"Graph %d, order %d.\n",n,g->nnodes);
  fprintf(f,"%d %d\n",n,g->nedges);
  for (i=0; i<g->nnodes; i++) if (g->d[i]>0) { // i has neighbours
    for (j=0; j<g->d[i]; j++) if (j>i) fprintf(f,"%d %d  ",i,j);
  }
  fprintf(f,"\n");
  fclose(f);
  printf("%s written:\n",fn);
}

void graph_to_dimacs(graph_t g, char* fn) {
  // ftp://dimacs.rutgers.edu/pub/challenge/graph/doc/ccformat.tex
  // 1-indexed
  int i,j;
  FILE* f=fopen(fn,"w");
  fprintf(f,"c written by graph_to_dimacs\n");
  fprintf(f,"p edge %d %d\n",g->nnodes,g->nedges);
  for (i=0; i<g->nnodes; i++) if (g->d[i]>0) { // i has neighbours
    for (j=0; j<g->d[i]; j++) if (j>i) fprintf(f,"e %d %d\n",i+1,j+1);
  }
  fclose(f);
  printf("%s written\n",fn);
}

void graph_to_theta(graph_t g, char* fn) { // FIXME
  // 1-indexed, complemented, for DSDP Lovasz theta
  int i,j,n=g->nnodes;
  FILE* f=fopen(fn,"w");
  int u[g->nnodes];
  fprintf(f,"%d %d\n",n,n*(n-1)/2-g->nedges);
  for (i=0; i<n; i++) if (g->d[i]>0) { // i has neighbours
    bzero(u,n*sizeof(int));
    for (j=0; j<g->d[i]; j++) u[g->a[i][j]]=1;
    for (j=i+1; j<n; j++) if (!u[j]) fprintf(f,"%d %d\n",i+1,j+1);
  }
  fclose(f);
  printf("%s written: now do \"theta %s\"\n",fn,fn);
}

graph_t graph_line_graph(graph_t g) {
  int i,j,k,m=0,n=g->nedges;
  int e0[n],e1[n];
  graph_t l=graph_new(n);
  for (i=0; i<g->nnodes; i++) {
    for (k=0; k<g->d[i]; k++) { // loop over edges of g
      j=g->a[i][k];
      if (i<j) { // i--j is an edge; do each edge once only
        assert(m<n);
        if (m>n-1) { fprintf(stderr,"Arghhh!!"); exit(1); }
        e0[m]=i;
        e1[m]=j;
        m++;
      }
    }
  }
  assert(m==n);
  for (i=0; i<n; i++) { // loop over possible edges of line graph
    for (j=0; j<n; j++) if (i!=j) {
      if (e0[i]==e0[j] ||
          e1[j]==e1[i] ||
          e1[i]==e0[j] ||
          e0[i]==e1[j]) { // common end-point
            graph_add_edge(l,i,j);
       }
    }
  }
  return l;
}

void graph_make_dotfile(graph_t g, char* fn) {
  int i,j,k;
  FILE* f=fopen(fn,"w");
  fprintf(f,"graph %s {\n","graph_make_dotfile");
  fprintf(f,"  center=true; unflatten=true;\n ratio=1;\n /* mclimit=5;\n nslimit=500;\n nslimit1=500;\n remincross=true; */\n size=\"8,8\";\n node  [shape=\"circle\",size=\"0.1\"]\n node [color=green,fontcolor=red,fontsize=8]\n edge [color=blue];\n graph [concentrate=\"false\"];\n");
  for (i=0; i<g->nnodes; i++) if (g->d[i]>0) { // i has neighbours
    for (k=0; k<g->d[i]; k++) if ((j=g->a[i][k])>i) fprintf(f,"  %d--%d;\n",i,j);
  }
  fprintf(f,"}\n");
  fclose(f);
  printf("%s written: now do \"neato -Tps %s | gv -\"\n",fn,fn);
}

void rgb2hexstring(double r, double g, double b, char *s) {
  double gi=1.0/1.0; // inverse gamma
  r=floor(255.9999999*pow(r,gi));
  g=floor(255.9999999*pow(g,gi));
  b=floor(255.9999999*pow(b,gi));
  snprintf(s,7,"%02x%02x%02x",(int)r,(int)g,(int)b);
}

static double logpoisson(int k, double l) {
  // log Poisson
  // log Pr[x=k]
  return -l+k*log(l)-lgamma(k+1);
}

// http://www.graphviz.org/doc/info/colors.html
void graph_make_dotfile_colored(graph_t g, char* fn) {
  int i,j,k,m=-1;
  char s[7];
  double rd,gr,bl,c;
  FILE* f=fopen(fn,"w");
  fprintf(f,"graph %s {\n","graph_make_dotfile_colored");
  for (i=0; i<g->nnodes; i++) if (g->c[i]>m) m=g->c[i];
  fprintf(f,"  center=true; unflatten=true;\n ratio=1;\n /* mclimit=5;\n nslimit=500;\n nslimit1=500;\n remincross=true; */\n size=\"8,10\";\n node [shape=\"circle\",size=\"0.03\"]\n node [fontsize=10,style=filled]\n edge [style=\"setlinewidth(1)\",color=black];\n graph [concentrate=\"false\",label=\"%d colors\"];\n",m);
  for (i=0; i<g->nnodes; i++) {
    c=g->c[i]/(double)m;
    rainbow(c,1.0,1.0,&rd,&gr,&bl);
    rgb2hexstring(rd,gr,bl,s);
    //fprintf(stderr,"%d %d %f %s\n",i,g->c[i],(g->c[i])/(double)m,s);
    fprintf(f,"  %d [color=\"#%s\",fillcolor=\"#%s\"];\n",i,s,s);
  }
  for (i=0; i<g->nnodes; i++) if (g->d[i]>0) { // i has neighbours
    for (k=0; k<g->d[i]; k++) if ((j=g->a[i][k])>i) fprintf(f,"  %d--%d;\n",i,j);
  }
  fprintf(f,"}\n");
  fclose(f);
  printf("%s written: now do \"neato -Tps %s | gv -\"\n",fn,fn);
}

// non-recursive depth-first search - recommended over graph_dfs...
static void graph_dfs_nonrec(graph_t g, node_t s, node_t l) {
  size_t sz=g->nnodes,usz=1024;
  node_t i,j,nhbr,*stack=malloc(sz*sizeof(node_t));
  int *u,c,top=0,utop,ok;
  u=malloc(usz*sizeof(int));
  stack[0]=s;
  while (top>=0) {
    j=stack[top--]; // pop
    cluster(g,j)=l; // cluster number
    visited(g,j)=1; // we have now visited node i
    utop=0;
    for (i=0; i<degree(g,j); i++) {
      nhbr=neighbour(g,j,i);
      if (!visited(g,nhbr)) { // push
        if (utop>=usz-1) u=realloc(u,(usz+=1024)*sizeof(int));
        u[utop++]=color(g,nhbr);
        if (top>=sz-1) stack=realloc(stack,(sz+=1024)*sizeof(node_t));
        stack[++top]=nhbr; // push
      }
    }
    // find smallest available color for node j
    c=0;
    while (1) {
      ok=1;
      for (i=0; i<utop; i++) if (c==u[i]) { ok=0; break; }
      if (ok) {
        color(g,j)=c;
        break;
      }
      c++;
    }
  }
  free(u);
  free(stack);
}

// recursive depth-first search - may get stack overflow problems...
// OBSOLETE
static void graph_dfs(graph_t g, node_t s, node_t l) {
  node_t i;
  int m=-1;
  node_t *gas=g->a[s];
  //fprintf(stderr,"s=%d l=%d kkk=%lu\n",s,l,kkk++);
  cluster(g,s)=l; // cluster number
  visited(g,s)=1; // we have now visited node s
  // find biggest color used by neighbours
  for (i=0; i<degree(g,s); i++) {
    int j;
    if ((j=color(g,gas[i]))>=0) m=j>m?j:m;
  }
  if (m==-1)      // no neighbours are yet colored
    color(g,s)=0;
  else if (m==0)  // one neighbours colored 0
    color(g,s)=1;
  else { // find colors used by neighbours - they are in the range [0,m]
    // find smallest color not used by neighbours
    int j;
    //char *u=calloc(m+2,sizeof(char)); // extra 0 element as sentinel
    char u[m+2]; bzero(u,(m+2)*sizeof(char)); // extra 0 element as sentinel
    for (i=0; i<degree(g,s); i++) if ((j=color(g,gas[i]))>=0) u[j]=1;
    i=-1; while (u[++i]);
    color(g,s)=i;
    //free(u);
  }
  // visit neighbours
  for (i=0; i<degree(g,s); i++) if (!visited(g,gas[i])) graph_dfs(g,gas[i],l);
}

// recursive depth-first search - may get stack overflow problems...
// OBSOLETE
static void graph_dfs_stack(graph_t g, node_t s, int l) {
  // undocumented - uses stack for temporary arrays (small graphs only)
  int i,j,m=-1;
  node_t *gas=g->a[s];
  g->l[s]=l; // cluster number
  g->v[s]=1; // we have just visited node s
  // find biggest color used by neighbours
  for (i=0; i<g->d[s]; i++) if ((j=color(g,gas[i]))>=0) m=j>m?j:m;
  if (m==-1)      // no neighbours are yet colored
    g->c[s]=0;
  else if (m==0)  // one neighbours colored 0
    g->c[s]=1;
  else { // find colors used by neighbours - they are in the range [0,m]
    int u[m+2]; bzero(u,(m+2)*sizeof(int)); // extra 0 element as sentinel
    for (i=0; i<g->d[s]; i++) if ((j=color(g,gas[i]))>=0) u[j]=1;
    // find smallest color not used by neighbours
    i=0; while (u[i++]); color(g,s)=i-1;
  }
  // visit neighbours
  for (i=0; i<g->d[s]; i++) if (!visited(g,gas[i])) graph_dfs_stack(g,gas[i],l);
}

int graph_greedy_color(graph_t g, int perm[]) {
  // greedy coloring
  int i,pi,c=0;
  for (i=0; i<g->nnodes; i++) { g->v[i]=0; g->c[i]=-1; } // -1=not colored
  for (i=0; i<g->nnodes; i++) {
    if (perm) pi=perm[i]; else pi=i;
    if (!visited(g,pi)) graph_dfs_nonrec(g,pi,c++);
  }
  c=-1;
  for (i=0; i<g->nnodes; i++) if (color(g,i)>c) c=color(g,i);
  return c+1; // colors used
}

int graph_sequential_color(graph_t g, int perm[], int ub) {
  // N[i][j] => node i has a neighbour colored j
  // colored in order given by perm
  // return -1 as soon as more colors than ub used
  int i,j=0,k,m=-1,pi,**N=Malloc(g->nnodes*sizeof(int*));
  int *Npi;
  node_t *gapi;
  for (i=0; i<g->nnodes; i++) N[i]=Calloc(g->nnodes,sizeof(int));
  for (i=0; i<g->nnodes; i++) {
    if (perm) pi=perm[i]; else pi=i;
    Npi=N[pi];
    gapi=g->a[pi];
    j=-1; while (Npi[++j]);
    for (k=0; k<g->d[pi]; k++) N[gapi[k]][j]=1;
    g->c[pi]=j;
    if (j+1>ub) { m=-2; goto quit; };
    m=j>m?j:m;
  }
  quit: for (i=0; i<g->nnodes; i++) Free(N[i],g->nnodes*sizeof(int));
  Free(N,g->nnodes*sizeof(int*));
  return m+1; // number of colors used
}

//! return an upper bound to chi
int graph_sequential_color_repeat(graph_t g, int n) {
  int k,perm[g->nnodes],j,chi_ub=INT_MAX;
  for (j=0; j<g->nnodes; j++) perm[j]=j;
  for (j=0; j<n; j++) {
    k=graph_sequential_color(g,perm,chi_ub);
    if (k>0 && k<chi_ub) chi_ub=k; // got a better bound
    shuffle(perm,g->nnodes,g->nnodes);
  }
  return chi_ub;
}

int graph_check_coloring(graph_t g) {
  int i,j,k,ci;
  for (i=0; i<g->nnodes; i++) {
    ci=g->c[i];
    j=0;
    for (k=0; k<g->d[i]; k++) // for each neighbour
      if (ci==g->c[g->a[i][k]]) { // something wrong
          //for (j=0; j<g->nnodes; j++) printf("color(%d)=%d\n",j,g->c[j]);
        return 0;
      }
  }
  return 1;  // ok
}

void graph_gnp_slow(graph_t g, double p) {
  int i,j;
  graph_empty(g);
  for (i=0; i<g->nnodes; i++) for (j=i+1; j<g->nnodes; j++)
    if (rand()<p*RAND_MAX) graph_append_edge(g,i,j);
}

double binomial_deviate(double pp, int n) {
  int j;
  static int nold=-1;
  double am,em,g,p,bnl,sq,t,y;
  static double pold=-1.0,pc,plog,pclog,en,oldg;
  p=(pp<=0.5?pp:1.0-pp);
  am=n*p;
  if (n<30) {
    bnl=0.0; 
    for (j=1; j<=n; j++) if (UNIFORM_01<p) ++bnl;
  } else if (am<1.0) {
    g=exp(-am);
    for (t=1.0,j=0; j<=n; j++) { t*=UNIFORM_01; if (t<g) break; }
    bnl=(j<=n?j:n);
  } else {
  if (n!=nold) {
    en=n; oldg=lgamma(en+1.0); nold=n;
  } 
  if (p!=pold) {
    pc=1.0-p; plog=log(p); pclog=log(pc); pold=p;
  }
  sq=sqrt(2.0*am*pc);
  do {
    do {
      y=tan(M_PI*UNIFORM_01); em=sq*y+am;
    } while (em<0.0 || em>=(en+1.0));
    em=floor(em);
    t=1.2*sq*(1.0+y*y)*exp(oldg-lgamma(em+1.0)-lgamma(en-em+1.0)+em*plog+(en-em)*pclog);
    } while (UNIFORM_01>t);
    bnl=em; 
  }
  if (p!=pp) return n-bnl;
  return bnl;
}

void graph_gnp(graph_t g, double p) {
  if (p>0.5) // dense - FIXME what is the best cross-over value??
    graph_gnp_slow(g,p);
  else { // sparse
    int i,k=0,n2=nnodes(g)*(nnodes(g)-1)/2;
    double log1mp=log1p(-p);
    graph_empty(g);
    while (1) {
      if ((i=(int)floor(log1p(-UNIFORM_01)/log1mp))) { // geometric deviate
        k+=i; 
        if (k>=n2) return; 
        graph_add_edge1(g,k);
      }
    }
  }
}

void graph_gnm_slow(graph_t g, unsigned long m) {
/* KMB 2006 Dec 04
  OBSOLETE - use graph_gnm
  Return the random graph G_{n,m}.
  Gives a graph picked uniformly at randomly out of the set of all graphs
  with n nodes and m edges.
  Parameters:
   - m: the number of edges
  Algorithm by Keith M. Briggs - python original Mar 31, 2006.
  Inspired by Knuth's Algorithm S (Selection sampling technique),
  in section 3.4.2 of The Art of Computer Programming by Donald E. Knuth
  Volume 2 / Seminumerical algorithms; Third Edition, Addison-Wesley, 1997.
*/
  int u=0,v=1,n=nnodes(g);
  unsigned long k=0,t=0;
  double mmax=(n*(n-1.0))/2; // beware overflow
  if (n==1) return;
  if (m>=mmax) {
    graph_complete(g);
    return;
  }
  graph_empty(g);
  while (1) {
    if ((mmax-t)*(rand()/(1.0+RAND_MAX))<m-k) {
      assert(u<n);
      assert(v<n);
      graph_append_edge(g,u,v);
      if (++k==m) return;
    }
    t++;
    v++;
    if (v==n) {
      u++;
      v=u+1;
    }
  }
}

// helper functions for Vitter's fast random sampling method
// ACM Trans Math Soft 13, 58 (1987)

inline static void graph_add_edge1(graph_t g, double k) {
  // add an edge, counted by k in 0..n*(n-1)/2-1
  // we will map k to (i,j) as follows, to get an entry in the lower triangle
  // of the adjacency matrix.
  // example for nnodes=5, i=row, j=column, k=...
  // x x x x x
  // 0 x x x x
  // 1 2 x x x
  // 3 4 5 x x
  // 6 7 8 9 x
  // k=i*(i-1)/2+j, so the inverse function is...
  // i=floor((1+sqrt(1+8*k))/2)
  // j=k-i*(i-1)/2
  // could use iteration i<-floor((i^2+2*k)/(2*i-1)) instead of sqrt here...
  node_t i=(node_t)floor((1.0+sqrt(1.0+8*k))/2.0),j=(node_t)(k-(i*(i-1))/2);
  //fprintf(stderr,"%.0f->(%lu,%lu)\n",k,(long unsigned int)i,(long unsigned int)j);
  graph_append_edge(g,i,j); // can append - each possible k will occur at most once
}

static void Method_A_double(double nr, double Nr, double offset, graph_t g) {
  // ACM Trans Math Soft 13, 58 (1987)
  // Method A page 65.  For n/N large.
  // offset is added to each output
  // do everything in double to get 2^53 range
  double v,quot,top=Nr-nr,s,k=-1.0;
  while (nr>1) {
    v=VN_GRAPH_UNIFORM_01;
    s=0.0;
    quot=top/Nr;
    while (quot>v) {
      s++;
      top--;
      Nr--;
      quot*=top/Nr;
    }
    //if (s==0.0) fprintf(stderr,"!! Nr=%g nr=%g k=%g quot=%g v=%g\n",Nr,nr,k,quot,v),exit(1);
    k+=s+1; // skip
    graph_add_edge1(g,k+offset);
    Nr--;
    nr--;
  }
  k+=1+floor(Nr*VN_GRAPH_UNIFORM_01);
  graph_add_edge1(g,k+offset);
}

static void Method_D_double(double nr, double Nr, graph_t g) {
  // ACM Trans Math Soft 13, 58 (1987)
  // Method D page 66.  For n/N small.
  // do everything in double to get 2^53 range
  double ninv=1.0/nr,nmin1inv,u,x,vprime,y1,y2,top,bot,qu1r=Nr-nr+1.0,k=-1.0,s,t,limit,negalphainv=-13.0,threshold=-negalphainv*nr;
  vprime=pow(VN_GRAPH_UNIFORM_01,ninv);
  while (nr>1 && threshold<Nr) {
    nmin1inv=1.0/(nr-1.0);
    while (1) {
      while (1) {
        x=Nr*(1.0-vprime);
        s=floor(x);
        if (s<qu1r) break;
        vprime=pow(VN_GRAPH_UNIFORM_01,ninv);
      }
      u=VN_GRAPH_UNIFORM_01;
      y1=pow(u*Nr/qu1r,nmin1inv);
      vprime=y1*(1.0-x/Nr)*(qu1r/(qu1r-s));
      if (vprime<=1.0) break;
      y2=1.0;
      top=Nr-1.0;
      if (nr-1>s) {
        bot=Nr-nr;
        limit=Nr-s;
      } else {
        bot=Nr-s-1.0;
        limit=qu1r;
      }
      for (t=Nr-1; t>=limit; t--) {
        y2*=top/bot;
        top--;
        bot--;
      }
      if (Nr/(Nr-x)>=y1*pow(y2,nmin1inv)) {
        vprime=pow(VN_GRAPH_UNIFORM_01,nmin1inv);
        break;
      }
      vprime=pow(VN_GRAPH_UNIFORM_01,ninv);
    }
    k+=s+1; // skip
    graph_add_edge1(g,k);
    Nr-=s+1.0;
    nr--;
    ninv=nmin1inv;
    qu1r-=s;
    threshold+=negalphainv;
  }
  if (nr>1) {
    //fprintf(stderr,"A: Nr=%g k=%g nr=%g offset=%g\n",Nr,k,nr,k);
    if (k==-1) Method_A_double(nr,Nr,0,g); // outside threshold
    else       Method_A_double(nr,Nr-k,k,g);
  } else {
    s=1+floor(Nr*vprime);
    k+=s; // skip
    graph_add_edge1(g,k);
  }
}

void graph_gnm(graph_t g, unsigned long m) {
/* KMB 2006 Dec 11
  Return a random graph from G_{n,m}.
  Gives a graph picked uniformly at random out of the set of all graphs
  with n nodes and m edges.
  Parameters:
   - m: the number of edges
*/
  graph_empty(g);
  double nr=m,Nr=g->nnodes*(g->nnodes-1.0)/2.0;
  //Method_A_double(nr,Nr,0,g); // slow
  Method_D_double(nr,Nr,g); // fast!
}

void graph_grg(graph_t g, double r) {
  // slow
  int i,j,n=g->nnodes;
  double r2=r*r,dx,dy,x[n],y[n];
  graph_empty(g);
  for (i=0; i<n; i++) {
    x[i]=UNIFORM_01;
    y[i]=UNIFORM_01;
  }
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) {
      dx=x[i]-x[j];
      dy=y[i]-y[j];
      if (dx*dx+dy*dy<r2) graph_append_edge(g,i,j);
  }
}

typedef struct { double x,y; } xy;

static int cmp(const void *a, const void *b) {
  xy *aa=(xy*)a;
  xy *bb=(xy*)b;
  return aa->x<bb->x?-1:(aa->x>bb->x?1:0);
}

static double torus_distance2(double dx, double dy) {
  if (dx>0.5) dx=1.0-dx;
  if (dy>0.5) dy=1.0-dy;
  return dx*dx+dy*dy;
}

void graph_grg_torus(graph_t g, double r) {
  int i,j,n=g->nnodes;
  double r2=r*r,x[n],y[n];
  graph_empty(g);
  for (i=0; i<n; i++) {
    x[i]=UNIFORM_01;
    y[i]=UNIFORM_01;
  }
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++)
      if (torus_distance2(fabs(x[i]-x[j]),fabs(y[i]-y[j]))<r2) graph_append_edge(g,i,j);
}

void graph_lognormal_grg_torus(graph_t g, double r, double alpha) {
  // nodes d apart are joined with probability (1-erf(log(d/r)/alpha))/2
  // so alpha=0 corresponds to GRG(n,r), but alpha>0 makes the cutoff fuzzier.
  // log(d/r)/alpha = log(sqrt(d^2)/r)/alpha = (log(d^2)/2-log(r))/alpha
  int i,j,n=g->nnodes;
  if (alpha==0.0) {
    graph_grg_torus(g,r);
    return;
  }
  double d2,x[n],y[n],c=log(r);
  graph_empty(g);
  for (i=0; i<n; i++) {
    x[i]=UNIFORM_01;
    y[i]=UNIFORM_01;
  }
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) {
      d2=torus_distance2(fabs(x[i]-x[j]),fabs(y[i]-y[j]));
      if (UNIFORM_01<0.5*(1.0-erf((log(d2)/2-c)/alpha)))
        graph_append_edge(g,i,j);
  }
}

void graph_lognormal_grg(graph_t g, double r, double alpha) {
  // nodes d apart are joined with probability (1-erf(log(d/r)/alpha))/2
  // so alpha=0 corresponds to GRG(n,r), but alpha>0 makes the cutoff fuzzier.
  // log(d/r)/alpha = log(sqrt(d^2)/r)/alpha = (log(d^2)/2-log(r))/alpha
  int i,j,n=g->nnodes;
  if (alpha==0.0) {
    graph_grg(g,r);
    return;
  }
  double d2,x[n],y[n],c=log(r);
  graph_empty(g);
  for (i=0; i<n; i++) {
    x[i]=UNIFORM_01;
    y[i]=UNIFORM_01;
  }
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) {
      d2=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]);
      if (UNIFORM_01<0.5*(1.0-erf((log(d2)/2-c)/alpha)))
        graph_append_edge(g,i,j);
  }
}

void graph_grg_torus_sorted(graph_t g, double r) {
  int i,j,m=0,n=g->nnodes;
  double r2=r*r;
  xy z[n];
  for (i=0; i<n; i++) {
    z[i].x=UNIFORM_01;
    z[i].y=UNIFORM_01;
  }
  qsort(z,n,sizeof(xy),cmp); // sort on x
  graph_empty(g);
  for (i=0; i<n; i++) {
    for (j=i+1; j<n; j++)
      if (torus_distance2(fabs(z[i].x-z[j].x),fabs(z[i].y-z[j].y))<r2) graph_append_edge(g,i,j),m++;
  }
}

int graph_grg_torus_iterator(unsigned int n, double r, int f(graph_t,int,void*), void* cd) {
  /* 
  Iterate over the ensemble grg(n,r).
  Strategy: divide unit square into nb*nb boxes of width>=r.  Two nodes 
  can only be neighbours if they are in the same or abutting boxes. 
  We then only need on average to look at 
    (n/(nb*nb)*(n/(nb*nb)-1)/2)*(5*nb*nb)
  pairs of nodes, compared to n*(n-1)/2 the slow way.
  This method wins when (n/(nb*nb)*(n/(nb*nb)-1)/2)*(5*nb*nb)>n*(n-1)/2,
  which gives nb>sqrt(5) when n is large.
  Since nb~1/r, the method checks about 5/2*n^2*r^2 pairs.
  f() is called for each graph, and should return a zero value to
  terminate the iteration.  Client data may be passed to f through cd.
  Revised 2007 Feb 02.
  */
  const int timing=0;
  int i,j,k,l,m,ig=0;
  double ave,stdev,*x,*y,r2;
  size_t initial_length;
  unsigned int nb,***boxes,**count;
  size_t **alloc;
  /* timing  */
  struct tms buffer;
  int tps; 
  clock_t current_time,start_time;
  // internal functions...
  int i0,i1,k0,k1;
  double dx,dy;
  graph_t g=graph_new(n);
  if (r<=0.0) {
    graph_empty(g);
    do { // iterate over graphs
    } while ((l=f(g,ig++,cd)));
    return l;
  }
  inline void check_box(int i, int j) {
    // check all pairs in box (i,j)
    for (i0=0; i0<count[i][j]; i0++) {
      k0=boxes[i][j][i0];
      for (i1=0; i1<count[i][j]; i1++) {
        k1=boxes[i][j][i1];
        if (k0<k1) {
  	  dx=fabs(x[k0]-x[k1]);
          dy=fabs(y[k0]-y[k1]);
          if (torus_distance2(dx,dy)<r2) graph_append_edge(g,k0,k1);
        }
      }
    }
  }
  inline void check_pairs(int i, int j, int l, int m) {
    // check all pairs, one in box (i,j), one in box (l,m)
    // (l,m) will always be one of the boxes above, above right, right, 
    // or below right of (i,j).   Thus we get each potential edge once only. 
    for (i0=0; i0<count[i][j]; i0++) {
      k0=boxes[i][j][i0];
      for (i1=0; i1<count[l][m]; i1++) {
        k1=boxes[l][m][i1];
        if (1) {
  	  dx=fabs(x[k0]-x[k1]);
          dy=fabs(y[k0]-y[k1]);
          if (torus_distance2(dx,dy)<r2) graph_append_edge(g,k0,k1);
        }
      }
    }
  }
  nb=(int)floor(1.0/r);
  r2=r*r;
  x=malloc(n*sizeof(double));
  y=malloc(n*sizeof(double));
  //> lprint(solve((n/(nb*nb)*(n/(nb*nb)-1)/2)*(5*nb*nb)=n*(n-1)/2,nb));
  //1/(4+n)*5^(1/2)*((4+n)*n)^(1/2), -1/(4+n)*5^(1/2)*((4+n)*n)^(1/2)
  if (nb<4 || nb<1.0/(4.0+n)*sqrt(5*(4.0+n)*n)) { // use 'slow' method
    do { // iterate over graphs
      graph_empty(g);
      for (k=0; k<n; k++) {
        x[k]=UNIFORM_01;
        y[k]=UNIFORM_01;
      }
      for (i=0; i<n; i++)
        for (j=i+1; j<n; j++)
          if (torus_distance2(fabs(x[i]-x[j]),fabs(y[i]-y[j]))<r2)
            graph_append_edge(g,i,j);
    } while ((l=f(g,ig++,cd)));
    free(x);
    free(y);
    graph_clear(g);
    return l;
  }
  if (timing) {
    tps=sysconf(_SC_CLK_TCK); // clock ticks per second
    times(&buffer);
    start_time=buffer.tms_utime;
  }
  boxes=malloc(nb*sizeof(unsigned int**));
  count=malloc(nb*sizeof(unsigned int*));
  alloc=malloc(nb*sizeof(size_t*));
  for (i=0; i<nb; i++) {
    boxes[i]=malloc(nb*sizeof(unsigned int*));
    count[i]=calloc(nb,sizeof(unsigned int)); // initially zero
    alloc[i]=malloc(nb*sizeof(size_t)); // allocated length of list boxes[i][j]
  }
  ave=n/((double)nb*nb);  // expected contents of each box
  stdev=1.0/sqrt(ave);    // typical standard deviation of box contents
  initial_length=(size_t)ceil(ave*(1.0+6.0*stdev)); // hope to never have to reallocate
  if (initial_length<8) initial_length=8;
  for (i=0; i<nb; i++) for (j=0; j<nb; j++) {
    alloc[i][j]=initial_length;
    boxes[i][j]=malloc(initial_length*sizeof(unsigned int));
  }
  if (timing) {
    times(&buffer);
    current_time=buffer.tms_utime;
    printf("setup time=%g\n",(current_time-start_time)/(double)tps);
    start_time=buffer.tms_utime;
  }
  do { // iterate over graphs
    graph_empty(g);
    for (k=0; k<n; k++) { // generate data and count box contents
      x[k]=UNIFORM_01;
      y[k]=UNIFORM_01;
      i=(int)floor(nb*x[k]);
      j=(int)floor(nb*y[k]); // node k is in box (i,j)
      if (count[i][j]>=alloc[i][j]) {
        alloc[i][j]=4+count[i][j];
        boxes[i][j]=realloc(boxes[i][j],alloc[i][j]*sizeof(unsigned int));
      }
      boxes[i][j][count[i][j]++]=k;
    }
    // loop over all pairs of abutting boxes
    for (i=0; i<nb; i++) for (j=0; j<nb; j++) {
      check_box(i,j);                   // self
      check_pairs(i,j,i,m=(j+1)%nb);    // above
      check_pairs(i,j,l=(i+1)%nb,m);    // top right
      check_pairs(i,j,l,j);             // right
      check_pairs(i,j,l,(j+nb-1)%nb);   // bottom right
    }
    for (i=0; i<nb; i++) for (j=0; j<nb; j++) count[i][j]=0;
  } while ((l=f(g,ig++,cd)));
  if (timing) {
    times(&buffer);
    current_time=buffer.tms_utime;
    printf("iterate time=%g\n",(current_time-start_time)/(double)tps);
    printf("fast: %g seconds for %d graphs\n",(current_time-start_time)/(double)tps,ig);
  }
  for (i=0; i<nb; i++) {
    for (j=0; j<nb; j++) free(boxes[i][j]);
    free(boxes[i]);
    free(count[i]);
    free(alloc[i]);
  }
  free(count);
  free(boxes);
  free(alloc);
  free(x);
  free(y);
  graph_clear(g);
  return l; // pass back last return value of f
}

// the next two functions have no real speed advantage (unlike graph_grg_torus_iterator) but are provided for uniformity of interface

int graph_gnm_iterator(unsigned int n, unsigned long m, int f(graph_t,int,void*), void* cd) {
  int l,ig=0;
  graph_t g=graph_new(n);
  do { // iterate over graphs
    graph_gnm(g,m);
  } while ((l=f(g,ig++,cd)));
  graph_clear(g);
  return l;
}

int graph_gnp_iterator(unsigned int n, double p, int f(graph_t,int,void*), void* cd) {
  int l,ig=0;
  graph_t g=graph_new(n);
  do { // iterate over graphs
    graph_gnp(g,p);
  } while ((l=f(g,ig++,cd)));
  graph_clear(g);
  return l;
}

int graph_random_line_graph_iterator(unsigned int nlines, int f(graph_t,int,void*), void* cd) {
  // model: 
  //  1. place nlines dots uniformly in the unit square
  //  2. put a line of random slope through each dot
  //  3. intersections of these lines are the modes of the graph
  //  4. line segments between intersections are the edges of the graph
  //  See R E Miles PNAS 52, 901 (1964), PNAS 52, 1157 (1964)
  // python originals: regular_degree4.py, random_planar.py
  typedef struct {
    double x;
    int n;
  } xs_cell;
  int cmp_xs(const void *ap, const void *bp) {
    double a=(*(xs_cell*)ap).x;
    double b=(*(xs_cell*)bp).x;
    if (a<b) return -1;
    if (a>b) return +1;
    return 0;
  }
  int done,i,j,l,k,lq,n,len_xs;
  double x0i,y0i,si,x0j,y0j,sj,x,y;
  double *lines_x0,*lines_y0,*lines_s;
  int **line_contains,*len_lc;
  xs_cell *xs,*q;
  graph_t g;
  // these over-allocate to save the bother of dynamically allocated stacks
  // and lists.  But as we don't expect more than about nlines=20, it should 
  // not matter...
  xs=malloc(nlines*(nlines-1)/2*sizeof(xs_cell));
  q=malloc((nlines-1)*sizeof(xs_cell));
  line_contains=malloc(nlines*sizeof(int*));
  len_lc=malloc(nlines*sizeof(int));
  for (i=0; i<nlines; i++) line_contains[i]=malloc(nlines*sizeof(int));
  lines_x0=malloc(nlines*sizeof(double));
  lines_y0=malloc(nlines*sizeof(double));
  lines_s =malloc(nlines*sizeof(double));
  k=0;
  while (1) {
    // get lines
    for (i=0; i<nlines; i++) {
      lines_x0[i]=UNIFORM_01;
      lines_y0[i]=UNIFORM_01;
      lines_s[i] =tan(M_PI*UNIFORM_01);
      len_lc[i]=0;
    }
    // get intersections...
    for (n=len_xs=i=0; i<nlines; i++) {
      x0i=lines_x0[i];
      y0i=lines_y0[i];
      si =lines_s[i];
      for (j=i+1; j<nlines; j++) {
        x0j=lines_x0[j];
        y0j=lines_y0[j];
        sj =lines_s[j];
        x=(sj*x0j-si*x0i+y0i-y0j)/(sj-si); // x coordinate of intersection
        if (0.0<=x && x<1.0) {
          y=si*(x-x0i)+y0i; // y coordinate of intersection
          if (0.0<=y && y<1.0) {
            line_contains[i][len_lc[i]++]=n; // lines i and j contain point n
            line_contains[j][len_lc[j]++]=n; // lines i and j contain point n
            xs[len_xs].x=x;
            xs[len_xs].n=n;
            len_xs++; // FIXME do we need both len_xs and n?
            n++;
          }
        }
      }
    }
    g=graph_new(n);
    for (l=0; l<nlines; l++) {
      for (lq=i=0; i<len_lc[l]; i++) q[lq++]=xs[line_contains[l][i]];
      qsort(q,lq,sizeof(xs_cell),cmp_xs);
      for (i=0; i<lq-1; i++) graph_add_edge(g,q[i].n,q[i+1].n);
    }
    done=(f(g,k++,cd)==0);
    graph_clear(g);
    if (done) break;
  }
  free(lines_x0);
  free(lines_y0);
  free(lines_s);
  free(q);
  free(xs);
  free(len_lc);
  for (i=0; i<nlines; i++) free(line_contains[i]);
  free(line_contains);
  return 0;
}

int graph_random_line_graph_poisson_iterator(double tau, int f(graph_t,int,void*), void* cd) {
  // See R E Miles PNAS 52, 901 (1964), PNAS 52, 1157 (1964)
  // The mean number of lines in the unit square is 4*tau/pi (Miles thm 3)
  // first allocate enough for largest expected n. 
  // model: 
  //  0. get n=Poisson(4*tau/pi)
  //  1. place n dots uniformly in the unit square
  //  2. put a line of random slope through each dot
  //  3. intersections of these lines are the modes of the graph
  //  4. line segments between intersections are the edges of the graph
  typedef struct {
    double x;
    int n;
  } xs_cell;
  int cmp_xs(const void *ap, const void *bp) {
    double a=(*(xs_cell*)ap).x;
    double b=(*(xs_cell*)bp).x;
    if (a<b) return -1;
    if (a>b) return +1;
    return 0;
  }
  int nlines,nlinesmax;
  int done,i,j,l,k,lq,n,len_xs;
  double x0i,y0i,si,x0j,y0j,sj,x,y;
  double *lines_x0,*lines_y0,*lines_s;
  int **line_contains,*len_lc;
  xs_cell *xs,*q;
  nlinesmax=1+(int)tau;
  while (logpoisson(nlinesmax,tau)>-14) nlinesmax++; // P[n>nlinesmax]<10^-6
  fprintf(stderr,"tau=%g nlinesmax=%d\n",tau,nlinesmax);
  graph_t g;
  // these over-allocate to save the bother of dynamically allocated stacks
  // and lists.  But as we don't expect more than about nlines=20, it should 
  // not matter...
  xs=malloc(nlinesmax*(nlinesmax-1)/2*sizeof(xs_cell));
  q=malloc((nlinesmax-1)*sizeof(xs_cell));
  line_contains=malloc(nlinesmax*sizeof(int*));
  len_lc=malloc(nlinesmax*sizeof(int));
  for (i=0; i<nlinesmax; i++) line_contains[i]=malloc(nlinesmax*sizeof(int));
  lines_x0=malloc(nlinesmax*sizeof(double));
  lines_y0=malloc(nlinesmax*sizeof(double));
  lines_s =malloc(nlinesmax*sizeof(double));
  k=0;
  while (1) {
    do nlines=poisson_deviate(4*tau/M_PI); while (nlines>nlinesmax);
    //if (nlines>0) fprintf(stderr,"tau=%g nlines=%d\n",tau,nlines);
    // get lines
    for (i=0; i<nlines; i++) {
      lines_x0[i]=UNIFORM_01;
      lines_y0[i]=UNIFORM_01;
      lines_s[i] =tan(M_PI*UNIFORM_01);
      len_lc[i]=0;
    }
    // get intersections...
    for (n=len_xs=i=0; i<nlines; i++) {
      x0i=lines_x0[i];
      y0i=lines_y0[i];
      si =lines_s[i];
      for (j=i+1; j<nlines; j++) {
        x0j=lines_x0[j];
        y0j=lines_y0[j];
        sj =lines_s[j];
        x=(sj*x0j-si*x0i+y0i-y0j)/(sj-si); // x coordinate of intersection
        if (0.0<=x && x<1.0) {
          y=si*(x-x0i)+y0i; // y coordinate of intersection
          if (0.0<=y && y<1.0) {
            line_contains[i][len_lc[i]++]=n; // lines i and j contain point n
            line_contains[j][len_lc[j]++]=n; // lines i and j contain point n
            xs[len_xs].x=x;
            xs[len_xs].n=n;
            len_xs++; // FIXME do we need both len_xs and n?
            n++;
          }
        }
      }
    }
    g=graph_new(n);
    for (l=0; l<nlines; l++) {
      for (lq=i=0; i<len_lc[l]; i++) q[lq++]=xs[line_contains[l][i]];
      qsort(q,lq,sizeof(xs_cell),cmp_xs);
      for (i=0; i<lq-1; i++) graph_add_edge(g,q[i].n,q[i+1].n);
    }
    done=(f(g,k++,cd)==0);
    graph_clear(g);
    if (done) break;
  }
  free(lines_x0);
  free(lines_y0);
  free(lines_s);
  free(q);
  free(xs);
  free(len_lc);
  for (i=0; i<nlinesmax; i++) free(line_contains[i]);
  free(line_contains);
  return 0;
}

int graph_ncolors(graph_t g) {
  int i,mc=-1;
  for (i=0; i<g->nnodes; i++) if (g->c[i]>mc) mc=g->c[i];
  return mc+1;
}

int graph_nclusters(graph_t g) {
  int i,ml=-1;
  graph_greedy_color(g,NULL);
  for (i=0; i<g->nnodes; i++) if (g->l[i]>ml) ml=g->l[i];
  return ml+1;
}

int graph_connected(graph_t g) {
  return graph_nclusters(g)==1?1:0;
}

static int int_sort_down(const void *ap, const void *bp) {
  int a=*(int*)ap;
  int b=*(int*)bp;
  return a>b?-1:(a==b?0:1);
}

int* graph_cluster_sizes(graph_t g) {
  // return a sorted list of cluster sizes, biggest first.
  // A zero terminates the list.
  int i,nc=graph_nclusters(g);
  int *cs=calloc(nc+1,sizeof(int));
  for (i=0; i<nnodes(g); i++) cs[cluster(g,i)]++;
  qsort(cs,nc,sizeof(int),int_sort_down);
  cs[nc]=0; // sentinel
  return cs;
}

int graph_max_cluster(graph_t g) {
  // return size of largest cluster
  int i,m=0,nc=graph_nclusters(g);
  int cs[nc];
  bzero(cs,nc*sizeof(int));
  for (i=0; i<nnodes(g); i++) cs[cluster(g,i)]++;
  for (i=0; i<nc; i++) if (cs[i]>m) m=cs[i];
  return m;
}

//! in-place local_complementation at node i; i.e. subgraph N_i(g) is complemented, where N_i(g) is the set of neighbours of i
void graph_local_complement(graph_t g, node_t i) {
  const int verbose=0;
  int a,d,nd=0,na=0,n=g->d[i]*(g->d[i]-1)/2;
  node_t j,k,nj,nk,addj[n],addk[n],delj[n],delk[n];
  for (j=0; j<degree(g,i); j++) {
    nj=g->a[i][j];
    for (k=0; k<degree(g,i); k++) {
      nk=g->a[i][k];
      if (nj<nk) { 
        if (graph_has_edge(g,nj,nk)) { 
          delj[nd]=nj; delk[nd]=nk; nd++;
        } else { 
          addj[na]=nj; addk[na]=nk; na++;
        }
      }
    }
  }
  if (verbose) {
    printf("N_%d(g): ",i);
    for (d=0; d<nd; d++) printf("%d--%d ",delj[d],delk[d]);
    printf("\nna=%d nd=%d\n",na,nd);
  }
  for (a=0; a<na; a++) graph_add_edge(g,addj[a],addk[a]);
  for (d=0; d<nd; d++) graph_del_edge(g,delj[d],delk[d]);
}

// histogram functions

histogram_t histogram_new(char* lbl) {
  int sz=1024;
  histogram_t h=Malloc(sizeof(struct _histogram));
  h->h =Calloc(sz,sizeof(double));
  h->m =Calloc(sz,sizeof(double));
  h->ss=Calloc(sz,sizeof(double));
  h->a=sz;
  h->n=0;
  h->lbl=lbl;
  return h;
}

void histogram_empty(histogram_t h) {
  memset(h->h, 0,h->a*sizeof(double));
  memset(h->m, 0,h->a*sizeof(double));
  memset(h->ss,0,h->a*sizeof(double));
  h->n=0;
}

void histogram_add(histogram_t h, int i) {
  double om,x;
  if (i<0) return;
  if (i>=h->a-1) {
    int j=h->a;
    h->a=i+1;
    h->h =Realloc(h->h, i*sizeof(double),h->a*sizeof(double));
    h->m =Realloc(h->m, i*sizeof(double),h->a*sizeof(double));
    h->ss=Realloc(h->ss,i*sizeof(double),h->a*sizeof(double));
    for (; j<h->a; j++) h->h[j]=h->m[j]=h->ss[j]=0.0;
  }
  h->n++;
  h->h[i]++;
  // get the mean proportion so far in bin i...
  x=h->h[i]/(double)h->n; // the value being included in the mean and ss
  om=h->m[i];
  h->m[i]+=(x-om)/h->h[i];
  h->ss[i]+=(x-h->m[i])*(x-om);
}

void histogram_clear(histogram_t h) {
  Free(h->h,0);
  Free(h->m,0);
  Free(h->ss,0);
  Free(h,0);
}

int histogram_max(histogram_t h) {
  // maximum value in h
  int i=h->a;
  while (i>0 && 0==h->h[--i]);
  return i;
}

void histogram_write(FILE* f, histogram_t h) {
  int i;
  int i0,i1,imx=-1,mx=-1;
  double se,sum=0.0,ave=0.0,ss=0.0;
  i0=-1;   while (i0<h->a-1 && !h->h[++i0]);
  i1=h->a; while (i1>i0     && !h->h[--i1]);
  for (i=i0; i<=i1; i++) {
    if ((int)h->h[i]>mx) { mx=h->h[i]; imx=i; }
    sum+=h->h[i];
    ave+=h->h[i]*i;
  }
  if (sum) {
    ave/=sum;
    fprintf(f,"%s %ld items; ",h->lbl,(long)sum);
    fprintf(f,"mode=%d",imx);
    fprintf(f," mean=%.6f",ave);
    for (i=i0; i<=i1; i++) ss+=h->h[i]*(i-ave)*(i-ave);
    fprintf(f," var=%.4g",ss/(sum-1));
    fprintf(f," stdev=%.4g",sqrt(ss/(sum-1)));
    fprintf(f," stderr=%.4g\n",sqrt(ss/(sum-1))/sqrt(sum));
    for (i=i0; i<=i1; i++) {
      se=sqrt(h->ss[i]/(h->h[i]-1))/sqrt(h->h[i]); // standard error of the mean
      fprintf(f,"%2d %10lu   %8.6f   %6.3f%% se=%.3g\n",
        i,(unsigned long)h->h[i],h->h[i]/(double)sum,100*h->h[i]/(double)sum,se);
    }
  }
}

void histogram_show(histogram_t h) {
  histogram_write(stdout,h);
}

double histogram_quantile(const histogram_t h, double y) {
  double cum(double x) { // cumulative distribution from h
    int i=0;
    double s=0.0,x0,x1,y0,y1;
    x0=floor(x-0.5)+0.5;
    if (x0>=h->a) x0=h->a-1;
    for (i=0; i<x0; i++) s+=h->h[i];
    x1=x0+1.0;
    y0=s;
    y1=s+h->h[i];
    return ((y1-y0)/(x1-x0)*(x-x0)+y0)/h->n; // linear interpolation
  }
  // bisect
  double md,dn=0,up=histogram_max(h);
  md=0.5*(up+dn);
  while (fabs(up-dn)>1e-9) {
    if (cum(md)>y) up=md; else dn=md;
    md=0.5*(up+dn);
  }
  return md;
}

double histogram_mean(const histogram_t h) {
  unsigned int i;
  double m=0.0;
  for (i=0; i<h->a; i++) m+=h->h[i]*i;
  return m/h->n;
}

double histogram_variance(const histogram_t h) { // FIXME use h->ss
  unsigned int i;
  double m=0.0,s=0.0;
  for (i=0; i<h->a; i++) {
    m+=h->h[i]*i;
    s+=h->h[i]*i*i;
  }
  return (s-m*m/h->n)/(h->n-1);
}

unsigned int histogram_mode(const histogram_t h) {
  unsigned int i,imx=0;
  double mx=0.0;
  for (i=0; i<h->a; i++) if (h->h[i]>mx) { mx=h->h[i]; imx=i; }
  return imx;
}

int histogram_biggest_done(const histogram_t h, double eps) {
  // true if relative standard error of the mean of the highest bin is<eps
  unsigned int i,imx=0;
  double mx=0.0;
  for (i=0; i<h->a; i++) if (h->h[i]>mx) { mx=h->h[i]; imx=i; }
  return sqrt(h->ss[imx])<eps*(1+mx);
}

int histogram_biggest2_done(const histogram_t h, double eps) {
  // true if relative standard error of the mean of the two highest bins is<eps
  unsigned int i,imx0=0,imx1=0;
  double mx0=0.0,mx1=0.0;
  for (i=0; i<h->a; i++) if (h->h[i]>mx0) { mx0=h->h[i]; imx0=i; }
  if (sqrt(h->ss[imx0])<eps*(1+mx0)) { // maximum converged
    for (i=0; i<h->a; i++) // find 2nd biggest
      if (i!=imx0 && h->h[i]>mx1) { mx1=h->h[i]; imx1=i; }
    return sqrt(h->ss[imx1])<eps*(1+mx1);
  }
  return 0;
}

int histogram_mean_done(const histogram_t h, double eps) {
  // heuristic: true if standard error of the mean is<eps*(1+mean)
  unsigned int i;
  double x,m=0.0,s=0.0,v;
  if (h->n<10) return 0;
  for (i=0; i<h->a; i++) {
    x=h->h[i]*i;
    m+=x;
    s+=x*i;
  }
  //fprintf(stderr,"m=%g s=%g %g; ",m,s,m*m/h->n);
  m/=h->n;                 // mean
  v=(s-h->n*m*m)/(h->n-1); // variance
  s=sqrt(v/h->n);          // standard error of the mean
  //fprintf(stderr,"%g %g %g\n",v,s/sqrt(h->n),eps*(1.0+m));
  return s<eps*(1.0+m);
}

int graph_chromatic_number(graph_t g, clock_t timeout) {
  int c,i,j,k,l,place,*valid,*clique;
  for (i=0; i<g->nnodes; i++) g->c[i]=-1; // in case already colored
  if (g->nedges==0) { // empty graph
    for (i=0; i<g->nnodes; i++) g->c[i]=0;
    return 1;
  }
  if (g->nedges==1) {
    for (i=0; i<g->nnodes; i++) g->c[i]=0;
    for (i=0; i<g->nnodes; i++) { // find the edge
      if (g->d[i]>0) {
        g->c[i]=1;
        return 2;
      }
    }
  }
  if (2*g->nedges==g->nnodes*(g->nnodes-1)) { // complete graph
    for (i=0; i<g->nnodes; i++) g->c[i]=i;
    return g->nnodes;
  }
  #ifdef MALLOC_FREE_CHECK
    alloc=0;
  #endif
  cg.g=g;
  cg.num_node=g->nnodes;
  cg.adj=Malloc(cg.num_node*sizeof(int*));
  cg.ColorAdj=Malloc(cg.num_node*sizeof(int*));
  for (i=0; i<cg.num_node; i++) {
    cg.ColorAdj[i]=Malloc((cg.num_node+1)*sizeof(int));
    cg.adj[i]=Calloc(g->nnodes,sizeof(int));
    for (j=0; j<g->d[i]; j++) cg.adj[i][g->a[i][j]]=1;
  }
  valid        =Malloc(cg.num_node*sizeof(int));
  clique       =Calloc(cg.num_node,sizeof(int));
  cg.ColorClass=Malloc(cg.num_node*sizeof(int));
  cg.Order     =Malloc(cg.num_node*sizeof(int));
  cg.Handled   =Malloc(cg.num_node*sizeof(int));
  cg.ColorCount=Malloc(cg.num_node*sizeof(int));
  for (i=0; i<cg.num_node; i++)
  cg.prob_count=0;
  for (i=0; i<cg.num_node; i++) {
    for (j=0; j<=cg.num_node; j++) cg.ColorAdj[i][j]=0;
    for (j=0; j<cg.num_node; j++) if (cg.adj[i][j]) cg.ColorAdj[i][0]++;
    cg.ColorCount[i]=0;
    cg.Handled[i]=0;
    valid[i]=1;
  }
  cg.BestColoring=cg.num_node+1;
  cg.best_clique=0;
  cg.lb=max_w_clique(valid,clique,0,cg.num_node);
  cg.num_prob=0;
  cg.max_prob=10000;
  // timing...
  cg.tps=sysconf(_SC_CLK_TCK); // clock ticks per second
  times(&cg.buffer);
  cg.start_time=cg.buffer.tms_utime;
  cg.timeout=timeout;
  if (1) {
    place=0;
    for (i=0; i<cg.num_node; i++) if (clique[i]) {
      cg.Order[place]=i;
      cg.Handled[i]=1;
      place++;
      AssignColor(i,place);
      for (j=0; j<cg.num_node; j++)
        if (i!=j && clique[j] && !cg.adj[i][j]) {
          fprintf(stderr,"graph_chromatic_number internal error: result is not a clique!\n");
          graph_show(g);
          for (l=0; l<cg.num_node; l++) printf("clique[%d]=%d\n",l,clique[l]);
          for (l=0; l<cg.num_node; l++) {
            for (k=0; k<cg.num_node; k++) printf("%d ",cg.adj[l][k]);
            printf("\n");
          }
          exit(2);
        }
    }
    c=Color(place,place);
  }
  for (i=0; i<cg.num_node; i++) {
    Free(cg.ColorAdj[i],(cg.num_node+1)*sizeof(int));
    Free(cg.adj[i],g->nnodes*sizeof(int));
  }
  Free(cg.adj,g->nnodes*sizeof(int*));
  Free(cg.ColorAdj,cg.num_node*sizeof(int*));
  Free(cg.ColorClass,cg.num_node*sizeof(int));
  Free(cg.Order,cg.num_node*sizeof(int));
  Free(cg.Handled,cg.num_node*sizeof(int));
  Free(cg.ColorCount,cg.num_node*sizeof(int));
  Free(valid,cg.num_node*sizeof(int));
  Free(clique,cg.num_node*sizeof(int));
  #ifdef MALLOC_FREE_CHECK
    check_alloc("graph_chromatic_number end")
  #endif
  if (c==-1) return -1; // timeout occurred
  return c;
}

int graph_chromatic_number_special(graph_t g, int lb, int ub, clock_t timeout) {
  // no cliques - assumes lb and ub for chi known.
  int c,i,j,place;//,*valid; //,*clique;
  for (i=0; i<g->nnodes; i++) g->c[i]=-1; // in case already colored
  if (g->nedges==0) { // empty graph
    for (i=0; i<g->nnodes; i++) g->c[i]=0;
    return 1;
  }
  if (g->nedges==1) {
    for (i=0; i<g->nnodes; i++) g->c[i]=0;
    for (i=0; i<g->nnodes; i++) { // find the edge
      if (g->d[i]>0) {
        g->c[i]=1;
        return 2;
      }
    }
  }
  if (2*g->nedges==g->nnodes*(g->nnodes-1)) { // complete graph
    for (i=0; i<g->nnodes; i++) g->c[i]=i;
    return g->nnodes;
  }
  #ifdef MALLOC_FREE_CHECK
    alloc=0;
  #endif
  cg.g=g;
  cg.num_node=g->nnodes;
  cg.adj=Malloc(cg.num_node*sizeof(int*));
  cg.ColorAdj=Malloc(cg.num_node*sizeof(int*));
  for (i=0; i<cg.num_node; i++) {
    cg.ColorAdj[i]=Malloc((cg.num_node+1)*sizeof(int));
    cg.adj[i]=Calloc(g->nnodes,sizeof(int));
    for (j=0; j<g->d[i]; j++) cg.adj[i][g->a[i][j]]=1;
  }
  //valid        =Malloc(cg.num_node*sizeof(int));
  //clique       =Calloc(cg.num_node,sizeof(int));
  cg.ColorClass=Malloc(cg.num_node*sizeof(int));
  cg.Order     =Malloc(cg.num_node*sizeof(int));
  cg.Handled   =Malloc(cg.num_node*sizeof(int));
  cg.ColorCount=Malloc(cg.num_node*sizeof(int));
  for (i=0; i<cg.num_node; i++)
  cg.prob_count=0;
  for (i=0; i<cg.num_node; i++) {
    for (j=0; j<=cg.num_node; j++) cg.ColorAdj[i][j]=0;
    for (j=0; j<cg.num_node; j++) if (cg.adj[i][j]) cg.ColorAdj[i][0]++;
    cg.ColorCount[i]=0;
    cg.Handled[i]=0;
    //valid[i]=1;
  }
  cg.BestColoring=ub;
  //cg.best_clique=0;
  cg.lb=lb;
  cg.num_prob=0;
  cg.max_prob=10000;
  // timing...
  cg.tps=sysconf(_SC_CLK_TCK); // clock ticks per second
  times(&cg.buffer);
  cg.start_time=cg.buffer.tms_utime;
  cg.timeout=timeout;
  if (1) {
    place=0;
    //for (i=0; i<cg.num_node; i++) if (clique[i]) {
    //  cg.Order[place]=i;
    //  cg.Handled[i]=1;
    //  place++;
    //  AssignColor(i,place);
    //  for (j=0; j<cg.num_node; j++)
    //    if (i!=j && clique[j] && !cg.adj[i][j]) {
    //      fprintf(stderr,"graph_chromatic_number internal error: result is not a clique!\n");
    //      graph_show(g);
    //      for (l=0; l<cg.num_node; l++) printf("clique[%d]=%d\n",l,clique[l]);
    //      for (l=0; l<cg.num_node; l++) {
    //        for (k=0; k<cg.num_node; k++) printf("%d ",cg.adj[l][k]);
    //        printf("\n");
    //      }
    //      exit(2);
    //    }
    //}
    c=Color(place,place);
  }
  for (i=0; i<cg.num_node; i++) {
    Free(cg.ColorAdj[i],(cg.num_node+1)*sizeof(int));
    Free(cg.adj[i],g->nnodes*sizeof(int));
  }
  Free(cg.adj,g->nnodes*sizeof(int*));
  Free(cg.ColorAdj,cg.num_node*sizeof(int*));
  Free(cg.ColorClass,cg.num_node*sizeof(int));
  Free(cg.Order,cg.num_node*sizeof(int));
  Free(cg.Handled,cg.num_node*sizeof(int));
  Free(cg.ColorCount,cg.num_node*sizeof(int));
  //Free(valid,cg.num_node*sizeof(int));
  //Free(clique,cg.num_node*sizeof(int));
  #ifdef MALLOC_FREE_CHECK
    check_alloc("graph_chromatic_number end")
  #endif
  return c;
}

int graph_edge_chromatic_number(graph_t g, clock_t timeout) {
  int chi_dash,lb,ub;
  if (g->nedges==0) return 0;
  if (g->nedges==1) return 1;
  lb=graph_max_degree(g);
  ub=lb+1;
  graph_t l=graph_line_graph(g);
  //fprintf(stderr,"before l=\n"); graph_show(l); // !!!!!!!!
  chi_dash=graph_chromatic_number_special(l,lb,ub,timeout);
  //fprintf(stderr,"after l=\n"); graph_show(l); // !!!!!!!!
  graph_clear(l);
  return chi_dash;
}

int graph_clique_number(graph_t g) {
  int i,j;
  int* u;
  clqdata.matrix_size=g->nnodes;
  u=malloc(clqdata.matrix_size*sizeof(int));
  clqdata.deg=Calloc(clqdata.matrix_size,sizeof(int));
  clqdata.IS=Calloc(clqdata.matrix_size,sizeof(int));
  clqdata.adj_list=Malloc(clqdata.matrix_size*sizeof(int*));
  clqdata.bucket=Malloc(clqdata.matrix_size*sizeof(int*));
  for (i=0; i<clqdata.matrix_size; i++) {
    clqdata.adj_list[i]=Malloc(clqdata.matrix_size*sizeof(int));
    clqdata.bucket[i]=Calloc(clqdata.matrix_size,sizeof(int));
    for (j=0; j<clqdata.matrix_size; j++) clqdata.adj_list[i][j]=1;
  }
  for (i=0; i<g->nnodes; i++) {
    bzero(u,clqdata.matrix_size*sizeof(int));
    for (j=0; j<g->d[i]; j++) u[g->a[i][j]]=1;
    for (j=0; j<clqdata.matrix_size; j++) if (!u[j] && i!=j) { // i!-j
      clqdata.adj_list[i][clqdata.deg[i]++]=j; // complement of g
    }
  }
  free(u);
  clqdata.maxsz=-1;
  clqdata.mode=0;
  clqdata.verbose=0;
  backtrack(0);
  for (i=0; i<clqdata.matrix_size; i++) {
    Free(clqdata.adj_list[i],0);
    Free(clqdata.bucket[i],0);
  }
  Free(clqdata.deg,0);
  Free(clqdata.IS,0);
  Free(clqdata.adj_list,0);
  Free(clqdata.bucket,0);
  return clqdata.maxsz;
}

//! read and process "showg -l0 -e" output
histogram_t graph_showg_reader(FILE* f, int op(graph_t), char* lbl) {
  int i,j,k,verbose=0,num_node,num_graph,m;
  graph_t g;
  histogram_t h=histogram_new(lbl);
  while (!feof(f)) { // loop over graphs in "showg -e -l0" format
    i=scanf(" Graph %d, order %d.",&num_graph,&num_node);
    if (i<=0) break; // done
    if (verbose) printf("Graph %d with %d nodes, ",num_graph,num_node);
    scanf("%d %d",&num_node,&m);
    if (verbose) printf("%d edges\n",m);
    g=graph_new(num_node);
    for (k=0; k<m; k++) {
      scanf("%d %d",&i,&j);
      if (i<num_node && j<num_node) graph_add_edge(g,i,j);
      else {
        fprintf(stderr,"graph_showg_reader: bad edge: %d--%d, n=%d\n",i,j,num_node);
        exit(1);
      }
    }
    histogram_add(h,op(g));
    graph_clear(g);
  }
  return h;
}

//! read and process "geng" output
histogram_t graph_geng_reader(FILE* f, int op(graph_t), char* lbl) {
  char *s;
  graph_t g;
  histogram_t h=histogram_new(lbl);
  while ((s=geng_getline(stdin))) {
    g=geng_stringtograph(s);
    //graph_show(g);
    histogram_add(h,op(g));
    if (s) free(s);
    graph_clear(g);
  }
  return h;
}

unsigned int xcng(void) {
  static unsigned int x=12345689,y=521289,z,w;
  z=x+y; x=y; w=(z<<8); y=w+17*(z>>24);
  if (z<x) y+=4352; if (y<w) y+=17; return y;
}

double uniform01(void) {
  return 2.3283064365386962891e-10*xcng();
}

unsigned int randint(int upto) {
  return (unsigned int)upto*((double)xcng()/(double)UINT_MAX);
  //return (int)upto*uniform01();
}

static void shuffle(int p[], int n, int m) {
  // uses a different RNG so we get the same sequence of random graphs
  int i,j,k,l,p0=p[0];
  for (k=0; k<m; k++) {
    i=randint(n);
    while ((j=randint(n))==i);
    l=p[i]; p[i]=p[j]; p[j]=l;
  }
  if (p0==p[0]) { // make sure 0 is not fixed
    while ((j=randint(n))==0);
    l=p[0]; p[0]=p[j]; p[0]=l;
  }
}

double geometric_deviate(double p) {
  // number of failures before a success in a Bernoulli(p) process
  if (p==0.0) return DBL_MAX;
  if (p==1.0) return 0.0;
  return floor(log1p(-UNIFORM_01)/log1p(-p));
}

double poisson_deviate(double xm) {
  static double sq,alxm,g,oldm=-1.0;
  double em,t,y;
  if (xm<12.0) {
    if (xm!=oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em=-1;
    t=1.0;
    do {
      em+=1.0;
      t*=UNIFORM_01;
    }
    while (t>g);
  } else {
    if (xm!=oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-lgamma(xm+1.0);
    }
    do {
        do {
          y=tan(M_PI*UNIFORM_01);
          em=sq*y+xm;
        } while (em<0.0);
        em=floor(em);
        t=0.9*(1.0+y*y)*exp(em*alxm-lgamma(em+1.0)-g);
    } while (UNIFORM_01>t);
  } 
  return em;
}

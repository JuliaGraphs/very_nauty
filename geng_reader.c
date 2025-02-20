// KMB 2006 Jul 18
// contains code from nauty22/showg.c by Brendan McKay

#include <stdio.h>
#include <stdlib.h>

#define BIAS6 63
#define MAXBYTE 126
#define SMALLN 62
#define TOPBIT6 32
#define C6MASK 63
#define TIMESWORDSIZE(w) ((w)<<5)
#define ADDELEMENT(i,j) graph_add_edge(g,i,j);
#define GRAPHROW(v) (v)

#define DYNALL(type,name,name_sz)\
  type *name=NULL;\
  size_t name_sz=0

#define DYNALLOC1(type,name,name_sz,sz,msg)\
  if ((size_t)(sz) > name_sz) {\
    if (name_sz) free(name);\
    name_sz = (sz);\
    if ((name=(type*)calloc(sz,sizeof(type))) == NULL)\
      fprintf(stderr,"%s\n",msg),exit(1);\
  }

#define DYNREALLOC(type,name,name_sz,sz,msg) {\
  if ((size_t)(sz) > name_sz) {\
    if ((name = (type*)realloc(name,(sz)*sizeof(type))) == NULL) \
      fprintf(stderr,"%s\n",msg),exit(1);\
    else\
      name_sz = (sz);\
    }\
  }

#define DYNFREE(name,name_sz)\
  if (name_sz) { free(name); name_sz = 0; }

char* geng_getline(FILE *f) {   /* read a line with error checking */
/* includes \n (if present) and \0.
                  Immediate EOF causes NULL return. */
  DYNALL(char,s,s_sz);
  int i,c;
  DYNALLOC1(char,s,s_sz,32,"getline");
  i = 0;
  while ((c = getc(f)) != EOF && c != '\n') {
    if (i == s_sz-2) DYNREALLOC(char,s,s_sz,s_sz+32,"getline");
    s[i++] = c;
  }
  if (i == 0 && c == EOF) return NULL;
  if (c == '\n') s[i++] = '\n';
  s[i] = '\0';
  return s;
}

int geng_graphsize(char *s) {
/* Get size of graph out of graph6 or sparse6 string. */
  char *p;
  int n;
  if (s[0] == ':')
    p = s + 1;
  else
    p = s;
  n = *p++ - BIAS6;
  if (n > SMALLN) {
    n = *p++ - BIAS6;
    n = (n << 6) | (*p++ - BIAS6);
    n = (n << 6) | (*p++ - BIAS6);
  }
  return n;
}

graph_t geng_stringtograph(char *s) {
  char *p;
  int gi, gj, n, i, j, k, v, x=0, nb;
  n = geng_graphsize(s);
  graph_t g=graph_new(n);
  p = s + 1 + (s[0] == ':');
  if (n > SMALLN)
    p += 3;
  if (s[0] != ':') {                /* graph6 format */
    k = 1;
    for (j = 1; j < n; ++j) {
      gj = GRAPHROW(j);
      for (i = 0; i < j; ++i) {
        if (--k == 0) {
          k = 6;
          x = *(p++) - BIAS6;
        }
        if (x & TOPBIT6) {
          gi = GRAPHROW(i);
          ADDELEMENT(gi, j);
          ADDELEMENT(gj, i);
        }
        x <<= 1;
      }
    }
  } else {                        /* sparse6 format */
    for (i = n - 1, nb = 0; i != 0; i >>= 1, ++nb) { }
    k = 1;
    v = 0;
    for (;;) {
      if (--k == 0) {
        k = 6;
        if (*p == '\n' || *p == '\0')
          break;
        else
          x = *p - BIAS6;
        ++p;
      } else
        x <<= 1;
      if (x & TOPBIT6)
        ++v;
      j = 0;
      for (i = 0; i < nb; ++i) {
        if (--k == 0) {
          k = 6;
          if (*p == '\n' || *p == '\0')
            break;
          else
            x = *p - BIAS6;
          ++p;
        } else
          x <<= 1;
        if (x & TOPBIT6)
          j = (j << 1) | 1;
        else
          j <<= 1;
      }
      if (i < nb)
        break;
      if (j > v)
        v = j;
      else if (v < n) {
        ADDELEMENT(GRAPHROW(v),j);
        ADDELEMENT(GRAPHROW(j),v);
      }
    }
  }
  return g;
}

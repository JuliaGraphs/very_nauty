#!/usr/bin/python
# KMB 2007 Jan 08 was regular_degree4.py, now with 'showg -l0 -e' output
# mcopy -o -v random_planar.py a:

# e.g.:
# ./random_planar.py 15 1
# d random_planar.png
# cat foo.showg | ~/graph-theory/showg-e2dot.py | neato -Tps | gv -

# KMB 2005 Sep 09
# Make a random 4-regular graph by dropping lines onto a plane.
# The intersetions of the lines define the nodes.
# ./random_planar.py && neato -Tpng foo.dot >| random_planar.png && display random_planar.png
# mcopy -o -v random_planar.py random_planar.png foo.dot a:

# idea: 
# 
# 1. randomly distribute n points in a region
# 2. through each draw an infinite line with angle distributed uniformly in [0,2pi)
# 3. Put nodes at every intersection of lines
# 4. Throw away the original points.
# 
# Generically you get a regular degree 4 (except at the edges) planar graph
# which is naturally embedded in the plane, but is not isomorphic to a subgraph
# of the square lattice.
# 
# The only problem is to extend it beyond the original region without altering the
# degree 4 part.

from random import random
from math import pi,tan
from sys import stderr,argv

#seed=1.0/3.0
#def random():
#  global seed
#  seed*=2.23606797749978969640917366
#  return seed-int(seed)

npoints=8
ngraphs=100
if len(argv)>1: npoints=int(argv[1])
if len(argv)>2: ngraphs=int(argv[2])

def intersections(lines):
  points=[] # x,y pairs
  where={}
  xs=[]     # x coordinates of nodes and their index
  line_contains=dict([(i,[]) for i in range(len(lines))])
  n=0
  for i in xrange(len(lines)):
    x0i,y0i,si=lines[i]
    for j in xrange(i+1,len(lines)):
      x0j,y0j,sj=lines[j]
      x=(sj*x0j-si*x0i+y0i-y0j)/(sj-si) # x coordinate of intersection
      if 0.0<=x<=1.0:
        y=si*(x-x0i)+y0i
        if 0.0<=y<=1.0: # intersection is in the unit square
          points.append((x,y))
          where[n]=(i,j) # point n lies on lines i and j
          line_contains[i].append(n) # lines i and j contain point n
          line_contains[j].append(n)
          xs.append((x,n))
          n+=1
  return points,line_contains,xs

f=open('foo.showg','w')
for ng in xrange(ngraphs):
  lines=[]
  for i in xrange(npoints):
    x0,y0=random(),random()
    a=pi*random() # angle
    s=tan(a)      # slope - random line is now y=s*(x-x0)+y0
    lines.append((x0,y0,s))
    if abs(s)<1.0:  # shallow slope
      x1,y1=0.0,s*(   -x0)+y0
      x2,y2=1.0,s*(1.0-x0)+y0
    else:
      x1,y1=x0+(   -y0)/s,0.0
      x2,y2=x0+(1.0-y0)/s,1.0
    #draw.line([plot(x1,y1),plot(x2,y2)],blue)
  points,line_contains,xs=intersections(lines)
  edges=set([])
  nodes=set([])
  for l in line_contains:
    nis=line_contains[l] # node indices on line l
    q=[xs[ni] for ni in nis]
    q.sort()
    for ni,nj in zip(q[:-1],q[1:]):
      nodes.add(ni[1])
      nodes.add(nj[1])
      e=[ni[1],nj[1]]
      e.sort()
      edges.add(tuple(e))
      #f.write('  %d--%d;\n'%(ni[1],nj[1]))
  f.write('Graph %d, order %d.\n'%(ng,len(points)))
  f.write('%d %d\n'%(len(points),len(edges)))
  for e in edges: f.write('%d %d  '%(e[0],e[1]))
  f.write('\n\n')

f.close()

if ngraphs==1:
  import Image,ImageDraw,ImageFont
  class Palette:
    def __init__(self):
      self.palette=[]
    def __call__(self,r,g,b): # map rgb tuple to colour index
      rgb=r,g,b
      try:
        return self.palette.index(rgb)
      except:
        i=len(self.palette)
        if i>=256: raise RuntimeError, "all palette entries are used"
        self.palette.append(rgb)
        return i
    def getpalette(self): # return flattened palette
      palette=[]
      for r,g,b in self.palette: palette=palette+[r,g,b]
      return palette
  rgb=Palette()
  black=rgb(0,0,0)
  red=rgb(255,0,0)
  blue=rgb(0,0,255)
  green=rgb(0,255,0)
  width,height=500,500
  font=ImageFont.load_default()
  img=Image.new('P',(width,height),rgb(255,255,255))
  draw=ImageDraw.Draw(img)
  def plot(x,y): return int(width*x),int(height-height*y)
  for n,point in enumerate(points): 
    x,y=plot(point[0],point[1])
    draw.rectangle([(x-2,y-2),(x+2,y+2)],red)
    draw.text((x+2,y),str(n),font=font,fill=black)
  for i in xrange(npoints):
    x0,y0,s=lines[i]
    if abs(s)<1.0:  # shallow slope
      x1,y1=0.0,s*(   -x0)+y0
      x2,y2=1.0,s*(1.0-x0)+y0
    else:
      x1,y1=x0+(   -y0)/s,0.0
      x2,y2=x0+(1.0-y0)/s,1.0
    draw.line([plot(x1,y1),plot(x2,y2)],blue)
  img.putpalette(rgb.getpalette())
  img.save('random_planar.png','png')
  # deduce the topology...
  #hdr='''  graph [fontcolor=black,fontsize=16];
  #  center=true;
  #  unflatten=true;
  #  /* ratio=1;
  #  mclimit=5;
  #  nslimit=500;
  #  nslimit1=500;
  #  remincross=true; */
  #  size="8,8";
  #  node [shape="box",fixedsize="true"]
  #  node [color=blue,fontcolor=red,fontsize=18,height="0.6",width="0.6"]
  #  edge [color=blue,fontsize=10];\n'''
  #f.write('graph g {\n')
  #f.write(hdr)
  #f.close()


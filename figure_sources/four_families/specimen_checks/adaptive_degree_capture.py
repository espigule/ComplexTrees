#!/usr/bin/env python3
"""Adaptive exact outer-disk degree certificate for every binary zipper type.

Each polygon edge and the corresponding true curve subarc lie in the same
convex disk which excludes zero, including the entire requested parameter
disk. Thus the polygon and true boundary map have the same winding number.
No fixed-depth polygon or floating-point decision enters verification.
"""
from fractions import Fraction as Q
import json,argparse,math
from pathlib import Path
from zipper_degree_capture import spec,exact_base,gmul,gadd,gsub,gconj,separated_intervals,source_interval,rational_upper_sqrt

# Map=(Gaussian numerator translation, coefficient, conjugation, denominator,
# number of 0 letters, number of 1 letters, parameter orientation reversal).
def extend(F,G,D,i,e):
 t,a,k,den,c0,c1,o=F;s,b,l=G
 return (gadd((t[0]*D,t[1]*D),gmul(a,gconj(s) if k else s)),gmul(a,gconj(b) if k else b),k^l,den*D,c0+(i==0),c1+(i==1),o^e)

def exact_map(f,w,D,e):
 F=((0,0),(1,0),0,1,0,0,0)
 for c in w:
  i=int(c);F=extend(F,f[i],D,i,e[i])
 return F

def endpoint(F,which):
 t,a,k,den,*_=F
 return (gadd(t,a) if which else t),den

def difference(P,Qp):
 z,d=P;w,h=Qp
 return (gsub((z[0]*h,z[1]*h),(w[0]*d,w[1]*d)),d*h)

def degree(points):
 w=0
 for (a,da),(b,db) in zip(points,points[1:]):
  cross=a[0]*b[1]-a[1]*b[0]
  if a[1]<=0<b[1] and cross>0:w+=1
  elif b[1]<=0<a[1] and cross<0:w-=1
 return w

def verify_adaptive(name,x,y,u,v,radius=Q(0),max_leaves=200000,max_depth=160):
 k,e=spec(name)
 if radius<0:raise ValueError('Radius must be nonnegative.')
 if not u or not v or set(u+v)-{'0','1'}:raise ValueError('Nonempty binary words required.')
 if not separated_intervals(name,u,v):raise ValueError('The source intervals must be strictly disjoint.')
 D,f=exact_base(name,x,y);U=exact_map(f,u,D,e);V=exact_map(f,v,D,e)
 # Boundary variables are the local s,t of gamma, so start traversal parity
 # at zero independently of the orientation of their source cylinders.
 U=U[:6]+(0,);V=V[:6]+(0,)
 r1=rational_upper_sqrt(x*x+y*y);r2=rational_upper_sqrt((1-x)**2+y*y);r=max(r1,r2)
 if r+radius>=1:raise ValueError('Parameter disk must lie in the contraction lens.')
 R=max(r2/(2*(1-r1)),r1/(2*(1-r2)))
 motion=2*(R+Q(1,2))*radius/(1-r-radius)
 powers0=[Q(1)];powers1=[Q(1)]
 def scale(c0,c1):
  while len(powers0)<=c0:powers0.append(powers0[-1]*r1)
  while len(powers1)<=c1:powers1.append(powers1[-1]*r2)
  return powers0[c0]*powers1[c1]
 visits=0;leaves=[];minimum=None;blocks=[];failed=None
 # Boundary order: (s,0), (1,t), (1-s,1), (0,1-t).
 for block,(W,anchor,sign,reverse,w0) in enumerate([(U,endpoint(V,0),1,False,u),(V,endpoint(U,1),-1,False,v),(U,endpoint(V,1),1,True,u),(V,endpoint(U,0),-1,True,v)]):
  stack=[(W,w0)];points=[];block_leaves=[]
  while stack:
   F,w=stack.pop();visits+=1;t,a,rev,den,c0,c1,orientation=F
   center=(gadd((2*t[0],2*t[1]),a),2*den)
   z,zd=difference(center,anchor)
   distance2=Q(z[0]*z[0]+z[1]*z[1],zd*zd)
   bound=R*scale(c0,c1)+motion
   gap=distance2-bound*bound
   if gap>0:
    p=difference(endpoint(F,orientation),anchor)
    points.append(p)
    leaves.append((block,w))
    block_leaves.append(w)
    minimum=gap if minimum is None else min(minimum,gap)
    if len(leaves)>max_leaves:failed='leaf budget';break
   else:
    if len(w)>=max_depth:failed='maximum word depth';break
    order=(0,1) if orientation==0 else (1,0)
    for i in order[::-1]:stack.append((extend(F,f[i],D,i,e[i]),w+str(i)))
  if failed:break
  points.append(difference(endpoint(W,1-W[6]),anchor))
  if reverse:points.reverse()
  if sign<0:points=[((-z[0],-z[1]),d) for z,d in points]
  if blocks:
   assert difference(blocks[-1],points[0])[0]==(0,0)
   blocks.extend(points[1:])
  else:blocks=points
 result={'schema':'adaptive-zipper-degree-capture-v1','family':name,'q':{'real':str(x),'imag':str(y)},'parameter_disk_radius':str(radius),
   'words':{'u':u,'v':v},'source_intervals':[[str(z) for z in source_interval(name,u)],[str(z) for z in source_interval(name,v)]],
   'arithmetic':'Exact Python integers and fractions.Fraction; upward square-root bounds proved by integer square comparisons.',
   'contraction_upper_bounds':[str(r1),str(r2)],'invariant_disk_radius_upper':str(R),'uniform_boundary_motion_upper':str(motion),
   'visited_cylinders':visits,'accepted_cylinders':len(leaves),'maximum_word_depth':max((len(w) for _,w in leaves),default=0)}
 if failed:
  result.update(status='unresolved',failure=failed);return result
 assert difference(blocks[0],blocks[-1])[0]==(0,0)
 d=degree(blocks)
 result.update(status='proved_nonembedded_parameter_disk' if d else 'unresolved_zero_degree',degree=d,
    minimum_squared_disk_clearance=str(minimum),approximate_minimum_squared_disk_clearance=float(minimum),
    boundary_cover=[{'side':b,'word':w} for b,w in leaves])
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('family');p.add_argument('real');p.add_argument('imag');p.add_argument('u');p.add_argument('v');p.add_argument('--radius',default='0');p.add_argument('--output')
 a=p.parse_args();r=verify_adaptive(a.family,Q(a.real),Q(a.imag),a.u,a.v,Q(a.radius));s=json.dumps(r,indent=2)+'\n'
 if a.output:Path(a.output).write_text(s);print(json.dumps({k:v for k,v in r.items() if k not in ('boundary_cover','minimum_squared_disk_clearance')},indent=2))
 else:print(s)

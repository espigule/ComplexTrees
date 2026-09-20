#!/usr/bin/env python3
"""Adaptive endpoint-chord displays of the recovered exact historical probes.

Stopping uses exact rational upper bounds on the similarity-product norm,
not estimated floating-point geometry. Coordinates themselves are displays.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from functools import lru_cache
import numpy as np
import math, json, hashlib, io, os, tempfile

ROOT=Path(__file__).resolve().parent
X=Q('0.34090875');Y=Q('0.43484625');EPS=Q(1,2000)

def upper_sqrt(x,den=10**12):
    n=math.isqrt(x.numerator*den*den//x.denominator)
    if Q(n*n,den*den)<x:n+=1
    assert Q(n*n,den*den)>=x
    return Q(n,den)

R1=upper_sqrt(X*X+Y*Y);R2=upper_sqrt((1-X)**2+Y*Y)
DELTA=upper_sqrt((X-Q(1,2))**2+Y*Y)
R=max(R1,R2);GLOBAL=DELTA/(1-R)
def sharpen_global_bound(name,depth=14):
    from zipper_degree_capture import exact_curve
    points,den,_,_=exact_curve(name,X,Y,depth)
    steps=2**depth
    maximum=max((px*steps-j*den)**2+(py*steps)**2 for j,(px,py) in enumerate(points))
    polygon_deviation=upper_sqrt(Q(maximum,(den*steps)**2))
    remainder=R**depth*GLOBAL
    result=polygon_deviation+remainder
    assert result<GLOBAL
    return result,{'reference_polygon_depth':depth,'polygon_identity_deviation_upper':str(polygon_deviation),
        'reference_curve_tail_upper':str(remainder),'sharpened_global_deviation_upper':str(result),
        'sharpened_global_deviation_decimal':float(result),
        'proof':'The convex norm of P_N(t)-t attains a maximum at a polygon vertex on each source interval. Add r^N delta/(1-r) to its exact vertex norm bound.'}

def build(name):
    sharp,sharp_record=sharpen_global_bound(name)
    @lru_cache(None)
    def bound(m,n):return R1**m*R2**n*sharp
    e=(int(name[2]),int(name[3]));ks=(0,int(name[1]=='O'))
    q=complex(float(X),float(Y));b=1-q
    fs=[(q if e[0] else 0j,-q if e[0] else q,ks[0]),(1+0j if e[1] else q,-b if e[1] else b,ks[1])]
    # t, a, spatial conjugation, source orientation, zero count, one count.
    stack=[(0j,1+0j,0,0,0,0)];points=[0j];counts=Counter();gaps=[];maxgap=0.;maxbound=Q(0)
    while stack:
        t,a,k,reverse,m,n=stack.pop()
        if bound(m,n)<=EPS:
            start=t+a*(1 if reverse else 0)
            finish=t+a*(0 if reverse else 1)
            maxgap=max(maxgap,abs(start-points[-1]));points.append(finish)
            counts[m,n]+=1;maxbound=max(maxbound,bound(m,n))
        else:
            assert m+n<100
            order=(0,1) if reverse==0 else (1,0)
            for i in reversed(order):
                ft,fa,fk=fs[i]
                stack.append((t+a*(ft.conjugate() if k else ft),a*(fa.conjugate() if k else fa),k^fk,reverse^e[i],m+(i==0),n+(i==1)))
    assert abs(points[0])<1e-14 and abs(points[-1]-1)<1e-14
    assert maxgap<1e-12
    edges=len(points)-1
    # Exact source-interval lengths partition [0,1].
    partition=sum((Q(cnt,2**(m+n)) for (m,n),cnt in counts.items()),Q(0))
    assert partition==1
    record={'family':name,'probe':'C','q_exact':[str(X),str(Y)],'requested_uniform_error':str(EPS),
      'certified_upper_bound':str(maxbound),'uniform_error_upper_decimal':float(maxbound),
      'r1_upper':str(R1),'r2_upper':str(R2),'delta_upper':str(DELTA),'sharp_global_bound':sharp_record,
      'vertices':len(points),'segments':edges,'max_source_depth':max(m+n for m,n in counts),
      'min_source_depth':min(m+n for m,n in counts),'maximum_display_endpoint_join_error':maxgap,
      'exact_source_partition_length':str(partition),
      'leaf_counts':[{'m':m,'n':n,'count':count} for (m,n),count in sorted(counts.items())],
      'rule':'For each source cylinder w, stop only if R1^count0(w) R2^count1(w) sharpened_global_deviation_upper <= 1/2000. Endpoint order follows cumulative XOR of e_i, independently of spatial conjugation.',
      'display_arithmetic':'Float64 complex coordinates; exact rational stopping and source partition checks.'}
    arr=np.array(points,dtype=np.complex128)
    # First root branch ends at source t=1/2: exact lengths locate split.
    total=Q(0);idx=0
    # This is the image q; choose closest numerically, retaining all its possible secondary visits.
    candidates=np.flatnonzero(abs(arr-q)<1e-13)
    # Root partition appears midway through traversal's dyadic source length, but count halves vary by endpoint order.
    # Re-run source ordering counts cheaply to find exact t=1/2 index.
    cut=0
    for i in (0,):
        # count all leaves under first root branch, independent of endpoint orientation
        s=[(1,0)]
        while s:
            m,n=s.pop()
            if bound(m,n)<=EPS:cut+=1
            else:s.extend([(m+1,n),(m,n+1)])
    assert abs(arr[cut]-q)<1e-12
    record['first_level_join_index']=cut
    data=io.BytesIO();np.savez_compressed(data,curve=arr,first_level_join_index=np.array(cut))
    target=ROOT/f'{name}_C_adaptive_curve.npz'
    tmp=target.with_suffix('.npz.tmp');tmp.write_bytes(data.getvalue());os.replace(tmp,target)
    record['array_file']=target.name;record['array_sha256']=hashlib.sha256(target.read_bytes()).hexdigest()
    print(name,len(points),maxgap,float(maxbound),flush=True)
    return record

if __name__=='__main__':
    records=[build(name) for name in ['DD10','DD11','DO10','DO11']]
    (ROOT/'adaptive_curve_metadata.json').write_text(json.dumps({'schema':'adaptive-zipper-display-v1','builder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'results':records},indent=2)+'\n')

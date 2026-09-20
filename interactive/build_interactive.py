#!/usr/bin/env python3
"""Build the offline four-family atlas from certified q-cells.

Only all-positive 2 x 2 blocks are retained. Conjugation reflects them to the
lower q half-plane. The smooth map of each rectangle is approximated by two
display triangles; neither meshing nor curve drawing supplies a certificate.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
import base64, hashlib, json, math, os
import numpy as np

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'figure_sources' / 'zippers'
FAMILIES = ('DD10', 'DD11', 'DO10', 'DO11')
PROBES = {'A': ('0.5','0.2'), 'B': ('0.5','0.4'), 'C': ('0.34090875','0.43484625')}

def sqrt_upper(x, den=10**12):
    n = math.isqrt(x.numerator*den*den//x.denominator)
    if F(n*n,den*den)<x: n+=1
    return F(n,den)

def curve(family, key):
    x,y = map(F,PROBES[key]); q=complex(float(x),float(y)); e=(1,int(family[-1]))
    r1=sqrt_upper(x*x+y*y); r2=sqrt_upper((1-x)**2+y*y)
    delta=sqrt_upper((x-F(1,2))**2+y*y); glob=delta/(1-max(r1,r2))
    epsilon=F(1,200)
    @lru_cache(None)
    def bound(m,n):return r1**m*r2**n*glob
    fs=[(q,-q,0),(1+0j if e[1] else q,(-1)**e[1]*(1-q),int(family[1]=='O'))]
    stack=[(0j,1+0j,0,0,0,0)]; pts=[0j]; counts=Counter(); maxbound=F(0); join=0
    while stack:
        t,a,k,reverse,m,n=stack.pop()
        if bound(m,n)<=epsilon:
            start=t+a*reverse; end=t+a*(1-reverse)
            assert abs(start-pts[-1])<1e-12
            pts.append(end); counts[m,n]+=1; maxbound=max(maxbound,bound(m,n))
        else:
            order=(0,1) if not reverse else (1,0)
            for i in reversed(order):
                ft,fa,fk=fs[i]
                stack.append((t+a*(ft.conjugate() if k else ft),a*(fa.conjugate() if k else fa),k^fk,reverse^e[i],m+(i==0),n+(i==1)))
    stack=[(1,0)]
    while stack:
        m,n=stack.pop()
        if bound(m,n)<=epsilon:join+=1
        else:stack.extend([(m+1,n),(m,n+1)])
    assert abs(pts[join]-q)<1e-12
    assert sum((F(n,2**(m+k)) for (m,k),n in counts.items()),F(0))==1
    values=np.column_stack((np.real(pts),np.imag(pts))).astype('<f4')
    rounding=float(np.max(np.abs(values.astype(float)-np.column_stack((np.real(pts),np.imag(pts))))))
    stable=key!='C' or family in ('DD11','DO11')
    status='Certified stable' if stable else ('Self-intersection proved' if family=='DO10' else 'Curve example')
    return {'array':base64.b64encode(values.tobytes()).decode(),'count':len(pts),'join':join,
        'q':[float(x),float(y)],'q_exact':[str(x),str(y)],'stable':stable,'status':status,
        'display_error_bound':float(maxbound)+math.sqrt(2)*rounding,
        'source_partition':'1','max_depth':max(m+n for m,n in counts)}, values

def main():
    data={'families':{},'probes':PROBES,'schema':'four-family-offline-atlas-v1',
          'colormap':['#253A78','#255BA1','#338CA7','#77B7AE','#E4C55E']}
    qa={'schema':'four-family-offline-atlas-build-qa-v1','only_selected_families':list(FAMILIES),
        'projection':'((4+rho*cos(arg(a)))*cos(arg(b)), (4+rho*cos(arg(a)))*sin(arg(b)), rho*sin(arg(a)))',
        'fourth_coordinate':'w = |a|/(|a|+|b|), encoded by color and weight filter',
        'positive_block_rule':'All four constituent q-cells must have a positive certificate.',
        'cell_edge':1/256,'fourth_coordinate_filter':'Displayed faces have centre-weight within +/- 0.025; this is a finite band, not an exact slice.',
        'families':{},'curves':{}}
    for family in FAMILIES:
        source=SOURCE/(family+'_mask.npy'); mask=np.load(source)
        positive=mask.reshape(128,2,256,2).all(axis=(1,3))
        cells=int(positive.sum()); qadata={'source':str(source.relative_to(ROOT.parent)),
            'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_positive_cells':int(mask.sum()),
            'retained_upper_rectangles':cells,'displayed_rectangles_with_conjugates':2*cells}
        data['families'][family]={'mask':base64.b64encode(np.packbits(positive,bitorder='little').tobytes()).decode(),'count':2*cells,'curves':{}}
        qa['families'][family]=qadata
        for key in PROBES:
            record,coords=curve(family,key)
            data['families'][family]['curves'][key]=record
            qa['curves'][family+'-'+key]={k:v for k,v in record.items() if k!='array'}
            assert record['display_error_bound']<0.0051
    island=np.genfromtxt(SOURCE/'DD11_island_zoom.csv',delimiter=',',names=True)
    assert len(island)==4096 and np.all(island['status']==1)
    data['island']={'bounds':[float(island['x0'].min()),float(island['x1'].max()),float(island['y0'].min()),float(island['y1'].max())],
                    'cells':4096,'sha256':hashlib.sha256((SOURCE/'DD11_island_zoom.csv').read_bytes()).hexdigest()}
    template=(ROOT/'atlas_template.html').read_text()
    assert template.count('__ATLAS_DATA__')==1
    html=template.replace('__ATLAS_DATA__',json.dumps(data,separators=(',',':')))
    path=ROOT/'Four_Family_Stable_Atlas.html'; tmp=path.with_suffix('.html.tmp')
    tmp.write_text(html);os.replace(tmp,path)
    qa['html_bytes']=path.stat().st_size
    qa['html_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    qa['offline']=True; qa['external_dependencies']=[]
    (ROOT/'INTERACTIVE_BUILD_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
    print(json.dumps({'file':str(path),'bytes':qa['html_bytes'],'rectangles':sum(a['count'] for a in data['families'].values()),'curve_count':12}))

if __name__=='__main__':main()

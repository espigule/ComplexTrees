#!/usr/bin/env python3
"""Recompute the adaptive C display's global reference bound independently.

Common-denominator integer iteration is written here independently of the
producer's exact_curve; only saved records and point arrays are consumed.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'specimen_checks'
QX,QY,D=272727,347877,800000


def independently_iterate(name,depth):
    pts=[(0,0),(1,0)]; den=1
    e2=int(name[-1]); opposite=name[:2]=='DO'
    for _ in range(depth):
        one=[];two=[]
        for px,py in reversed(pts):
            # F1(z)=q-qz and source orientation is reversed.
            one.append((QX*den-QX*px+QY*py,QY*den-QX*py-QY*px))
        iterator=reversed(pts) if e2 else pts
        bx,by=(QX-D,QY) if e2 else (D-QX,-QY)
        tx,ty=(D,0) if e2 else (QX,QY)
        for px,py in iterator:
            if opposite:py=-py
            two.append((tx*den+bx*px-by*py,ty*den+bx*py+by*px))
        assert one[-1]==two[0]
        pts=one+two[1:];den*=D
    assert pts[0]==(0,0) and pts[-1]==(den,0)
    return pts,den


def main():
    records=json.loads((SOURCE/'adaptive_curve_metadata.json').read_text())
    output=[]
    for rec in records['results']:
        name=rec['family'];sharp=rec['sharp_global_bound'];depth=sharp['reference_polygon_depth']
        pts,den=independently_iterate(name,depth);steps=2**depth
        maxdev2=max(F((px*steps-j*den)**2+(py*steps)**2,(den*steps)**2) for j,(px,py) in enumerate(pts))
        upper=F(sharp['polygon_identity_deviation_upper'])
        assert upper*upper>=maxdev2
        r1,r2=F(rec['r1_upper']),F(rec['r2_upper']);delta=F(rec['delta_upper']);r=max(r1,r2)
        assert r1*r1>=F(QX*QX+QY*QY,D*D)
        assert r2*r2>=F((D-QX)**2+QY*QY,D*D)
        assert delta*delta>=F((QX-D//2)**2+QY*QY,D*D)
        tail=r**depth*delta/(1-r)
        assert tail==F(sharp['reference_curve_tail_upper'])
        global_upper=upper+tail
        assert global_upper==F(sharp['sharpened_global_deviation_upper'])
        maxleaf=max(r1**leaf['m']*r2**leaf['n']*global_upper for leaf in rec['leaf_counts'])
        assert maxleaf==F(rec['certified_upper_bound'])
        assert maxleaf<F(1,2000)
        assert sum((F(leaf['count'],2**(leaf['m']+leaf['n'])) for leaf in rec['leaf_counts']),F(0))==1
        path=SOURCE/rec['array_file']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['array_sha256']
        arr=np.load(path);z=arr['curve'];join=int(arr['first_level_join_index'])
        assert len(z)==rec['vertices'] and join==rec['first_level_join_index']
        assert z[0]==0 and abs(z[-1]-1)<1e-12 and abs(z[join]-complex(QX/D,QY/D))<1e-12
        assert np.all(np.isfinite(z))
        output.append({'family':name,'independent_integer_reference_depth':depth,'reference_vertices':len(pts),
            'reference_maximum_deviation_squared':str(maxdev2),'strict_uniform_chord_error_upper':str(maxleaf),
            'strict_uniform_chord_error_upper_decimal':float(maxleaf),'segments':len(z)-1,
            'first_level_join_index':join,'array_sha256':rec['array_sha256'],'passed':True})
    result={'passed':True,'method':'Independent common-denominator integer polygon iteration and Fraction bounds; producer geometry modules not imported.',
            'scope':'Exact error bounds apply to mathematical endpoint chords; float64 display vertices add unrounded numerical evaluation error, not an interval rendering guarantee.',
            'records':output}
    (ROOT/'independent_curve_bound_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'families':[{k:r[k] for k in ('family','segments','strict_uniform_chord_error_upper_decimal','first_level_join_index')} for r in output]},indent=2))


if __name__=='__main__':main()

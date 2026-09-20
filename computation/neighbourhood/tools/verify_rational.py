#!/usr/bin/env python3
"""Exact one-sided disconnection test for rational semilinear maps and a neighbourhood.

Standard library only. Input multipliers must be decimal/fraction STRINGS,
not floats or angles. Success supplies a finite cover and a closed product of
complex parameter disks. An angle box needs a separate containment argument.
"""
from fractions import Fraction as Q
from math import isqrt
import argparse
import hashlib
import json

VERSION = '2.0.0'

def parse_q(s):
    if not isinstance(s, str) or len(s) > 1200:
        raise ValueError('Coefficients must be exact decimal or fraction strings.')
    q = Q(s)
    if max(q.numerator.bit_length(), q.denominator.bit_length()) > 4096:
        raise ValueError('Coefficient exceeds 4096-bit input bound.')
    return q

def multiply(a,b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def conjugate(a,e):
    return (a[0], -a[1] if e else a[1])

def compose(a,b):
    # a o b; tuple is (translation, multiplier, conjugation parity, radius bound).
    t,m,e,r=a; u,n,f,s=b
    mu=multiply(m,conjugate(u,e))
    return ((t[0]+mu[0],t[1]+mu[1]), multiply(m,conjugate(n,e)), e^f,r*s)

def norm2(a):
    return a[0]*a[0]+a[1]*a[1]

def sqrt_upper(x,bits=80):
    scale=1<<bits
    k=isqrt((x.numerator*scale*scale)//x.denominator)
    if k*k*x.denominator < x.numerator*scale*scale:
        k+=1
    return Q(k,scale)

def verify(spec,max_depth=14,max_pairs=100000):
    if not 1 <= max_depth <= 30 or not 1 <= max_pairs <= 1000000:
        raise ValueError('Use depth 1..30 and max_pairs 1..1000000.')
    coefficients=spec.get('multipliers')
    parity=spec.get('orientation')
    if not isinstance(coefficients,list) or len(coefficients)!=2 or not isinstance(parity,list) or len(parity)!=2 or any(type(x) is not int or x not in (0,1) for x in parity):
        raise ValueError('Two multipliers and two integer orientation bits required.')
    values=[]
    for z in coefficients:
        if not isinstance(z,list) or len(z)!=2: raise ValueError('Each multiplier is [real, imaginary].')
        values.append(tuple(map(parse_q,z)))
    if any(not 0 < norm2(z) < 1 for z in values):
        raise ValueError('Nonzero strict contractions required.')
    bounds=[sqrt_upper(norm2(z)) for z in values]
    if max(bounds)>=1: raise ValueError('80-bit enclosure cannot resolve contraction; increase precision in reviewed code.')
    radius=1/(1-max(bounds))
    branches=[((Q(-1 if i==0 else 1),Q(0)),values[i],parity[i],bounds[i]) for i in range(2)]
    input_record={'multipliers':[[str(x) for x in z] for z in values],'orientation':parity}
    report={'schema':'rational-disconnection/2','verifier_version':VERSION,'input':input_record,
      'input_sha256':hashlib.sha256(json.dumps(input_record,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
      'arithmetic':'exact rational; rational outward norm bound','claim_scope':'rational centre and, on success, the explicit product of complex parameter disks',
      'status':'unresolved','connectedness_certified':False,'disconnection_certified':False,
      'max_depth':max_depth,'max_pairs':max_pairs,'enclosing_radius':str(radius),
      'norm_upper_bounds':[str(q) for q in bounds],'tested_pairs':0,'pruned_pairs':0,'depth_reached':0}
    current=[(branches[0],branches[1],'-','+')]
    minimum=None
    gap_lower=None
    leaves=[]
    for depth in range(1,max_depth+1):
        survivors=[]
        report['depth_reached']=depth
        for a,b,u,v in current:
            if report['tested_pairs']>=max_pairs:
                report['reason']='pair budget exhausted'; return report
            report['tested_pairs']+=1
            d=(a[0][0]-b[0][0],a[0][1]-b[0][1])
            margin=norm2(d)-(radius*(a[3]+b[3]))**2
            if margin>0:
                report['pruned_pairs']+=1
                minimum=margin if minimum is None else min(minimum,margin)
                # Rational lower bound for |centre difference| - radius sum.
                gap=margin/(sqrt_upper(norm2(d))+radius*(a[3]+b[3]))
                gap_lower=gap if gap_lower is None else min(gap_lower,gap)
                leaves.append([u,v])
            else:
                survivors.append((a,b,u,v))
        if not survivors:
            q=max(bounds)
            qbar=(1+q)/2
            nonzero=min(max(abs(z[0]),abs(z[1])) for z in values)/2
            epsilon=min((1-q)/2,gap_lower*(1-qbar)**2/4,nonzero)
            report.update(status='certified-disconnected',disconnection_certified=True,
              reason='all first-level cylinder pairs strictly separated',minimum_squared_margin=str(minimum),
              contact_gap_lower_bound=str(gap_lower),
              separation_cover={'alphabet':'-+','composition':'outermost first','leaves':leaves},
              parameter_neighbourhood={'metric':'max(|lambda_minus-mu_minus|, |lambda_plus-mu_plus|)',
                'closed_radius':str(epsilon),'contraction_upper_bound':str(qbar),
                'contact_gap_lower_bound':str(gap_lower/2),'orientation':parity,
                'claim':'Every parameter in this closed product of complex disks is disconnected.'})
            return report
        if depth==max_depth:
            report.update(reason='enclosure pairs survive at requested depth',surviving_pairs=len(survivors)); return report
        # Budget caps growth before allocating a potentially exponential frontier.
        if len(survivors)*4 > max_pairs-report['tested_pairs']:
            report.update(reason='pair budget cannot cover next frontier',surviving_pairs=len(survivors)); return report
        current=[(compose(a,i),compose(b,j),u+'-+'[ii],v+'-+'[jj])
                 for a,b,u,v in survivors for ii,i in enumerate(branches) for jj,j in enumerate(branches)]
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input');p.add_argument('--depth',type=int,default=14);p.add_argument('--max-pairs',type=int,default=100000);p.add_argument('--output')
    a=p.parse_args()
    with open(a.input) as f: spec=json.load(f)
    try: result=verify(spec,a.depth,a.max_pairs)
    except (ValueError,KeyError,TypeError,ZeroDivisionError) as e: p.error(str(e))
    out=json.dumps(result,indent=2)+'\n'
    if a.output:
        with open(a.output,'w') as f:f.write(out)
    else: print(out,end='')

if __name__=='__main__':main()

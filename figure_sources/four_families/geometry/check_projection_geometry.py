#!/usr/bin/env python3
"""Independent exact-rational and numerical audit of the display chart.

Exact rational complex arithmetic and a real 2x2 linear solve verify the
normalization formulas; binary64 tests check coordinate reconstruction and the
global periodic projection. This does not rerun the stability certificates.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
import json
import numpy as np
from projection_geometry import FAMILIES, canonical_coefficients, atlas_coordinates, coefficient_torus, inverse_coefficient_torus


@dataclass(frozen=True)
class G:
    x: Q = Q(0)
    y: Q = Q(0)
    def __post_init__(self):
        object.__setattr__(self, 'x', Q(self.x))
        object.__setattr__(self, 'y', Q(self.y))
    @staticmethod
    def cast(z): return z if isinstance(z,G) else G(Q(z))
    def __add__(self,z):
        z=G.cast(z); return G(self.x+z.x,self.y+z.y)
    __radd__=__add__
    def __neg__(self): return G(-self.x,-self.y)
    def __sub__(self,z): return self+-G.cast(z)
    def __rsub__(self,z): return G.cast(z)+-self
    def __mul__(self,z):
        z=G.cast(z); return G(self.x*z.x-self.y*z.y,self.x*z.y+self.y*z.x)
    __rmul__=__mul__
    def __truediv__(self,z):
        z=G.cast(z); r=z.norm(); return self*z.conj()*G(1/r)
    def conj(self): return G(self.x,-self.y)
    def norm(self): return self.x*self.x+self.y*self.y
    def complex(self): return complex(float(self.x),float(self.y))


def real_matrix(a, opposite):
    return ((a.x,a.y), (a.y,-a.x)) if opposite else ((a.x,-a.y),(a.y,a.x))


def independent_normalize(q,family):
    e2=int(family[-1]); opposite=family.startswith('DO')
    t1,a1=q,-q; t2,a2=q+e2*(1-q),(-1)**e2*(1-q)
    L1=real_matrix(a1,False); L2=real_matrix(a2,opposite)
    A=[[Q(int(i==j))-(L1[i][j]+L2[i][j])/2 for j in range(2)] for i in range(2)]
    T=(t1+t2)*G(Q(1,2))
    det=A[0][0]*A[1][1]-A[0][1]*A[1][0]
    assert det>0
    m=G((T.x*A[1][1]-A[0][1]*T.y)/det,(A[0][0]*T.y-T.x*A[1][0])/det)
    d=t2+a2*(m.conj() if opposite else m)-m
    assert d.norm()>0
    a=a1; b=a2*d.conj()/d if opposite else a2
    assert (t1+a1*m-m)/d==G(-1)
    assert (t2+a2*(m.conj() if opposite else m)-m)/d==G(1)
    # Check the whole conjugacy at three off-axis rational probes.
    for z in (G(Q(2,3),Q(-4,7)),G(Q(-11,5),Q(3,8)),G(0)):
        original=m+d*z
        assert (t1+a1*original-m)/d==-1+a*z
        assert (t2+a2*(original.conj() if opposite else original)-m)/d==1+b*(z.conj() if opposite else z)
    return a,b,m,d


def main():
    exact_rows=[]; points=[]; maximum_float_error=0.
    for i in range(1,40):
        for j in range(-34,35):
            q=G(Q(i,40),Q(j,40))
            if q.norm()<1 and (1-q).norm()<1: points.append(q)
    for family in FAMILIES:
        count=0
        for q in points:
            a,b,m,d=independent_normalize(q,family)
            expect_b=(-1)**int(family[-1])*(1-q)
            if family.startswith('DO'):
                s=q.norm(); D=3*(1+2*q.x)
                N=G(3*q.x-2*q.y*q.y,q.y*(2*q.x-1)) if family=='DO10' else 1+q+s*(2*q-1)
                expected_m=(4*q+2*q.conj())/D if family=='DO10' else (1+3*q+2*s)/D
                assert m==expected_m and d==N/D
                expect_b=expect_b*N.conj()/N
            assert a==-q and b==expect_b
            assert a.norm()==q.norm() and b.norm()==(1-q).norm()
            af,bf=canonical_coefficients(q.complex(),family)
            maximum_float_error=max(maximum_float_error,abs(complex(af)-a.complex()),abs(complex(bf)-b.complex()))
            count+=1
        exact_rows.append({'family':family,'exact_rational_parameters':count,'normalization_and_conjugacy':'passed'})
    qs=np.asarray([q.complex() for q in points])
    phase_error=0.; inverse_error=0.; seam_error=0.; radial_error=0.
    for family in FAMILIES:
        a,b=canonical_coefficients(qs,family)
        d=atlas_coordinates(qs,family)
        aa=2*d['w']/d['rho']*np.exp(1j*(d['alpha']-d['gamma']))
        bb=2*(1-d['w'])/d['rho']*np.exp(-1j*(d['alpha']+d['gamma']))
        phase_error=max(phase_error,float(np.max(abs(a-aa))),float(np.max(abs(b-bb))))
        xyz=coefficient_torus(d['theta1'],d['theta2'],d['rho'])
        inv=inverse_coefficient_torus(xyz)
        inverse_error=max(inverse_error,float(np.max(abs(np.exp(1j*d['theta1'])-np.exp(1j*inv['theta1'])))),float(np.max(abs(np.exp(1j*d['theta2'])-np.exp(1j*inv['theta2'])))),float(np.max(abs(d['rho']-inv['rho']))))
        for k,l in ((1,1),(1,-1),(2,0),(-1,1)):
            alpha=d['alpha']+k*np.pi; gamma=d['gamma']+l*np.pi
            moved=coefficient_torus(alpha-gamma,-alpha-gamma,d['rho'])
            seam_error=max(seam_error,float(np.max(abs(xyz-moved))))
        c1,c2=1/a,1/b
        ww=abs(c2)/(abs(c1)+abs(c2)); rr=2*abs(c1)*abs(c2)/(abs(c1)+abs(c2))
        radial_error=max(radial_error,float(np.max(abs(ww-d['w']))),float(np.max(abs(rr-d['rho']))))
        assert np.all(abs(c1)>1) and np.all(abs(c2)>1)
    assert maximum_float_error<1e-12 and phase_error<1e-12 and inverse_error<1e-12 and seam_error<1e-12 and radial_error<1e-12
    report={'status':'passed','scope':'Exact rational chart and conjugacy checks plus numerical display checks; not stability certification.',
        'exact_checks':exact_rows,'rational_parameters_per_family':len(points),'total_exact_normalizations':4*len(points),
        'binary64_maximum_errors':{'coefficients':maximum_float_error,'atlas_roundtrip':phase_error,'torus_inverse':inverse_error,'phase_lattice_invariance':seam_error,'exterior_radial_coordinates':radial_error},
        'projection':'X=(4+rho*cos(theta1))*cos(theta2); Y=(4+rho*cos(theta1))*sin(theta2); Z=rho*sin(theta1)',
        'fourth_coordinate':'w, separately encoded; projection is injective after fixing w and orientation type on the displayed connectedness band 0<rho<=2 (not the unbounded full atlas).'}
    out=Path(__file__).with_name('projection_geometry_checks.json');out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()

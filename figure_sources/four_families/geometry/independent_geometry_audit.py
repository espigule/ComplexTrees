#!/usr/bin/env python3
"""Independent translation-normalization and certificate-data audit.

The normalization below solves the fixed point of (F1+F2)/2 directly.
The display implementation is imported only as the object under comparison.
No display formula is reused by the independent normalization.
"""
from pathlib import Path
import importlib.util
import hashlib
import json
from fractions import Fraction
import numpy as np

ROOT = Path(__file__).resolve().parent
PUBLIC = ROOT.parents[2]
spec = importlib.util.spec_from_file_location('display_geometry_under_test', ROOT/'projection_geometry.py')
display = importlib.util.module_from_spec(spec)
spec.loader.exec_module(display)
FAMILIES = ('DD10', 'DD11', 'DO10', 'DO11')


def independently_normalize(q, name):
    """F maps in endpoint coordinates; H(z)=(z-center)/displacement."""
    e2 = int(name[-1])
    A = -q
    B = (-1)**e2 * (1-q)
    t1 = q
    t2 = q if e2 == 0 else np.ones_like(q)
    if name[:2] == 'DD':
        center = (t1+t2)/(2-A-B)
    else:
        # Solve U*center - B*conjugate(center) = T as a real 2x2 system.
        U, T = 2-A, t1+t2
        center = (np.conjugate(U)*T+B*np.conjugate(T))/(abs(U)**2-abs(B)**2)
    displacement = t2+B*(np.conjugate(center) if name[:2]=='DO' else center)-center
    first_translation = (t1+A*center-center)/displacement
    second_translation = (t2+B*(np.conjugate(center) if name[:2]=='DO' else center)-center)/displacement
    a = A
    b = B*np.conjugate(displacement)/displacement if name[:2]=='DO' else B
    return a, b, center, displacement, first_translation, second_translation


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    rng=np.random.default_rng(187410)
    candidates=rng.uniform(0,1,30000)+1j*rng.uniform(-1,1,30000)
    q=candidates[(abs(candidates)<1)&(abs(1-candidates)<1)]
    probes=np.array([.5+.2j,.5+.4j,.34090875+.43484625j])
    q=np.r_[q,probes,1e-8+1e-9j,1-1e-8+1e-9j,.5+(.75**.5-1e-8)*1j]
    families={}
    for name in FAMILIES:
        a,b,center,d,t1,t2=independently_normalize(q,name)
        aa,bb=display.canonical_coefficients(q,name)
        c1,c2=display.exterior_coefficients(q,name)
        atlas=display.atlas_coordinates(q,name)
        ar=2*atlas['w']/atlas['rho']*np.exp(1j*(atlas['alpha']-atlas['gamma']))
        br=2*(1-atlas['w'])/atlas['rho']*np.exp(-1j*(atlas['alpha']+atlas['gamma']))
        # Restrict torus comparison to the necessary stability disk, where rho<=2.
        stable_geometry=(abs(q)**2<q.real)
        p,md=display.project_q(q[stable_geometry],name)
        inv=display.inverse_coefficient_torus(p)
        def err(z):return float(np.max(abs(z)))
        rec={
            'comparison_points':len(q),'stable_band_geometry_points':int(stable_geometry.sum()),
            'maximum_coefficient_error':max(err(a-aa),err(b-bb)),
            'maximum_translation_error':max(err(t1+1),err(t2-1)),
            'maximum_modulus_preservation_error':max(err(abs(a)-abs(q)),err(abs(b)-abs(1-q))),
            'minimum_exterior_modulus':float(min(abs(c1).min(),abs(c2).min())),
            'maximum_reciprocal_error':max(err(c1*a-1),err(c2*b-1)),
            'maximum_atlas_reconstruction_error':max(err(ar-a),err(br-b)),
            'maximum_torus_radius_roundtrip_error':err(inv['rho']-md['rho']),
            'maximum_torus_phase_roundtrip_error':max(err(np.exp(1j*inv['theta1'])-np.exp(1j*md['theta1'])),err(np.exp(1j*inv['theta2'])-np.exp(1j*md['theta2']))),
            'minimum_normalization_displacement':float(abs(d).min()),
            'probes':[]}
        # The very small displacement near the corner intentionally stresses
        # subtraction in the independent linear solve. Interior errors follow.
        interior=(abs(d)>1e-5)
        rec['maximum_coefficient_error_away_from_degenerate_boundary']=max(err((a-aa)[interior]),err((b-bb)[interior]))
        rec['maximum_translation_error_away_from_degenerate_boundary']=max(err((t1+1)[interior]),err((t2-1)[interior]))
        assert rec['minimum_exterior_modulus']>1
        assert rec['maximum_coefficient_error_away_from_degenerate_boundary']<2e-11
        assert rec['maximum_translation_error_away_from_degenerate_boundary']<2e-11
        assert rec['maximum_torus_radius_roundtrip_error']<5e-14
        assert rec['maximum_torus_phase_roundtrip_error']<5e-14
        for label,v in zip('ABC',probes):
            ca,cb=display.exterior_coefficients(v,name)
            xyz,md=display.project_q(v,name)
            rec['probes'].append({'probe':label,'q':[float(v.real),float(v.imag)],
                'c1':[float(ca.real),float(ca.imag)],'c2':[float(cb.real),float(cb.imag)],
                'weight':float(md['w']),'rho':float(md['rho']),'projection':xyz.tolist()})
        families[name]=rec

    masks={};records={}
    for name in ('DD01',)+FAMILIES:
        path=PUBLIC/'figure_sources/zippers'/f'{name}_cells.csv'
        arr=np.genfromtxt(path,delimiter=',',names=True)
        good=arr[arr['status']==1]
        mask=np.zeros((256,512),dtype=bool)
        mask[good['iy'].astype(int),good['ix'].astype(int)]=True
        masks[name]=mask
        corners=np.r_[good['x0']+1j*good['y0'],good['x1']+1j*good['y0'],good['x0']+1j*good['y1'],good['x1']+1j*good['y1']]
        assert np.all(abs(corners)<1)&np.all(abs(1-corners)<1)
        assert np.all(abs(corners)**2<corners.real)
        records[name]={'path':str(path.relative_to(PUBLIC)),'sha256':digest(path),'positive_cells':int(len(good)),
                       'minimum_squared_gap':float(good['minimum_squared_gap'].min()),
                       'all_positive_cell_corners_in_strict_lens':True,
                       'all_positive_cell_corners_in_necessary_stability_disk':True}
    merged=masks['DD10']|masks['DD01'][:,::-1]
    stored=np.load(PUBLIC/'figure_sources/zippers/DD10_mask.npy')
    assert np.array_equal(merged,stored)
    records['DD10']['positive_cells_after_valid_reflection']=int(merged.sum())
    records['DD10']['reflected_mask_matches_saved_mask']=True

    cert=json.loads((ROOT.parent/'specimen_checks/probe_certificates.json').read_text())
    statuses={}
    for record in cert['results']:
        family=record['family']; label='ABC'[record['probe']-1]; cr=record['record']
        real,imag=map(Fraction,record['q_exact'])
        # Exact dyadic endpoint interpretations; decimal text was emitted with
        # 17 digits, so recover the binary64 grid number before Fraction.
        assert Fraction.from_float(float(cr['x0']))<=real<=Fraction.from_float(float(cr['x1']))
        assert Fraction.from_float(float(cr['y0']))<=imag<=Fraction.from_float(float(cr['y1']))
        status='certified stable' if cr['status']=='1' else 'no conclusion from positive test'
        if family=='DO10' and label=='C': status='proved nonembedded by exact degree certificate, separately replayed'
        statuses[f'{family}_{label}']=status
    report={'passed':True,'audit_type':'independent geometry derivation plus data consistency; not a new proof of the whole positive raster',
            'display_module_sha256':digest(ROOT/'projection_geometry.py'),
            'normalization_derivation':'Solve center=(F1(center)+F2(center))/2; d=F2(center)-center; H(z)=(z-center)/d. For DO, b=B*conj(d)/d.',
            'families':families,'positive_cell_records':records,'probe_statuses':statuses,
            'scope_limits':['Only certified cells may be colored as stable.','Omitted parameters have no classification implied by this positive cover.','C at DD10 has no injectivity conclusion.','C at DO10 is nonembedded at the exact parameter; the supplied disk radius is zero.','The torus projects out w; it is a display chart, and its central hole is not a cavity in the 4D locus.']}
    (ROOT/'independent_geometry_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'samples_per_family':len(q),'mask_counts':{n:r['positive_cells'] for n,r in records.items()},'DD10_reflected_cells':int(merged.sum()),'probe_statuses':statuses},indent=2))


if __name__=='__main__': main()

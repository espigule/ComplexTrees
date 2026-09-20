#!/usr/bin/env python3
"""Replay the exact historical-probe checks from this directory alone.

python run_specimen_checks.py          # interval membership and degree replay
python build_adaptive_specimens.py    # the exact-stopping display arrays
"""
from pathlib import Path
from fractions import Fraction as F
import subprocess, csv, json, hashlib, tempfile
from adaptive_degree_capture import verify_adaptive
from verify_adaptive_certificates import verify as independent_verify
from zipper_degree_capture import rational_upper_sqrt

ROOT=Path(__file__).resolve().parent
QS=[(F(1,2),F(1,5)),(F(1,2),F(2,5)),(F('0.34090875'),F('0.43484625'))]
FLAGS=['-O2','-std=c++17','-fno-fast-math','-ffp-contract=off']

def main():
    results=[]
    with tempfile.TemporaryDirectory(prefix='zipper-specimen-') as temporary:
        binary=Path(temporary)/'certify_zipper_cells'
        subprocess.run(['g++',*FLAGS,str(ROOT/'certify_zipper_cells.cpp'),'-o',str(binary)],check=True)
        for family in ['DD10','DD11','DO10','DO11']:
            for j,(x,y) in enumerate(QS,1):
                d=2**34;xx=int(x*d);yy=int(y*d)
                bounds=[xx/d,(xx+1)/d,yy/d,(yy+1)/d]
                assert F.from_float(bounds[0])<=x<=F.from_float(bounds[1])
                assert F.from_float(bounds[2])<=y<=F.from_float(bounds[3])
                out=ROOT/f'{family}_probe{j}.csv'
                subprocess.run([str(binary),family,str(out),'1','2000000','100',*[repr(z) for z in bounds]],check=True,capture_output=True,text=True)
                record=list(csv.DictReader(out.open()))[0]
                r1=rational_upper_sqrt(x*x+y*y);r2=rational_upper_sqrt((1-x)**2+y*y);r=max(r1,r2)
                delta=rational_upper_sqrt((x-F(1,2))**2+y*y)
                tails={str(n):str(r**n*delta/(1-r)) for n in [14,20,22,24,26]}
                results.append({'family':family,'probe':j,'q_exact':[str(x),str(y)],'box_dyadic_bits':34,
                    'certificate_file':out.name,'record':record,
                    'uniform_curve_tail_bounds_exact':tails,
                    'curve_tail_estimates':{n:float(F(v)) for n,v in tails.items()}})
                print(family,j,record['status'],flush=True)
    report={'source_sha256':hashlib.sha256((ROOT/'certify_zipper_cells.cpp').read_bytes()).hexdigest(),
        'compile_flags':' '.join(FLAGS),'budget':2000000,'maximum_depth':100,'results':results}
    (ROOT/'probe_certificates.json').write_text(json.dumps(report,indent=2)+'\n')
    x,y=QS[2]
    degree=verify_adaptive('DO10',x,y,'00011100','101101111111',F(0))
    assert degree['status']=='proved_nonembedded_parameter_disk' and degree['degree']==-1
    record_path=ROOT/'DO10_L_adaptive_exact.json'
    record_path.write_text(json.dumps(degree,indent=2)+'\n')
    # Separate implementation, exact Fractions, winding counted on another ray.
    independent=independent_verify(record_path)
    (ROOT/'DO10_L_independent_replay.json').write_text(json.dumps({'verified':True,'certificates':[independent]},indent=2)+'\n')
    print('DO10 C exact secondary contact independently verified.')

if __name__=='__main__':main()

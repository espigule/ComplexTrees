#!/usr/bin/env python3
"""Render whole-cell stable-slice certificates and canonical zipper curves.

The C++ interval verifier supplies five displayed inner-cover rasters.
The exact OO00, OO01 and OO10 disks and real interval are theorem-backed overlays. Gray is unresolved, never instability.
Curves are polygonal approximants only; their classifications come from the
separate interval records. All12 canonical curves use q*=3/8+i/4.
"""
from pathlib import Path
import json, hashlib, csv, io, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Polygon
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
STABLE='#16738E'; UNKNOWN='#E5EAED'; EXCLUDED='#F3DCD7'
INK='#193548'; GRAY='#627783'; RED='#B3262D'; BLUE='#164D8F'
N=512
TYPES=['DD','DO','OO']; SIGS=['00','01','10','11']
TESTED={'DD01','DD10','DD11','DO10','DO11','OO01','OO10'}
EXACT={'OO00','OO01','OO10'}
INNER_COVERS=TESTED-EXACT
STAR=.375+.25j
L=.34090875+.43484625j

def atomic_pdf(fig,filename,title):
    stream=io.BytesIO()
    fig.savefig(stream,format='pdf',metadata={'Title':title,'Author':'Bernat Espigulé','CreationDate':None,'ModDate':None})
    data=stream.getvalue()
    assert data.rstrip().endswith(b'%%EOF')
    assert len(PdfReader(io.BytesIO(data)).pages)==1
    final=ROOT/filename
    tmp=final.with_suffix('.pdf.tmp')
    tmp.write_bytes(data)
    assert tmp.stat().st_size==len(data)
    os.replace(tmp,final)

def read_mask(name):
    p=ROOT/(name+'_cells.csv')
    a=np.genfromtxt(p,delimiter=',',names=True)
    mask=np.zeros((N//2,N),dtype=bool)
    take=a['status']==1
    mask[a['iy'][take].astype(int),a['ix'][take].astype(int)]=True
    return mask, {'file':p.name,'certified_cells':int(take.sum()),
        'evaluated_cells':len(a),'maximum_word_depth':int(a['maxdepth'][take].max()),
        'minimum_squared_gap':float(a['minimum_squared_gap'][take].min()),
        'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def curve(kind,sig,q=STAR,depth=16):
    rev={'DD':(0,0),'DO':(0,1),'OO':(1,1)}[kind]
    e=[int(s) for s in sig]
    p=np.array([0j,1+0j])
    for _ in range(depth):
        a=p[::-1] if e[0] else p
        b=p[::-1] if e[1] else p
        if rev[0]:a=a.conjugate()
        if rev[1]:b=b.conjugate()
        a=e[0]*q+(-1)**e[0]*q*a
        b=q+e[1]*(1-q)+(-1)**e[1]*(1-q)*b
        assert abs(a[-1]-q)<1e-12 and abs(b[0]-q)<1e-12
        p=np.concatenate([a,b[1:]])
    assert abs(p[0])<1e-12 and abs(p[-1]-1)<1e-12
    return p

def draw_curve(ax,kind,sig,fontsize=10,q=STAR,depth=16):
    p=curve(kind,sig,q,depth)
    mid=(len(p)-1)//2
    ax.plot(p[:mid+1].real,p[:mid+1].imag,color=RED,lw=.55)
    ax.plot(p[mid:].real,p[mid:].imag,color=BLUE,lw=.55)
    ax.scatter([0,1],[0,0],s=6,color=INK,zorder=5)
    ax.scatter([q.real],[q.imag],s=8,color=INK,zorder=5)
    ax.set_aspect('equal',adjustable='datalim');ax.margins(.07);ax.axis('off')
    return p

def draw_slice(ax,kind,sig,masks,fontsize=8,show_y=False,mark=True):
    ax.set_facecolor('white')
    # Only the strict contraction lens consists of these IFS parameters.
    # Its exterior must not be labelled a nonembedded attractor family.
    yy=np.linspace(0,.52,600)
    xl=1-np.sqrt(1-yy*yy); xr=np.sqrt(1-yy*yy)
    poly=np.column_stack([np.r_[xl,xr[::-1]],np.r_[yy,yy[::-1]]])
    ax.add_patch(Polygon(poly,facecolor=EXCLUDED,edgecolor='none',zorder=-1))
    ax.plot(xl,yy,color='#C6CFD4',lw=.55,zorder=0)
    ax.plot(xr,yy,color='#C6CFD4',lw=.55,zorder=0)
    disk=Circle((.5,0),.5,facecolor=UNKNOWN,edgecolor='none',zorder=0)
    ax.add_patch(disk)
    name=kind+sig
    if name in EXACT:
        ax.add_patch(Circle((.5,0),.5,facecolor=STABLE,edgecolor='none',zorder=1))
    elif name in masks:
        image=np.zeros((N//2,N,4))
        from matplotlib.colors import to_rgba
        image[masks[name]]=to_rgba(STABLE)
        ax.imshow(image,extent=(0,1,0,.5),origin='lower',interpolation='none',zorder=2)
    theta=np.linspace(0,np.pi,501)
    ax.plot(.5+.5*np.cos(theta),.5*np.sin(theta),color='#B27669',lw=.7,zorder=3)
    ax.plot([0,1],[0,0],color=STABLE,lw=2.2,zorder=4,solid_capstyle='butt')
    ax.scatter([0,1],[0,0],s=fontsize*.75,facecolor='white',edgecolor=STABLE,lw=.55,zorder=8,clip_on=False)
    if mark:
        ax.scatter([STAR.real],[STAR.imag],s=fontsize*.9,color=INK,edgecolor='white',lw=.45,zorder=6)
    if name=='DD11':
        ax.scatter([L.real],[L.imag],marker='*',s=fontsize*3.7,color='#5B3572',edgecolor='white',lw=.4,zorder=7)
        ax.annotate('$L$',(L.real,L.imag),xytext=(-9,4),textcoords='offset points',fontsize=fontsize,color='#5B3572',ha='right')
    if name in EXACT:
        ax.text(.53,.135,'exact disk',ha='center',fontsize=fontsize*.86,color='white')
    elif name not in TESTED:
        ax.text(.52,.115,'unresolved',ha='center',fontsize=fontsize*.81,color=GRAY)
    ax.set_xlim(0,1);ax.set_ylim(-.004,.52);ax.set_aspect('equal',adjustable='box')
    ax.set_xticks([0,.5,1],['0','½','1'],fontsize=fontsize*.83)
    ax.set_yticks([0,.25,.5],['0','¼','½'] if show_y else [],fontsize=fontsize*.83)
    ax.tick_params(length=2.5,pad=2,colors=GRAY)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color('#A4B2BB')

def legend(fig,y,fs):
    handles=[Rectangle((0,0),1,1,fc=STABLE),Rectangle((0,0),1,1,fc=UNKNOWN),Rectangle((0,0),1,1,fc=EXCLUDED),Rectangle((0,0),1,1,fc='white',ec='#C6CFD4')]
    labels=['Proved stable','Unresolved','Proved nonembedded','Outside contraction lens']
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.53,y),ncol=4,frameon=False,fontsize=fs,
               handlelength=1.05,columnspacing=1.1)

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans','pdf.fonttype':42,
       'ps.fonttype':42,'figure.facecolor':'white','path.simplify':False})
    masks={};info={}
    for name in sorted(TESTED):masks[name],info[name]=read_mask(name)
    # The endpoint-reflection theorem, followed by conjugation, maps q to
    # 1-conjugate(q). Unite independently certified reflected companions.
    for a,b in [('DD01','DD10'),('OO01','OO10')]:
        merged=masks[a]|masks[b][:,::-1]
        masks[a]=merged;masks[b]=merged[:,::-1]
    for name in masks:
        np.save(ROOT/(name+'_mask.npy'),masks[name])
        info[name]['certificate_cells_after_reflection']=int(masks[name].sum())
        info[name]['display_mode']='exact_disk' if name in EXACT else 'certified_cell_cover'
    for name in sorted(EXACT):
        info.setdefault(name,{})
        info[name].update({'display_mode':'exact_disk','full_stable_set':'|q-1/2|<1/2',
          'proof':'Invariant triangle with first-piece intersection exactly q; common necessary-disk theorem.',
          'proof_label':'prop:oo00-stable-disk' if name=='OO00' else 'thm:oo01-oo10-stable-disks'})
    # Poster: parameter panel plus one true family-specific curve in each cell.
    fig=plt.figure(figsize=(60/2.54,32/2.54))
    left=.066;cw=.229;top=.88;rh=.248
    for c,sig in enumerate(SIGS):
        fig.text(left+c*cw+.100,.957,r'$e=('+sig[0]+','+sig[1]+')$',ha='center',fontsize=30,color=INK)
    for r,kind in enumerate(TYPES):
        cy=top-r*rh
        fig.text(.021,cy-.06,kind,fontsize=29,color=INK,fontweight='bold',va='center',rotation=90)
        for c,sig in enumerate(SIGS):
            x=left+c*cw
            draw_slice(fig.add_axes((x,cy-.159,.162,.205)),kind,sig,masks,fontsize=26,show_y=c==0)
            draw_curve(fig.add_axes((x+.167,cy-.14,.057,.154)),kind,sig,fontsize=20)
            if r==0:fig.text(x+.194,cy+.036,'$K_{q_*}$',ha='center',fontsize=22,color=GRAY)
    legend(fig,.096,25)
    fig.text(.5,.073,r'Upper half-plane; reflect across the real axis.  $q_*=3/8+i/4$ links each dot to its curve.',ha='center',fontsize=22,color=GRAY)
    fig.text(.5,.028,'Three exact stable disks; five further slices have certified blue cells. Grey remains unresolved.  L: Koch–wiggle island.',ha='center',fontsize=21,color=GRAY)
    atomic_pdf(fig,'zipper_stability_atlas.pdf','Twelve marked zipper families and certified stable slices')
    fig.savefig(ROOT/'zipper_stability_atlas.png',dpi=130)
    plt.close(fig)
    # Paper: same data, larger relative axes typography, no tiny curve insets.
    fig=plt.figure(figsize=(15.6/2.54,12.8/2.54))
    left=.095;cw=.221;top=.83;rh=.247
    for c,sig in enumerate(SIGS):fig.text(left+c*cw+.098,.947,r'$e=('+sig[0]+','+sig[1]+')$',ha='center',fontsize=9,color=INK)
    for r,kind in enumerate(TYPES):
        cy=top-r*rh
        fig.text(.024,cy-.017,kind,fontsize=9,fontweight='bold',rotation=90,color=INK,va='center')
        for c,sig in enumerate(SIGS):draw_slice(fig.add_axes((left+c*cw,cy-.115,.197,.223)),kind,sig,masks,fontsize=8.3,show_y=c==0)
    legend(fig,.016,7.2)
    atomic_pdf(fig,'zipper_stability_atlas_paper.pdf','Stable slices of twelve marked binary zippers')
    fig.savefig(ROOT/'zipper_stability_atlas_paper.png',dpi=230);plt.close(fig)
    # Curve gallery, same exact dyadic q across all12 presentations.
    fig=plt.figure(figsize=(15.6/2.54,13.3/2.54))
    for c,sig in enumerate(SIGS):fig.text(.17+c*.227,.95,r'$e=('+sig[0]+','+sig[1]+')$',ha='center',fontsize=9,color=INK)
    for r,kind in enumerate(TYPES):
        fig.text(.026,.805-r*.266,kind,rotation=90,fontsize=10,fontweight='bold',color=INK,va='center')
        for c,sig in enumerate(SIGS):
            draw_curve(fig.add_axes((.069+c*.227,.69-r*.266,.211,.218)),kind,sig)
            status='proved stable' if kind+sig in TESTED or kind+sig=='OO00' else 'unresolved'
            fig.text(.173+c*.227,.669-r*.266,status,fontsize=6.8,color=STABLE if status=='proved stable' else GRAY,ha='center')
    fig.text(.5,.025,r'Common parameter $q_*=3/8+i/4$; red and blue show the first two pieces.',ha='center',fontsize=7.6,color=GRAY)
    atomic_pdf(fig,'zipper_curve_gallery.pdf','The twelve canonical zipper curves at one common parameter')
    fig.savefig(ROOT/'zipper_curve_gallery.png',dpi=230);plt.close(fig)
    r=max(abs(STAR),abs(1-STAR));tail=r**16*abs(STAR-.5)/(1-r)
    record={'grid':{'Nx':N,'Ny':N//2,'q_rectangle':[0,1,0,.5],'cell_width':1/N,
       'status_semantics':{'stable':'OO00, OO01 and OO10 are exact open stable disks. Five further slices display closed cells proved embedded by outward interval tests; all real intervals are exact.',
          'unresolved':'No conclusion about injectivity. Four presentations were not covered by the seven-return table.',
          'nonembedded':'Inside the strict contraction lens and outside or on the common necessary diameter disk; these systems still have connected attractors.',
          'outside_contraction_lens':'White. Parameters not belonging to the strict contraction domain are not classified as attractors.'}},
       'families':info,'exact_complete_slices':sorted(EXACT),'displayed_computed_inner_covers':sorted(INNER_COVERS),
       'symmetry_completion':['DD01 <-> DD10 under q->1-conjugate(q)'],
       'retained_unused_display_records':['OO01/OO10 cell masks and reflection remain valid, but full disks are displayed by theorem.'],
       'specimens':{'q_exact':['3/8','1/4'],'depth':16,'vertices_per_curve':65537,'uniform_tail_bound_numerical':tail,
          'uniform_tail_bound_exact':'(sqrt(29)/8)^16*(sqrt(5)/8)/(1-sqrt(29)/8)',
          'interpretation':'Canonical polygonal approximants. Three labels follow from exact OO disk theorems and five from independent point certificates, not visual simplicity.'},
       'island':{'q_exact_decimal':['0.34090875','0.43484625'],'point_record':'DD11_island_point.csv',
          'source':'Danny Calegari, Wiggle Island, arXiv:2205.11442v2, §2.5; disconnected-component assertion is cited, not inferred from this grid.'},
       'implementation':{'interval_verifier':'certify_zipper_cells.cpp','compiler_flags':'-O2 -std=c++17 -fno-fast-math -ffp-contract=off',
          'arithmetic':'Outward binary64 intervals using nextafter after every elementary operation; all main-grid q endpoints are exact dyadics.'}}
    (ROOT/'zipper_stability_metadata.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Created all three PDF/PNG pairs; uniform curve-tail bound',tail)

if __name__=='__main__':main()

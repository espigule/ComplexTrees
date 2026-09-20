#!/usr/bin/env python3
"""Certified stable sheets in the four-dimensional exterior atlas.

Only DD10, DD11, DO10 and DO11 are displayed. Empty backgrounds carry no
classification. Stable pixels come from the supplied outward interval records;
floating-point projection is a display calculation, not a new certificate.
"""
from __future__ import annotations
from pathlib import Path
import argparse, hashlib, io, json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
from pypdf import PdfReader, PdfWriter
from PIL import Image

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'zippers'
sys.path.insert(0,str(ROOT))
from projection_geometry import project_q, coefficient_torus, exterior_coefficients

FAMILIES=('DD10','DD11','DO10','DO11')
PARAMS={'A':.5+.2j,'B':.5+.4j,'C':.34090875+.43484625j}
INK='#193548'; MUTED='#637681'; GRID='#D9E1E4'
RED='#B12B36'; BLUE='#174D8B'; VIOLET='#683A81'
COLORS=['#253A78','#255BA1','#338CA7','#77B7AE','#E4C55E']
CMAP=LinearSegmentedColormap.from_list('stable_weight',COLORS)
NORM=Normalize(0,1)
plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans',
    'pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white',
    'path.simplify':False,'axes.unicode_minus':True})

def atomic_bytes(path,data):
    path=Path(path); tmp=path.with_name(path.name+'.tmp')
    tmp.write_bytes(data)
    if tmp.stat().st_size!=len(data): raise IOError('Incomplete output')
    os.replace(tmp,path)

def save_figure(fig,stem,title,dpi=220):
    stream=io.BytesIO()
    fig.savefig(stream,format='pdf',dpi=300 if 'poster' in stem else 600,metadata={'Title':title,'Author':'Bernat Espigulé','CreationDate':None,'ModDate':None})
    data=stream.getvalue(); assert data.rstrip().endswith(b'%%EOF')
    assert len(PdfReader(io.BytesIO(data)).pages)==1
    atomic_bytes(ROOT/(stem+'.pdf'),data)
    stream=io.BytesIO(); fig.savefig(stream,format='png',dpi=dpi)
    data=stream.getvalue(); Image.open(io.BytesIO(data)).load()
    atomic_bytes(ROOT/(stem+'.png'),data)
    print(stem, len(data),flush=True)

def read_masks():
    masks={}; metadata={}
    for f in FAMILIES:
        a=np.genfromtxt(SOURCE/(f+'_cells.csv'),delimiter=',',names=True)
        raw=np.zeros((256,512),dtype=bool)
        take=a['status']==1
        raw[a['iy'][take].astype(int),a['ix'][take].astype(int)]=True
        sources=[f+'_cells.csv']
        if f=='DD10':
            b=np.genfromtxt(SOURCE/'DD01_cells.csv',delimiter=',',names=True)
            reflected=np.zeros_like(raw); yes=b['status']==1
            reflected[b['iy'][yes].astype(int),511-b['ix'][yes].astype(int)]=True
            raw|=reflected; sources.append('DD01_cells.csv')
        stored=np.load(SOURCE/(f+'_mask.npy'))
        assert np.array_equal(raw,stored),f
        masks[f]=raw
        metadata[f]={'certified_closed_q_cells':int(raw.sum()),
          'displayed_original_files':sources,
          'source_sha256':{s:hashlib.sha256((SOURCE/s).read_bytes()).hexdigest() for s in sources},
          'reflection_completion':'q -> 1-conjugate(q)' if f=='DD10' else None}
    return masks,metadata

def block_mesh(mask,family,block=2):
    """Use certified 2x2 blocks; planar display faces approximate curved images."""
    ny,nx=mask.shape
    valid=mask.reshape(ny//block,block,nx//block,block).all(axis=(1,3))
    iy,ix=np.nonzero(valid)
    x0=ix*block/512; y0=iy*block/512; step=block/512
    q=np.stack((x0+1j*y0,x0+step+1j*y0,x0+step+1j*(y0+step),x0+1j*(y0+step)),axis=-1)
    q=np.concatenate((q,q.conjugate()),axis=0)
    if family=='DD11':
        d=json.loads((SOURCE/'DD11_island_zoom_grid.json').read_text())
        u,v,s,t=d['bounds_binary64_exact']
        extra=np.array([[u+1j*s,v+1j*s,v+1j*t,u+1j*t]])
        q=np.concatenate((q,extra,extra.conjugate()),axis=0)
    xyz,atlas=project_q(q,family)
    faces=xyz; w=atlas['w'].mean(axis=1)
    return faces,w,{'positive_blocks':int(len(ix)),'block_side_q':step,
        'q_conjugation_added':True,'displayed_q_cells':int(len(ix)*block*block),'DD11_extra_island_square_added':family=='DD11'}

def torus_context(ax,fontsize=8):
    """The left-coefficient half of the universal rho=2 reference surface."""
    t1=np.linspace(np.pi/2,3*np.pi/2,85)
    t2=np.linspace(-np.pi,np.pi,140)
    for a in np.linspace(np.pi/2,3*np.pi/2,7):
        xyz=coefficient_torus(a+np.zeros_like(t2),t2,2+np.zeros_like(t2))
        ax.plot(*xyz.T,color=GRID,lw=.42,zorder=0,alpha=.65)
    for b in np.linspace(-np.pi,np.pi,13)[:-1]:
        xyz=coefficient_torus(t1,b+np.zeros_like(t1),2+np.zeros_like(t1))
        ax.plot(*xyz.T,color=GRID,lw=.42,zorder=0,alpha=.65)
    ax.set(xlim=(-4.3,4.3),ylim=(-4.3,4.3),zlim=(-2.12,2.12))
    ax.set_box_aspect((1,1,.53),zoom=1.52 if fontsize>15 else 1.03); ax.view_init(elev=27,azim=-57)
    ax.set_proj_type('ortho'); ax.set_axis_off()

def draw_projection(ax,families,meshes,fontsize=8,markers=True,labels=True):
    torus_context(ax,fontsize)
    for f in families:
        faces,w,_=meshes[f]
        coll=Poly3DCollection(faces,facecolors=CMAP(w),edgecolors='none',linewidths=0,
                              antialiased=False,rasterized=True,zsort='average')
        ax.add_collection3d(coll)
        if markers:
            for k,q in PARAMS.items():
                xyz,_=project_q(q,f)
                ax.scatter(*xyz,s=fontsize*2.2,facecolors='white',edgecolors=INK,linewidths=.7,
                           depthshade=False,zorder=100)
                if labels:
                    px,py,_=proj3d.proj_transform(*xyz,ax.get_proj())
                    if f.endswith('10'):
                        loc={'A':(.93,.80),'B':(.93,.58),'C':(.93,.36)}[k]
                    elif f=='DD11':
                        loc={'A':(.045,.82),'B':(.045,.60),'C':(.045,.38)}[k]
                    else:
                        loc={'A':(.045,.80),'B':(.94,.62),'C':(.94,.40)}[k]
                    ax.annotate(k,xy=(px,py),xycoords='data',xytext=loc,textcoords='axes fraction',
                        fontsize=fontsize*.88,ha='center',va='center',color=INK,
                        arrowprops={'arrowstyle':'-','lw':.60,'color':MUTED},zorder=210)
        if len(families)>1:
            xyz,_=project_q(.5+.14j,f)
            px,py,_=proj3d.proj_transform(*xyz,ax.get_proj())
            loc=(.15,.80) if f.endswith('11') else (.86,.37)
            ax.annotate(f[:2]+r' $('+f[-2]+','+f[-1]+')$',xy=(px,py),xycoords='data',
                xytext=loc,textcoords='axes fraction',ha='center',fontsize=fontsize,
                color=INK,arrowprops={'arrowstyle':'-','lw':.55,'color':MUTED},zorder=200)
    ax.text2D(.5,.008,r'$(\theta_1,\theta_2,\rho)$; colour records $w$',transform=ax.transAxes,
              ha='center',fontsize=fontsize*.91,color=MUTED)

def draw_exterior(ax,family,mask,fontsize=8,show_island=False):
    qx=np.arange(513)/512; qy=np.arange(257)/512
    qq=qx[None,:]+1j*qy[:,None]
    # q=0 touches no certified cell. A finite harmless coordinate avoids a
    # nonfinite pcolormesh corner; that corner remains masked throughout.
    qq[0,0]=1e-12
    cc=1/qq
    center=(qx[:-1][None,:]+1/1024)+1j*(qy[:-1,None]+1/1024)
    w=abs(center)/(abs(center)+abs(1-center))
    value=np.ma.array(w,mask=~mask)
    ax.pcolormesh(cc.real,cc.imag,value,cmap=CMAP,norm=NORM,shading='flat',
                  edgecolors='none',antialiased=False,rasterized=True,zorder=2)
    for k,q in PARAMS.items():
        c=1/q
        ax.scatter([c.real],[c.imag],s=fontsize*2.5,facecolor='white',edgecolor=INK,lw=.65,zorder=6)
        offset={'A':(6,4),'B':(6,3),'C':(6,-5)}[k]
        ax.annotate(k,(c.real,c.imag),xytext=offset,textcoords='offset points',fontsize=fontsize,
                    color=INK,zorder=7,ha='left',va='center')
    ax.set_xlim(.96,4.65);ax.set_ylim(-2.70,.13)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xticks([1,2,3,4]);ax.set_yticks([-2,-1,0])
    ax.tick_params(labelsize=fontsize*.9,length=2,pad=2,colors=MUTED)
    ax.set_xlabel(r'$\operatorname{Re}c$',fontsize=fontsize,labelpad=1,color=INK)
    ax.set_ylabel(r'$\operatorname{Im}c$',fontsize=fontsize,labelpad=1,color=INK)
    ax.spines[['top','right']].set_visible(False)
    for side in ('left','bottom'):ax.spines[side].set_color('#9BAAB2');ax.spines[side].set_linewidth(.55)
    ax.grid(color='#ECF0F2',linewidth=.35,zorder=-5)
    ax.text(.96,.96,r'$c=1/q$',transform=ax.transAxes,ha='right',va='top',fontsize=fontsize,color=INK)
    if show_island:
        ax.add_patch(Rectangle((.54,.025),.445,.775,transform=ax.transAxes,facecolor='white',edgecolor='#C0CBD1',lw=.5,zorder=20))
        ins=ax.inset_axes([.63,.20,.30,.41],zorder=21)
        draw_island(ins,max(8,fontsize*.88),compact=True)
        ax.text(.765,.70,'C: island',transform=ax.transAxes,ha='center',fontsize=fontsize*.9,color=VIOLET,zorder=25)

def island_data():
    d=json.loads((SOURCE/'DD11_island_zoom_grid.json').read_text());x0,x1,y0,y1=d['bounds_binary64_exact']
    x=np.linspace(x0,x1,d['N']+1);y=np.linspace(y0,y1,d['N']+1)
    q=x[None,:]+1j*y[:,None];c=1/q;c0=1/PARAMS['C']
    p=np.genfromtxt(SOURCE/'DD11_island_zoom.csv',delimiter=',',names=True)
    assert np.all(p['status']==1) and len(p)==64*64
    return q,c,c0

def draw_island(ax,fontsize=7,compact=False):
    q,c,c0=island_data()
    center=(q[:-1,:-1]+q[1:,1:])/2
    w=abs(center)/(abs(center)+abs(1-center))
    ax.pcolormesh((c.real-c0.real)*1e6,(c.imag-c0.imag)*1e6,w,cmap=CMAP,norm=NORM,
                  shading='flat',edgecolors='none',rasterized=True)
    ax.scatter([0],[0],marker='*',s=fontsize*4.5,color=VIOLET,edgecolor='white',lw=.5,zorder=10)
    ax.set_xlim(-2,2);ax.set_ylim(-2,2);ax.set_aspect('equal')
    ax.set_xticks([] if compact else [-1,1]);ax.set_yticks([-1,1]);ax.tick_params(labelsize=fontsize,length=1.5,pad=1,colors=MUTED)
    for s in ax.spines.values():s.set_linewidth(.5);s.set_color('#C0CBD1')
    ax.set_xlabel(r'$10^6(c-c_C)$',fontsize=fontsize,color=MUTED,labelpad=2)

def make_curve(family,q,depth):
    e=(1,int(family[-1]));conj=(False,family.startswith('DO'))
    p=np.array([0j,1+0j])
    for _ in range(depth):
        branches=[]
        for j in (0,1):
            x=p[::-1] if e[j] else p
            if conj[j]:x=x.conjugate()
            if j==0:x=e[0]*q+(-1)**e[0]*q*x
            else:x=q+e[1]*(1-q)+(-1)**e[1]*(1-q)*x
            branches.append(x)
        p=np.r_[branches[0],branches[1][1:]]
    return p

def load_curve(family,label,preview=False):
    file=(ROOT/'specimen_checks'/(family+'_C_adaptive_curve.npz')) if label=='C' else ROOT/'curve_data'/(family+'_'+label+'.npz')
    if file.exists():
        d=np.load(file)
        for key in ['points','z','curve','p']:
            if key in d:
                p=d[key];join=int(d['first_level_join_index']) if 'first_level_join_index' in d else (len(p)-1)//2
                return p,join
        raise KeyError(str(list(d.keys())))
    p=make_curve(family,PARAMS[label],17 if preview else 20)
    return p,(len(p)-1)//2

def draw_curve(ax,family,label,fontsize=8,preview=False):
    p,m=load_curve(family,label,preview)
    q=PARAMS[label]
    # Adaptive source traversal keeps the mandatory first join as a vertex.
    assert abs(p[m]-q)<1e-10
    ax.plot(p[:m+1].real,p[:m+1].imag,color=RED,lw=.47,rasterized=True)
    ax.plot(p[m:].real,p[m:].imag,color=BLUE,lw=.47,rasterized=True)
    ax.scatter([0,1],[0,0],s=fontsize*.75,color=INK,zorder=10)
    ax.scatter([q.real],[q.imag],s=fontsize*.9,facecolors='white',edgecolors=INK,lw=.5,zorder=10)
    ax.set_aspect('equal',adjustable='datalim');ax.margins(.065);ax.axis('off')
    if label=='C' and family=='DD11':
        ax.text(.5,-.055,'Wiggle Island',transform=ax.transAxes,ha='center',color=VIOLET,fontsize=fontsize*.91)
    elif label=='C' and family=='DO10':
        ax.text(.5,-.055,'Extra contact',transform=ax.transAxes,ha='center',color=MUTED,fontsize=fontsize*.85)
    elif label=='C' and family=='DO11':
        ax.text(.5,-.055,'certified stable',transform=ax.transAxes,ha='center',color=MUTED,fontsize=fontsize*.85)

def weight_bar(fig,rect,fontsize=8,label='Fourth coordinate: '+r'$w=|q|/(|q|+|1-q|)$'):
    bar=fig.colorbar(ScalarMappable(norm=NORM,cmap=CMAP),cax=fig.add_axes(rect),orientation='horizontal')
    bar.set_ticks([0,.25,.5,.75,1]);bar.ax.tick_params(labelsize=fontsize,length=2,pad=2,colors=MUTED)
    bar.outline.set_linewidth(.45);bar.outline.set_edgecolor('#B6C3CB')
    bar.set_label(label,fontsize=fontsize,color=INK,labelpad=2)

def build_projection(meshes):
    fig=plt.figure(figsize=(15.6/2.54,10.2/2.54))
    for i,kind in enumerate(('DD','DO')):
        ax=fig.add_axes([.005+.50*i,.18,.49,.77],projection='3d')
        draw_projection(ax,[kind+'10',kind+'11'],meshes,8.8,markers=False)
        fig.text(.25+.5*i,.958,'Direct / direct' if kind=='DD' else 'Direct / opposite',ha='center',fontsize=9.5,fontweight='bold',color=INK)
    fig.text(.5,.155,r'Canonical exterior coefficients: $|c_1|,|c_2|>1$',ha='center',fontsize=9,color=INK)
    weight_bar(fig,[.24,.098,.52,.022],8)
    save_figure(fig,'four_family_projection','Four zipper sheets in a 3D projection of the 4D atlas',260);plt.close(fig)

def build_family_paper(kind,masks,preview=False):
    fig=plt.figure(figsize=(15.6/2.54,14.5/2.54))
    columns=[.030,.355,.565,.775];widths=[.292,.192,.192,.192]
    fig.text(.162,.968,r'Exterior stable slice $c=1/q$',ha='center',fontsize=8.8,color=INK)
    for j,k in enumerate(PARAMS):
        fig.text(columns[j+1]+widths[j+1]/2,.968,k,ha='center',fontsize=11,color=INK,fontweight='bold')
    for r,f in enumerate((kind+'10',kind+'11')):
        top=.922-r*.405;bottom=top-.305
        fig.text(.03,top,r'$\mathrm{'+kind+r'}\ ('+f[-2]+','+f[-1]+')$',ha='left',fontsize=11,color=INK,fontweight='bold')
        draw_exterior(fig.add_axes([columns[0]+.035,bottom,.262,.26]),f,masks[f],9.0,show_island=f=='DD11')
        for j,k in enumerate(PARAMS):
            draw_curve(fig.add_axes([columns[j+1],bottom-.004,widths[j+1],.305]),f,k,8.8,preview)
    fig.text(.5,.140,r'$q_A=0.5+0.2i\qquad q_B=0.5+0.4i$',ha='center',fontsize=8.6,color=INK)
    fig.text(.5,.115,r'$q_C=L=0.34090875+0.43484625i$',ha='center',fontsize=8.6,color=INK)
    weight_bar(fig,[.245,.074,.51,.013],8.0,label=r'$w$ (the fourth atlas coordinate)')
    save_figure(fig,'four_family_'+kind,kind+' exterior stable slices and three shared zipper specimens',260);plt.close(fig)

def build_poster(masks,meshes,preview=False):
    fig=plt.figure(figsize=(60/2.54,42/2.54))
    fig.text(.135,.972,'IN THE 4D ATLAS',ha='center',fontsize=27,fontweight='bold',color=INK)
    fig.text(.38,.972,'EXTERIOR SLICE',ha='center',fontsize=27,fontweight='bold',color=INK)
    for x,k in zip((.598,.757,.913),PARAMS):
        fig.text(x,.974,k,ha='center',fontsize=31,fontweight='bold',color=INK)
    fig.text(.38,.945,r'$c=1/q$',ha='center',fontsize=24,color=MUTED)
    for x,s in zip((.598,.757,.913),(r'$q=0.5+0.2i$',r'$q=0.5+0.4i$',r'$q=L$')):
        fig.text(x,.945,s,ha='center',fontsize=22,color=MUTED)
    for r,f in enumerate(FAMILIES):
        top=.910-r*.200;bottom=top-.171
        fig.text(.02,top,f[:2]+' ('+f[-2]+','+f[-1]+')',ha='left',va='top',fontsize=29,fontweight='bold',color=INK)
        ax=fig.add_axes([.012,bottom-.010,.252,.18],projection='3d')
        draw_projection(ax,[f],meshes,23,markers=True,labels=True)
        draw_exterior(fig.add_axes([.288,bottom+.005,.205,.155]),f,masks[f],23,show_island=f=='DD11')
        for j,k in enumerate(PARAMS):
            draw_curve(fig.add_axes([.526+j*.157,bottom+.007,.145,.152]),f,k,23,preview)
    weight_bar(fig,[.075,.068,.365,.014],20,label=r'Fourth coordinate $w=|q|/(|q|+|1-q|)$')
    fig.text(.728,.076,r'$|c_1|,|c_2|>1$'+'  |  '+r'$L=0.34090875+0.43484625i$',ha='center',fontsize=23,color=INK)
    fig.text(.728,.044,'Curves: red and blue are the first two pieces.',ha='center',fontsize=21,color=MUTED)
    fig.text(.5,.008,'Coloured parameter regions: certified stable cover. A, B, C link to curve specimens.',ha='center',fontsize=21,color=MUTED)
    save_figure(fig,'four_family_poster','Four exterior stable slices, their 4D placement and twelve zipper specimens',145);plt.close(fig)

def build_island():
    fig=plt.figure(figsize=(7.6/2.54,7.5/2.54))
    fig.text(.5,.94,'The DD (1,1) island in exterior coordinates',ha='center',fontsize=8.6,color=INK)
    ax=fig.add_axes([.19,.30,.66,.56]);draw_island(ax,8)
    fig.text(.5,.07,r'$c_C=1/(0.34090875+0.43484625i)$',ha='center',fontsize=8.2,color=INK)
    fig.text(.5,.025,'4,096 certified cells; the star is C.',ha='center',fontsize=8,color=MUTED)
    save_figure(fig,'DD11_island_exterior','An exterior view of the certified DD11 island neighborhood',300);plt.close(fig)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');parser.add_argument('--only',choices=['all','projection','paper','poster','island'],default='all');args=parser.parse_args()
    masks,meta=read_masks();meshes={}
    for f in FAMILIES:meshes[f]=block_mesh(masks[f],f);meta[f]['projection_mesh']=meshes[f][2]
    if args.only in ('all','projection'):build_projection(meshes)
    if args.only in ('all','paper'):
        for k in ('DD','DO'):build_family_paper(k,masks,args.preview)
    if args.only in ('all','poster'):build_poster(masks,meshes,args.preview)
    if args.only in ('all','island'):build_island()
    writer=PdfWriter()
    for k in ('DD','DO'):writer.append(str(ROOT/('four_family_'+k+'.pdf')))
    stream=io.BytesIO();writer.write(stream);atomic_bytes(ROOT/'four_family_gallery_paper.pdf',stream.getvalue())
    record={'displayed_families':list(FAMILIES),'source_grid':[512,256],
      'positive_certificate_sources':meta,'mask_semantics':'Closed positive q-cell cover; empty background has no membership meaning. Planar faces and pcolormesh corner chords are piecewise-linear display approximations of the nonlinear transported cells.',
      'projection':{'coordinates':'X=(4+rho*cos(theta1))*cos(theta2),Y=(4+rho*cos(theta1))*sin(theta2),Z=rho*sin(theta1)',
         'theta1':'arg canonical a=-arg canonical c1','theta2':'arg canonical b=-arg canonical c2','fourth_coordinate':'w=abs(q)/(abs(q)+abs(1-q))',
         'ghost_reference':'rho=2; only theta1 in [pi/2,3pi/2] is drawn','not_a_geometric_hole':'The toroidal hole belongs to the coordinate display.',
         'quantities':'rho=2/(abs(q)+abs(1-q)); exterior abs(c1),abs(c2)>1',
         'DO_normalization':'projection_geometry.py contains the phase correction from endpoint to canonical charts'},
      'exterior_slice':{'coordinate':'unsigned endpoint c=1/q','window':[.96,4.65,-2.7,.13],
         'q_half':'Im q>=0, hence Im c<=0; other half follows by conjugation','full_slice_is_unbounded':True},
      'shared_specimens':{k:{'q_real':v.real,'q_imag':v.imag} for k,v in PARAMS.items()},
      'C_classification':{'DD10':'No membership assertion shown','DD11':'Certified stable; published Wiggle Island result','DO10':'Certified secondary contact by exact adaptive degree -1','DO11':'Certified stable'},
      'curve_approximation':'C: adaptive chord approximation with uniform error <5e-4, exact source traversal and first_level_join_index, no simplification; see specimen_checks/adaptive_curve_metadata.json. A/B: fixed depth '+str(17 if args.preview else 20)+'; preview flag is '+str(args.preview),
      'DD11_island':'The extra closed 2^-20 square consists of 64x64 certified closed cells, transported by c=1/q.',
      'figure_outputs':['four_family_projection.pdf','four_family_DD.pdf','four_family_DO.pdf','four_family_gallery_paper.pdf','four_family_poster.pdf','DD11_island_exterior.pdf']}
    atomic_bytes(ROOT/'four_family_figure_metadata.json',(json.dumps(record,indent=2)+'\n').encode())

if __name__=='__main__':main()

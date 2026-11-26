"""Render the current beta figure from fixed coordinates, distances and current test families.
No distance model or permutation test is fitted by this renderer.
"""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator,StrMethodFormatter
from matplotlib.transforms import Bbox
from paths import ROOT,INP,OUT,QA
ORDER=('B1','B2','B3','M1','M2','M3','A1','A2','A3')
CATEGORIES=('Low','Moderate','High')
COLOR=dict(zip(CATEGORIES,('#007C83','#D28E00','#9B4263')))
CATEGORY={s:CATEGORIES[i//3] for i,s in enumerate(ORDER)}
LOCALITY={s:('SI' if s in ('M2','M3','A3') else 'J') for s in ORDER}
SHAPE={'J':'o','SI':'^'}
PATHS={'s16_coords': 's16_coords.tsv', 's16_axes': 's16_axes.tsv', 's16_dispersion': 's16_dispersion.tsv', 's16_statistics': 's16_statistics.json', 'its_coords': 'its_coords.tsv', 'its_axes': 'its_axes.tsv', 'its_dispersion': 'its_dispersion.tsv', 'its_dispersion_statistics': 'its_dispersion_statistics.tsv'}

def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':8.5,
        'axes.titlesize':10,'xtick.labelsize':8.5,'ytick.labelsize':8.5,
        'legend.fontsize':8.5,'axes.spines.top':False,'axes.spines.right':False,
        'axes.linewidth':.65,'xtick.major.width':.6,'ytick.major.width':.6,
        'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white',
        'axes.facecolor':'white','savefig.facecolor':'white','axes.axisbelow':True})

def axes_grid(height, *, top=.85,bottom=.16,hspace=.78):
    fig, axes=plt.subplots(2,2,figsize=(180/25.4,height/25.4))
    fig.subplots_adjust(left=.13,right=.97,top=top,bottom=bottom,hspace=hspace,wspace=.48)
    for ax in axes.flat:
        ax.grid(axis='y',color='#E8E8E8',lw=.5)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    return fig,axes

def heading(fig,title,subtitle):
    fig.text(.13,.975,title,fontsize=12,fontweight='bold',va='top')
    fig.text(.13,.946,subtitle,fontsize=8.5,va='top')

def panel(ax,letter,title,annotation=''):
    ax.set_title(letter+'  '+title,loc='left',fontweight='bold',fontsize=10,pad=35 if annotation else 12)
    if annotation: ax.text(0,1.025,annotation,transform=ax.transAxes,va='bottom',fontsize=8.5)

def common_legend(fig):
    handles=[Line2D([],[],marker='s',ls='none',color=COLOR[c],markersize=5.5,label=c+' degradation') for c in CATEGORIES]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.55,.063),ncol=3,frameon=False,
               columnspacing=1.1,handletextpad=.4)
    shapes=[Line2D([],[],marker=SHAPE[loc],ls='none',color='#454545',markersize=5.5,label='Locality '+loc) for loc in SHAPE]
    shapes.append(Line2D([],[],color='#333333',lw=1.2,label='Category median'))
    fig.legend(handles=shapes,loc='lower center',bbox_to_anchor=(.55,.038),ncol=3,frameon=False,
               columnspacing=1.5,handletextpad=.5)

def overlap(a,b):
    return max(0,min(a.x1,b.x1)-max(a.x0,b.x0))*max(0,min(a.y1,b.y1)-max(a.y0,b.y0))

def labels(fig,ax,points):
    """Place IDs in display coordinates without changing sample positions."""
    fig.canvas.draw(); renderer=fig.canvas.get_renderer(); bounds=ax.get_window_extent(renderer)
    scale=fig.dpi/72
    xy={s:ax.transData.transform((x,y)) for s,x,y in points}
    boxes=[Bbox.from_extents(x-4*scale,y-4*scale,x+4*scale,y+4*scale) for x,y in xy.values()]
    # Keep ID text away from median bars; otherwise labels appear struck through.
    for line in ax.lines:
        if line.get_gid()=='category_median':
            locations=ax.transData.transform(np.column_stack((line.get_xdata(),line.get_ydata())))
            lo=locations.min(axis=0);hi=locations.max(axis=0)
            boxes.append(Bbox.from_extents(lo[0]-scale,lo[1]-2*scale,hi[0]+scale,hi[1]+2*scale))
    nearest={s:min(np.linalg.norm(c-v) for ss,v in xy.items() if ss!=s) for s,c in xy.items()}
    ordered=sorted(points,key=lambda r:(nearest[r[0]],ORDER.index(r[0])))
    candidates=[(6,6),(-6,6),(6,-6),(-6,-6),(10,0),(-10,0),(0,11),(0,-11),
        (14,10),(-14,10),(14,-10),(-14,-10),(10,16),(-10,16),(10,-16),(-10,-16),
        (23,0),(-23,0),(0,23),(0,-23),(25,15),(-25,15),(25,-15),(-25,-15)]
    placed=[]; records=[]
    for sid,x,y in ordered:
        ann=ax.annotate(sid,(x,y),xytext=(6,6),textcoords='offset points',fontsize=8.5,
            zorder=8,color='#242424',arrowprops={'arrowstyle':'-','color':'#8C8C8C','lw':.45,'shrinkA':1.5,'shrinkB':3})
        ann.set_gid('sample_label_'+sid)
        best=None
        for i,(dx,dy) in enumerate(candidates):
            ann.set_position((dx,dy)); ann.set_ha('left' if dx>0 else 'right' if dx<0 else 'center')
            ann.set_va('bottom' if dy>0 else 'top' if dy<0 else 'center'); ann.update_positions(renderer)
            bb=Text.get_window_extent(ann,renderer).expanded(1.12,1.12)
            outside=sum((max(0,bounds.x0-bb.x0),max(0,bb.x1-bounds.x1),max(0,bounds.y0-bb.y0),max(0,bb.y1-bounds.y1)))
            text_area=sum(overlap(bb,b) for b in placed); point_area=sum(overlap(bb,b) for b in boxes)
            score=10000*outside+200*text_area+100*point_area+abs(dx)+abs(dy)+i/1000
            if best is None or score<best[0]: best=(score,dx,dy,bb,text_area,point_area,outside)
        _,dx,dy,bb,ta,pa,outside=best
        ann.set_position((dx,dy));ann.set_ha('left' if dx>0 else 'right' if dx<0 else 'center')
        ann.set_va('bottom' if dy>0 else 'top' if dy<0 else 'center');ann.update_positions(renderer)
        placed.append(bb);records.append({'sample':sid,'dx_pt':dx,'dy_pt':dy,'label_overlap_px2':ta,'point_overlap_px2':pa,'outside_axes_px':outside})
    return records

def group_points(ax,values,*,technical=None,gid='point'):
    points=[]
    for i,sid in enumerate(ORDER):
        x=i//3+(-.19,0,.19)[i%3];y=float(values[sid])
        if technical:
            low,high=technical[sid]
            ax.plot([x,x],[low,high],color=COLOR[CATEGORY[sid]],lw=1,zorder=3)
        artist=ax.scatter([x],[y],s=37,color=COLOR[CATEGORY[sid]],marker=SHAPE[LOCALITY[sid]],
                          edgecolors='white',linewidths=.5,zorder=5)
        artist.set_gid(gid+'_'+sid);points.append((sid,x,y))
    for group in range(3):
        median=float(np.median([values[s] for s in ORDER[group*3:group*3+3]]))
        line,=ax.plot([group-.30,group+.30],[median,median],color='#333333',lw=1.1,zorder=2)
        line.set_gid('category_median')
    ys=np.asarray(list(values.values()),dtype=float);span=ys.max()-ys.min()
    ax.set_ylim(ys.min()-.22*span,ys.max()+.24*span);ax.set_xlim(-.56,2.56)
    ax.set_xticks(range(3),[c+'\n(n = 3)' for c in CATEGORIES]);ax.tick_params(axis='x',length=0)
    return points

def beta(data):
    fig,axes=axes_grid(240,top=.82,bottom=.17,hspace=.94)
    heading(fig,'Community composition and dispersion','Bray–Curtis ordinations and category dispersion retain marker-specific inputs.')
    rows=[];label_records=[]
    category={(r['Marker'],r['Endpoint']):r for r in data['category_tests']}
    by_its_disp={r['distance']:r for r in data['its_dispersion_statistics']}
    for col,marker in enumerate(('16S','ITS')):
        key='s16' if marker=='16S' else 'its';sid_key='SampleID' if marker=='16S' else 'sample'
        coords={r[sid_key]:r for r in data[key+'_coords']};ax=axes[0,col];letter='ab'[col]
        if marker=='16S':
            st=category[('16S','Bray-Curtis')];disp=data['s16_statistics']['permdisp']
            percent=[float(r['percent_corrected_positive_inertia']) for r in data['s16_axes'][:2]]
            note=f"Rarefied counts; partial R² = {float(st['R2_partial']):.3f}\nP = {float(st['P_raw']):.3f}; Holm₃ P = {float(st['P_Holm']):.3f}"
            dispersion_note=f"PERMDISP: F₂,₆ = {disp['F']:.3f}\nUnadjusted P = {disp['P_raw']:.3f}"
        else:
            st=category[('ITS','Bray-Curtis proportions')];disp=by_its_disp['bray_curtis_proportions']
            percent=[100*float(r['fraction_positive_inertia']) for r in data['its_axes'][:2]]
            note=f"Fungal proportions; partial R² = {float(st['R2_partial']):.3f}\nP = {float(st['P_raw']):.3f}; Holm₄ P = {float(st['P_Holm']):.3f}"
            dispersion_note=f"PERMDISP: F₂,₆ = {float(disp['F']):.3f}\nHolm₂ P = {float(disp['P_Holm_two_distances']):.3f}"
        panel(ax,letter,marker+' | Bray–Curtis',note)
        points=[]
        for sid in ORDER:
            x,y=float(coords[sid]['PCoA1']),float(coords[sid]['PCoA2'])
            artist=ax.scatter([x],[y],s=42,color=COLOR[CATEGORY[sid]],marker=SHAPE[LOCALITY[sid]],
                edgecolors='white',linewidths=.5,zorder=5);artist.set_gid('beta_'+letter+'_'+sid)
            points.append((sid,x,y))
            rows.append({'panel':letter,'marker':marker,'sample':sid,'category':CATEGORY[sid],
                'locality':LOCALITY[sid],'quantity':'PCoA coordinates','display_x':x,'value':y,
                'axis_1_percent_positive_inertia':percent[0],'axis_2_percent_positive_inertia':percent[1],
                'source':PATHS[key+'_coords']})
        xs=np.array([x for _,x,y in points]);ys=np.array([y for _,x,y in points])
        ax.set_xlim(xs.min()-.22*np.ptp(xs),xs.max()+.24*np.ptp(xs))
        ax.set_ylim(ys.min()-.22*np.ptp(ys),ys.max()+.24*np.ptp(ys))
        ax.set_aspect('equal',adjustable='datalim')
        ax.axhline(0,color='#CCCCCC',lw=.6,zorder=1);ax.axvline(0,color='#CCCCCC',lw=.6,zorder=1)
        ax.set_xlabel(f'PCoA 1 ({percent[0]:.2f}%)');ax.set_ylabel(f'PCoA 2 ({percent[1]:.2f}%)')
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4));ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        this_labels=labels(fig,ax,points)
        if marker=='16S':
            for sample,dx,dy in [('B1',-10,-6),('M3',10,9)]:
                ann=next(t for t in ax.texts if t.get_text()==sample)
                ann.set_position((dx,dy));ann.set_ha('left' if dx>0 else 'right')
                ann.set_va('bottom' if dy>0 else 'top')
                for r in this_labels:
                    if r['sample']==sample:
                        r.update(dx_pt=dx,dy_pt=dy,manual_reason='Separate nearby B1 and M3 leaders; sample coordinates unchanged.')
        for r in this_labels:label_records.append({'panel':letter,**r})
        ax=axes[1,col];letter='cd'[col];panel(ax,letter,marker+' | Dispersion',dispersion_note)
        field='distance_to_group_spatial_median' if marker=='16S' else 'distance_to_spatial_median'
        values={r[sid_key]:float(r[field]) for r in data[key+'_dispersion']}
        points=group_points(ax,values,gid='beta_'+letter)
        ax.set_ylabel('Distance to spatial median\n(bias-adjusted)')
        for r in labels(fig,ax,points):label_records.append({'panel':letter,**r})
        for sid,x,y in points:
            rows.append({'panel':letter,'marker':marker,'sample':sid,'category':CATEGORY[sid],
                'locality':LOCALITY[sid],'quantity':'bias-adjusted distance to spatial median','display_x':x,'value':y,
                'axis_1_percent_positive_inertia':'','axis_2_percent_positive_inertia':'','source':PATHS[key+'_dispersion']})
    common_legend(fig)
    fig.text(.55,.017,'Category tests condition on locality; dispersion models do not.',ha='center',fontsize=8.5)
    return fig,rows,label_records

def load():
    data={}
    for key,name in PATHS.items():
        p=INP/'beta'/name
        data[key]=json.loads(p.read_text(encoding='utf-8-sig')) if p.suffix=='.json' else list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
    data['category_tests']=list(csv.DictReader((INP/'alpha/CATEGORY_TESTS.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
    return data

if __name__=='__main__':
    data=load();style();fig,rows,offsets=beta(data)
    for t in fig.findobj(match=Text):
        if t.get_text()=='Locality J':t.set_text('J: Jaén')
        if t.get_text()=='Locality SI':t.set_text('SI: San Ignacio')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/'MAIN'/('Figure_03_Beta_Diversity_REVISED.'+ext),dpi=600)
    plt.close(fig)
    with (QA/'BETA_FIGURE_SOURCE.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (QA/'BETA_RENDER_VALIDATION.json').write_text(json.dumps({'status':'PASS','source_rows':len(rows),'geometry':'Frozen coordinates and spatial-median distances','category_family':'Current Holm3 and Holm4; dispersion separate','label_offsets':offsets},indent=2),encoding='utf-8')

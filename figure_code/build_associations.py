"""Render frozen verified tables. No inferential statistics are fitted here."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox

from paths import ROOT, INP, OUT, QA, TAB
SOIL=['pH','Organic matter','Total nitrogen','Total carbon','available phosphorus','conductivity']
METALS=['Al','Cr','Cu','Fe','Mn','Pb','Zn']; ENV=SOIL+METALS
SAMPLES=['B1','B2','B3','M1','M2','M3','A1','A2','A3']
INDICES=['Observed ASVs','Shannon (ln, nats)']
COL={'Low degradation':'#157E86','Moderate degradation':'#C49516','High degradation':'#934862'}
CMAP=LinearSegmentedColormap.from_list('blue_yellow',['#28669A','#F7F7F2','#D0AA20']); NORM=Normalize(-1,1)
SHORT={'Organic matter':'OM','Total nitrogen':'N','Total carbon':'C','available phosphorus':'P','conductivity':'EC'}
TITLE={'Organic matter':'Organic\nmatter','Total nitrogen':'Total\nnitrogen','Total carbon':'Total\ncarbon','available phosphorus':'Available\nphosphorus','conductivity':'Electrical\nconductivity'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
gs=pd.read_csv(TAB/'GENUS_ENVIRONMENT_SAMPLE_DATA.tsv',sep='\t')
ga=pd.read_csv(TAB/'GENUS_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv',sep='\t')
aa=pd.read_csv(TAB/'ALPHA_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv',sep='\t')
aps=pd.read_csv(TAB/'ALPHA_ENVIRONMENT_SAMPLE_DATA.tsv',sep='\t')
assert len(ga)==312 and len(aa)==52 and len(gs)==216 and len(aps)==36
assert all(ga.primary_family_size==156) and all(aa.primary_family_size==26)
assert set(ga.variable)==set(ENV)==set(aa.variable)
register=[]; plotted=[]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(fig,figid,name,caption,sources):
    folder='MAIN' if not figid.startswith('S') else 'SUPPLEMENTARY'
    base=OUT/folder/name
    # Move only annotation text; preserve measured coordinates and identify all units.
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for ax in fig.axes:
        anns=[a for a in ax.texts if isinstance(a,matplotlib.text.Annotation)]
        obstacles=[]
        for ann in anns:
            px,py=ax.transData.transform(ann.xy);obstacles.append(Bbox.from_extents(px-4,py-4,px+4,py+4))
        used=[]
        for ann in anns:
            choices=[(6,6),(-19,6),(6,-13),(-19,-13),(-6,14),(-6,-20),(17,0),(-29,0),(11,20),(-24,20),(11,-24),(-24,-24),(22,10),(-36,10),(22,-14),(-36,-14),(0,27),(0,-29)]
            options=[]
            for n,offset in enumerate(choices):
                ann.set_position(offset);box=ann.get_window_extent(renderer).expanded(1.40,1.25)
                count=sum(box.overlaps(b) for b in obstacles+used)
                outside=not (ax.bbox.x0<=box.x0 and box.x1<=ax.bbox.x1 and ax.bbox.y0<=box.y0 and box.y1<=ax.bbox.y1)
                options.append((100*count+1000*outside+n,offset,box))
            _,offset,box=min(options,key=lambda z:z[0]);ann.set_position(offset);used.append(box)
    for ext in ['pdf','svg','png']:
        kw={'dpi':400} if ext=='png' else {}
        if ext=='pdf':kw['metadata']={'Title':name.replace('_',' '),'Creator':'Matplotlib '+matplotlib.__version__,'CreationDate':None,'ModDate':None}
        fig.savefig(base.with_suffix('.'+ext),**kw)
    cap_source=INP/'FIGURE_CAPTIONS.json'
    if cap_source.exists():
        authoritative=json.loads(cap_source.read_text(encoding='utf-8'))[figid]
        caption='# '+authoritative['heading']+'\n\n'+authoritative['body']
    base.with_suffix('.caption.md').write_text(caption.strip()+'\n',encoding='utf-8')
    register.append({'figure_id':figid,'figure':name,'scope':folder,'status':'RENDERED_FROM_CURRENT_FROZEN_NUMERICAL_TABLES','source_id':'STATISTICAL_EXTENSION_20261008','source_evidence':'; '.join(sources),'caption_path':str(base.with_suffix('.caption.md').relative_to(OUT)), 'files':[{'path':str(base.with_suffix('.'+ext).relative_to(OUT)),'sha256':sha(base.with_suffix('.'+ext))} for ext in ['pdf','svg','png']]})
    plt.close(fig)

def gens(marker):return gs.loc[gs.marker==marker].sort_values('display_order').genus.drop_duplicates().tolist()
def assoc(marker,g,v):return ga.loc[(ga.marker==marker)&(ga.genus==g)&(ga.variable==v)].iloc[0]
def alpha(marker,index,v):return aa.loc[(aa.marker==marker)&(aa.alpha_index==index)&(aa.variable==v)].iloc[0]
def legend(fig,y):
    h=[Line2D([],[],marker='o',ls='',color=COL[k],markersize=5,label=k.replace(' degradation','')) for k in COL]
    h += [Line2D([],[],marker=m,ls='',color='#414141',markersize=5,label=l) for m,l in [('o','J (Jaén)'),('^','SI (San Ignacio)')]]
    fig.legend(handles=h,loc='lower center',bbox_to_anchor=(.53,y),ncol=5,frameon=False,fontsize=8.5,handletextpad=.5,columnspacing=1.1)

def gid_for(kind,marker,g,v):return kind+'_'+marker+'_'+g.replace(' ','_')+'_'+v.replace(' ','_')
def track(figid,gid,kind,source,x=None,y=None):
    plotted.append({'figure_id':figid,'svg_id':gid,'kind':kind,'source':source,'x':x,'y':y})

common_genus=('The fixed 12 nominal genus-rank lineages per marker match Supplementary Figure S25 and GENUS_SELECTION_CURRENT.tsv. '
 'Coefficients correlate globally ranked full-table CLR values and globally ranked recorded environmental variables after removing locality means. '
 'CLR uses all 1143 complete lineage groups for 16S or 1185 for ITS with a 0.5-count pseudocount. '
 'Two-sided P values enumerate all 4320 within-locality arrangements including identity once; BH adjustment uses 156 tests per marker (12 genera x 13 variables). '
 'This post hoc scope extension adds total nitrogen and total carbon. No association has BH156-adjusted P < 0.05. '
 'Conditioning is on locality, not jointly on degradation, vegetation or other environmental variables. Exchangeability within locality remains an assumption. '
 'These exploratory estimates establish neither causal effects, absolute abundance, contamination, tolerance nor bioremediation. '
 'J denotes Jaén and SI San Ignacio. Units and laboratory methods remain unconfirmed. BY156, pseudocount and leave-one-out diagnostics are reported in the current tables. '
 'Sensitivity ranges are not confidence intervals. Genus labels are reference annotations, not validated species identifications.')

def abundance_panel(ax,marker,figid):
    names=gens(marker)
    for i,g in enumerate(names):
        rr=gs.loc[(gs.marker==marker)&(gs.genus==g)].set_index('sample_id').loc[SAMPLES]
        mean=rr.relative_abundance_percent.mean()
        ax.barh(i,mean,height=.54,color='#CFDADF',edgecolor='#AABBC4',linewidth=.6,zorder=2)
        for j,(s,r) in enumerate(rr.iterrows()):
            y=i+(j-4)*.035
            p=ax.scatter(r.relative_abundance_percent,y,s=13 if r.locality=='J' else 17,marker='o' if r.locality=='J' else '^',c=COL[r.degradation],edgecolors='white',linewidths=.35,zorder=4)
            gid=gid_for('abundance',marker,g,s);p.set_gid(gid)
            track(figid,gid,'abundance',{'marker':marker,'genus':g,'sample_id':s},float(r.relative_abundance_percent),float(y))
    ax.set_yticks(range(12),[g.replace('_',' ') for g in names],fontstyle='italic',fontsize=9)
    ax.set_xlim(left=0);ax.xaxis.set_major_locator(MaxNLocator(3));ax.set_xlabel('Relative abundance\n(% of retained library)',fontsize=8)
    ax.set_title('Abundance',fontsize=10,pad=10);ax.spines['left'].set_visible(False);ax.grid(axis='x',color='#E5E9EB',linewidth=.55,zorder=0)

def lollipop(marker,variables,figid):
    name=f'Figure_{figid}_{marker}_'+('Genus_Metal_Associations' if len(variables)==7 else 'Soil_Associations_BH156')
    title=f'{marker} genus abundance and '+('soil-metal' if len(variables)==7 else 'soil-property')+' associations'
    fig=plt.figure(figsize=(12.8,8.4)); grid=fig.add_gridspec(1,len(variables)+1,left=.205,right=.985,bottom=.205,top=.855,width_ratios=[1.5]+[1]*len(variables),wspace=.3)
    axes=[fig.add_subplot(grid[0,i]) for i in range(len(variables)+1)]
    for ax in axes:
        ax.set_ylim(11.65,-.7);ax.tick_params(axis='y',length=0)
        for k in range(0,12,2):ax.axhspan(k-.48,k+.48,color='#F1F4F5',zorder=0)
        ax.set_yticks(range(12),[])
    abundance_panel(axes[0],marker,figid)
    for ax,v in zip(axes[1:],variables):
        ax.spines['left'].set_visible(False);ax.axvline(0,color='#AAB2B8',ls=':',linewidth=.9)
        for i,g in enumerate(gens(marker)):
            r=assoc(marker,g,v);x=float(r.rho_locality_adjusted)
            ax.hlines(i-.15,0,x,color='#909FA7',linewidth=1.1,zorder=2)
            p=ax.scatter(x,i-.15,s=27,c=[CMAP(NORM(x))],edgecolors='#637B89',linewidths=.65,zorder=3)
            gid=gid_for('rho',marker,g,v);p.set_gid(gid);track(figid,gid,'rho',{'marker':marker,'genus':g,'variable':v},x,float(i-.15))
            txt=ax.text(0,i+.23,f'P = {r.q_BH156:.3f}',ha='center',va='center',fontsize=8.5,color='#243740')
            txt.set_gid(gid+'_adjusted_p')
            track(figid,gid+'_adjusted_p','p_annotation',{'marker':marker,'genus':g,'variable':v,'value':float(r.q_BH156),'field':'q_BH156'},0,float(i+.23))
        ax.set_xlim(-1.12,1.12);ax.set_xticks([-1,0,1]);ax.tick_params(axis='x',labelsize=8)
        ax.set_title(TITLE.get(v,v),fontsize=10,pad=10)
    fig.text(.035,.958,title,fontsize=12,fontweight='bold');fig.text(.035,.915,'Twelve fixed genera; nine soil units. Numbers below estimates: BH-adjusted P (156 tests per marker).',fontsize=9,color='#4E626B')
    fig.text(.62,.133,'Locality-adjusted rank correlation',ha='center',fontsize=9);legend(fig,.068)
    fig.text(.035,.025,'Bars: sample means. Points: individual units. All adjusted P values exceed 0.05; raw P values are in the source table.',fontsize=8.8)
    cap=title+'. Bars are unweighted mean relative abundances across nine complete retained libraries; all nine sample proportions are shown without renormalizing the selected genera. Lollipops display all '+('84 metal associations' if len(variables)==7 else '72 soil-property associations')+' in this panel. '+common_genus
    save(fig,figid,name,cap,['GENUS_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv','GENUS_ENVIRONMENT_SAMPLE_DATA.tsv'])

for marker,a,b in [('16S','05','S16'),('ITS','06','S17')]:
    lollipop(marker,METALS,a);lollipop(marker,SOIL,b)

def matrix(ax,marker,variables,figid,numbers=False):
    names=gens(marker)
    for i,g in enumerate(names):
        for j,v in enumerate(variables):
            r=assoc(marker,g,v);x=float(r.rho_locality_adjusted)
            rect=Rectangle((j-.5,i-.5),1,1,facecolor=CMAP(NORM(x)),edgecolor='white',lw=.4)
            gid=gid_for('matrix',marker,g,v);rect.set_gid(gid);ax.add_patch(rect);track(figid,gid,'matrix',{'marker':marker,'genus':g,'variable':v},x,None)
            if numbers:
                label=f'{x:.2f}\nP = {r.q_BH156:.3f}' if numbers is True else f'{r.q_BH156:.3f}'
                txt=ax.text(j,i,label,ha='center',va='center',fontsize=8 if numbers is True else 7.6,color='white' if x<-.7 else '#22303A')
                txt.set_gid(gid+'_adjusted_p')
                track(figid,gid+'_adjusted_p','p_annotation',{'marker':marker,'genus':g,'variable':v,'value':float(r.q_BH156),'field':'q_BH156'},j,i)
    ax.set_xlim(-.5,len(variables)-.5);ax.set_ylim(11.5,-.5)
    ax.set_yticks(range(12),[g.replace('_',' ') for g in names],fontstyle='italic',fontsize=9)
    ax.set_xticks(range(len(variables)),[SHORT.get(v,v) for v in variables],fontsize=9);ax.tick_params(length=0)
    ax.set_title(marker,loc='left',fontweight='bold',fontsize=11)

for marker,figid in [('16S','S05'),('ITS','S06')]:
    fig,ax=plt.subplots(figsize=(9,8.3));fig.subplots_adjust(left=.28,right=.80,top=.86,bottom=.19)
    matrix(ax,marker,SOIL,figid,True)
    cb=fig.add_axes([.84,.29,.018,.48]);fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cb,label='Locality-adjusted rank correlation')
    fig.text(.035,.958,f'{marker} genus associations with six soil properties',fontsize=12,fontweight='bold')
    fig.text(.035,.112,'OM: organic matter. N: total nitrogen. C: total carbon.',fontsize=9)
    fig.text(.035,.075,'P: available phosphorus. EC: electrical conductivity.',fontsize=9)
    fig.text(.035,.036,'Cells: coefficient and BH156-adjusted P. All adjusted P values exceed 0.05.',fontsize=9)
    save(fig,figid,f'Figure_{figid}_Environment_Heatmap_{marker}',f'{marker} soil-property correlation matrix. Every one of the 72 coefficients from Figure '+('S16' if marker=='16S' else 'S17')+' is shown numerically; rows retain the fixed abundance order and columns retain the prespecified display order. The common scale is -1 to +1. No additional testing or significance selection is performed. '+common_genus,['GENUS_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv'])

fig,axes=plt.subplots(1,2,figsize=(17,8.5));fig.subplots_adjust(left=.135,right=.90,top=.82,bottom=.20,wspace=.64)
for ax,marker in zip(axes,['16S','ITS']):matrix(ax,marker,ENV,'S21',numbers='p_only')
cb=fig.add_axes([.925,.27,.012,.48]);fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cb,label='Locality-adjusted rank correlation')
fig.text(.03,.945,'Genus associations across the complete environmental family',fontsize=13,fontweight='bold')
fig.text(.03,.89,'Cell numbers: BH156-adjusted P. Colours: locality-adjusted rank correlation. 12 genera x 13 variables per marker.',fontsize=10)
fig.text(.03,.125,'OM: organic matter. N: total nitrogen. C: total carbon. P: available phosphorus. EC: electrical conductivity.',fontsize=9)
fig.text(.03,.08,'Rows follow mean abundance. Column order is identical in both markers. No association meets BH156-adjusted P < 0.05.',fontsize=9)
fig.text(.03,.035,'Raw P values and coefficients are retained in the source table. Cell numbers are adjusted probabilities, not coefficients.',fontsize=9)
save(fig,'S21','Figure_S21_All_Environmental_Associations','Complete genus-environment matrices. All 312 coefficients are shown, with 156 associations in each marker. Rows preserve the fixed abundance order, and columns include six soil properties and seven metals in the same order in both panels. These are the same estimates shown in Figures 5-6 and S16-S17; no additional test is performed. '+common_genus,['GENUS_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv'])

common_alpha=('The family comprises two alpha indices x 13 environmental variables (26 tests) separately for each marker. '
 'Coefficients correlate globally ranked alpha indices and globally ranked environmental values after removing locality means. '
 'Two-sided conditional P values enumerate all 4320 within-locality arrangements including identity once; BH adjustment uses all 26 tests per marker. '
 'Neither marker has a BH26-adjusted P < 0.05 association. The 16S values use one rarefaction realization at 174661 sequences; ITS values are means of 100 technical rarefaction draws at 450343 fungal sequences. '
 'These estimates use nine soil units, not technical draws as biological replicates. Category, vegetation and other environmental variables are not jointly controlled. '
 'The extension is exploratory and post hoc. Within-locality exchangeability remains assumed. No causal relationship or absence of an ecological effect is established. '
 'J denotes Jaén and SI San Ignacio. Laboratory units remain unconfirmed.')
fig,axes=plt.subplots(1,2,figsize=(9.8,7.8));fig.subplots_adjust(left=.23,right=.875,top=.85,bottom=.17,wspace=.75)
for ax,marker in zip(axes,['16S','ITS']):
    for i,v in enumerate(ENV):
        for j,index in enumerate(INDICES):
            r=alpha(marker,index,v);x=float(r.rho_locality_adjusted)
            rect=Rectangle((j-.5,i-.5),1,1,facecolor=CMAP(NORM(x)),edgecolor='white',lw=.6)
            gid=gid_for('alpha',marker,index,v);rect.set_gid(gid);ax.add_patch(rect);track('S03',gid,'alpha_matrix',{'marker':marker,'alpha_index':index,'variable':v},x,None)
            ax.text(j,i,f'{x:.2f}\n({r.q_BH26:.3f})',ha='center',va='center',fontsize=8.5,color='white' if x<-.7 else '#22303A')
    ax.set_xlim(-.5,1.5);ax.set_ylim(12.5,-.5);ax.set_xticks([0,1],['Observed ASVs','Shannon'],fontsize=9)
    ax.set_yticks(range(13),['Available phosphorus' if v=='available phosphorus' else 'Conductivity' if v=='conductivity' else v for v in ENV],fontsize=9)
    ax.tick_params(length=0);ax.set_title(marker,loc='left',fontweight='bold',fontsize=12)
cb=fig.add_axes([.915,.30,.016,.43]);fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cb,label='Locality-adjusted rank correlation')
fig.text(.035,.95,'Alpha-diversity associations with all environmental variables',fontsize=13,fontweight='bold')
fig.text(.035,.9,'Cells: coefficient, with BH26-adjusted P in parentheses',fontsize=10)
fig.text(.035,.10,'26 associations per marker. No association meets BH26-adjusted P < 0.05.',fontsize=10)
fig.text(.035,.052,'ITS uses technical rarefaction means; the analytical sample size remains nine units.',fontsize=9)
save(fig,'S03','Figure_S03_Alpha_Environment_Correlations','Exploratory alpha-diversity associations for both markers. All 52 coefficients and BH-adjusted P values are displayed, with 26 associations per marker. Parentheses show BH26-adjusted P, not uncertainty intervals. The common colour scale is -1 to +1. '+common_alpha,['ALPHA_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv'])

for marker,figid in [('16S','S04'),('ITS','S23')]:
    fig,axes=plt.subplots(6,2,figsize=(10,13.4));fig.subplots_adjust(left=.12,right=.96,top=.895,bottom=.105,wspace=.32,hspace=.88)
    for i,v in enumerate(SOIL):
        for j,index in enumerate(INDICES):
            ax=axes[i,j]; rr=aps.loc[(aps.marker==marker)&(aps.alpha_index==index)].set_index('sample_id').loc[SAMPLES]
            a=alpha(marker,index,v)
            for k,(s,r) in enumerate(rr.iterrows()):
                p=ax.scatter(r.alpha_value,r[v],s=28,marker='o' if r.locality=='J' else '^',c=COL[r.degradation],edgecolors='white',linewidth=.45,zorder=3)
                gid=gid_for('raw',marker,index,v)+'_'+s;p.set_gid(gid);track(figid,gid,'alpha_scatter',{'marker':marker,'alpha_index':index,'variable':v,'sample_id':s},float(r.alpha_value),float(r[v]))
                ax.annotate(s,(r.alpha_value,r[v]),xytext=(4,5 if k%2==0 else -9),textcoords='offset points',fontsize=7.5)
            ax.set_title(f'r = {a.rho_locality_adjusted:.3f}; BH26 P = {a.q_BH26:.3f}',loc='left',fontsize=9,pad=7)
            ax.set_ylabel(('Available phosphorus' if v=='available phosphorus' else 'Conductivity' if v=='conductivity' else v)+'\n'+('' if v=='pH' else '(recorded value)'),fontsize=9)
            ax.set_xlabel('Observed ASVs' if j==0 else 'Shannon (ln, nats)',fontsize=9)
            ax.margins(x=.15,y=.23);ax.xaxis.set_major_locator(MaxNLocator(4));ax.yaxis.set_major_locator(MaxNLocator(4));ax.grid(color='#E7ECEF',linewidth=.5);ax.set_axisbelow(True)
    fig.text(.07,.967,f'{marker} alpha diversity and six recorded soil properties',fontsize=14,fontweight='bold')
    fig.text(.07,.933,'Raw sample coordinates; annotations use locality-adjusted rank associations',fontsize=10)
    legend(fig,.048)
    fig.text(.07,.029,'No fitted lines: raw coordinates do not represent the locality-residualized coefficients.',fontsize=9)
    fig.text(.07,.013,'These 12 associations belong to the full 26-test marker-specific family shown in Figure S03.',fontsize=9)
    save(fig,figid,f'Figure_{figid}_{marker}_Alpha_Soil_Scatterplots',f'Raw observations underlying 12 {marker} alpha-diversity associations with six recorded soil properties. Every panel shows all nine labelled units. The displayed coefficients are locality-adjusted rank correlations; they are not slopes through the raw plotted coordinates. No regression lines or confidence bands are shown. These panels are a subset of the full 26-test family; all metal associations are retained in Figure S03 and the source table. '+common_alpha,['ALPHA_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv','ALPHA_ENVIRONMENT_SAMPLE_DATA.tsv'])

# Descriptive six-variable shared-information matrix: ranks then locality means.
env=aps.loc[(aps.marker=='16S')&(aps.alpha_index==INDICES[0])].set_index('sample_id').loc[SAMPLES]
ranks=env[SOIL].rank(method='average')
for locality in env.locality.unique():ranks.loc[env.locality==locality]-=ranks.loc[env.locality==locality].mean(axis=0)
corr=ranks.corr()
fig,ax=plt.subplots(figsize=(7.1,6.2));fig.subplots_adjust(left=.22,right=.85,top=.84,bottom=.20)
for i,v in enumerate(SOIL):
    for j,w in enumerate(SOIL):
        x=float(corr.loc[v,w]);rect=Rectangle((j-.5,i-.5),1,1,facecolor=CMAP(NORM(x)),edgecolor='white',lw=.5);ax.add_patch(rect)
        gid='environment_'+str(i)+'_'+str(j);rect.set_gid(gid);track('S07',gid,'descriptive_correlation',{'variable_x':v,'variable_y':w},x,None)
        ax.text(j,i,f'{x:.2f}',ha='center',va='center',color='white' if x<-.7 else '#22303A',fontsize=9)
ax.set_xlim(-.5,5.5);ax.set_ylim(5.5,-.5);ax.set_xticks(range(6),[SHORT.get(v,v) for v in SOIL]);ax.set_yticks(range(6),['pH','Organic matter','Total nitrogen','Total carbon','Available phosphorus','Conductivity']);ax.tick_params(length=0)
cb=fig.add_axes([.89,.30,.02,.44]);fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cb)
fig.text(.04,.944,'Shared variation among six recorded soil properties',fontsize=12,fontweight='bold')
fig.text(.04,.122,'Locality-adjusted rank correlations; nine soil units. No additional tests.',fontsize=9)
fig.text(.04,.078,'N: total nitrogen. C: total carbon. OM: organic matter.',fontsize=9)
fig.text(.04,.034,'P: available phosphorus. EC: electrical conductivity.',fontsize=9)
save(fig,'S07','Figure_S07_Soil_Variable_Correlations','Descriptive correlations among the six displayed soil properties, including total nitrogen and total carbon. Global midranks were residualized on intercept and locality and correlated across the same nine samples. Diagonal cells are self-correlations, and symmetric cells duplicate the same coefficient. No hypothesis test is added. Shared variation, including the strong OM/carbon relationship, limits interpreting individual associations as independent environmental effects. Coefficients are invariant to positive unit changes; they do not validate laboratory units or causal direction. J denotes Jaén and SI San Ignacio.',['ALPHA_ENVIRONMENT_SAMPLE_DATA.tsv','ENVIRONMENT_ALL13_CORRELATIONS_DESCRIPTIVE.tsv'])

(QA/'current_figure_register.json').write_text(json.dumps(register,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
(QA/'PLOTTED_DATA_REGISTER.json').write_text(json.dumps(plotted,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps({'created_figures':len(register),'recorded_graphical_data_objects':len(plotted)}))

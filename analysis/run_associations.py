"""Current Yanayacu genus/alpha-environment analyses (Python + NumPy).

Portable derivative of the verified source implementation. Historical-only
comparison tables and freeze commands are omitted; statistical calculations
and current output schema are unchanged and checked against released hashes.
"""
from pathlib import Path
import argparse, collections, csv, datetime, hashlib, itertools, json, platform, shutil, sys
import numpy as np

SAMPLES = ['B1','B2','B3','M1','M2','M3','A1','A2','A3']
MARKERS = ['16S','ITS']
VARIABLES = ['pH','Organic matter','Total nitrogen','Total carbon','available phosphorus','conductivity','Al','Cr','Cu','Fe','Mn','Pb','Zn']
METALS = ['Al','Cr','Cu','Fe','Mn','Pb','Zn']
FIELDS = {v: v+' (unit unconfirmed)' for v in VARIABLES}
FIELDS.update({'pH':'pH (dimensionless)','available phosphorus':'Available phosphorus (unit unconfirmed)','conductivity':'Conductivity (unit unconfirmed)'})
INDICES = ['Observed ASVs','Shannon (ln, nats)']

def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def dump(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def ranks(v):
    v=np.asarray(v); order=np.argsort(v,kind='stable'); out=np.empty(len(v),float); i=0
    while i<len(v):
        j=i+1
        while j<len(v) and v[order[i]]==v[order[j]]: j+=1
        out[order[i:j]]=(i+j+1)/2; i=j
    return out
def residualizer(loc):
    x=np.array([[1.,float(a=='SI')] for a in loc]); return np.eye(len(loc))-x@np.linalg.pinv(x)
def direction(x): return int(x>1e-12)-int(x< -1e-12)
def adjust(p,method):
    p=np.asarray(p,float); order=np.argsort(p,kind='stable'); m=len(p)
    factor=m*(sum(1/k for k in range(1,m+1)) if method=='BY' else 1.)
    out=np.empty(m); out[order]=np.minimum(1.,np.minimum.accumulate((p[order]*factor/np.arange(1,m+1))[::-1])[::-1]); return out
def perms(loc):
    blocks=[np.flatnonzero(np.asarray(loc)==v) for v in sorted(set(loc))]; out=[]
    for parts in itertools.product(*[list(itertools.permutations(b.tolist())) for b in blocks]):
        row=np.arange(len(loc))
        for block,part in zip(blocks,parts): row[block]=part
        out.append(row)
    out=np.asarray(out)
    assert len(out)==len(set(map(tuple,out)))==4320
    assert sum(np.array_equal(p,np.arange(9)) for p in out)==1
    assert all(np.array_equal(np.asarray(loc)[p],loc) for p in out)
    return out
def coefficient(a,b,loc):
    h=residualizer(loc); x=h@ranks(a); y=h@ranks(b); norm=np.linalg.norm(x)*np.linalg.norm(y)
    if norm<1e-12: raise ValueError('Undefined coefficient: zero residual variance')
    return float(x@y/norm)

def analyse(package, output):
    started=utc(); inp=package/'data/analysis_inputs'; results=output/'validation'; tables=output/'tables'
    receiptpath=results/'ASSOCIATION_EXTENSION_RECEIPT.json'
    if receiptpath.exists(): raise FileExistsError('Refusing to overwrite completed analysis; use a fresh replay release')
    planpath=inp/'DESIGN.json'; plan=json.loads(planpath.read_text()); assert plan['script_sha256']==sha(__file__)
    for item in plan['sources']: assert sha(inp/item['file'])==item['sha256']
    md={r['Sample ID']:r for r in read(inp/'S16_Study_metadata.tsv')}; loc=np.array([md[s]['Locality'] for s in SAMPLES]); mat=residualizer(loc); pi=perms(loc)
    env={v:np.array([float(md[s][FIELDS[v]]) for s in SAMPLES]) for v in VARIABLES}
    resid={v:mat@ranks(x) for v,x in env.items()}
    archived=read(inp/'GENUS_ENVIRONMENT_SAMPLE_DATA.tsv'); selections=read(inp/'GENUS_SELECTION.tsv')
    alpharows=read(inp/'Alpha_Primary_Estimates_EN.tsv'); assert len(alpharows)==18
    genus_tests=[]; alpha_tests=[]; gsamples=[]; asamples=[]; loorows=[]; clr_errors=[]; denominators=[]
    def fit(values,var,identity):
        x=resid[var]; y=mat@ranks(values); norm=np.linalg.norm(x)*np.linalg.norm(y); assert norm>1e-12
        rho=float(x@y/norm); null=x[pi]@y/norm; extreme=int(np.count_nonzero(np.abs(null)>=abs(rho)-1e-12)); loo=[]
        for i,sid in enumerate(SAMPLES):
            keep=np.arange(9)!=i; r=coefficient(env[var][keep],values[keep],loc[keep]); loo.append(r)
            loorows.append(dict(**identity,variable=var,omitted_sample=sid,n=8,rho_locality_adjusted=r))
        return dict(n=9,rho_locality_adjusted=rho,extreme_permutations=extreme,total_permutations=4320,p_two_sided=extreme/4320,
          loo_rho_min=min(loo),loo_rho_max=max(loo),loo_same_sign_count=sum(direction(r)==direction(rho) for r in loo),loo_total=9)
    for marker in MARKERS:
        full=read(inp/(marker+'_ALL_GENUS_RANK_LINEAGE_COUNTS_DERIVED.tsv')); keys=[r['feature_id'] for r in full]
        counts=np.array([[int(r[s]) for s in SAMPLES] for r in full]); denom=counts.sum(axis=0)
        assert len(keys)==len(set(keys))=={'16S':1143,'ITS':1185}[marker]
        assert int(counts.sum())=={'16S':1752112,'ITS':5388961}[marker]
        clr=np.log(counts+.5); clr-=clr.mean(axis=0); clr1=np.log(counts+1.); clr1-=clr1.mean(axis=0)
        chosen=[r for r in selections if r['marker']==marker]; assert len(chosen)==12
        eligible=[(r['feature_id'],float(np.mean(counts[i]/denom))) for i,r in enumerate(full) if r['nominal_genus'].lower()=='true' and np.count_nonzero(counts[i])>=6 and np.mean(counts[i]/denom)>=.001]
        eligible.sort(key=lambda t:(-t[1],t[0])); assert [r['feature_id'] for r in chosen]==[r[0] for r in eligible[:12]]
        for j,sid in enumerate(SAMPLES): denominators.append(dict(marker=marker,sample_id=sid,retained_library_count=int(denom[j]),full_reference_lineages=len(keys)))
        tests=[]
        for selected in chosen:
            key=selected['feature_id']; idx=keys.index(key); source={r['sample_id']:r for r in archived if r['marker']==marker and r['feature_id']==key}
            for j,sid in enumerate(SAMPLES):
                r=source[sid]; assert r['locality']==loc[j] and r['degradation']==md[sid]['Degradation']
                clr_errors.extend([abs(float(r['clr_pc0p5'])-clr[idx,j]),abs(float(r['clr_pc1'])-clr1[idx,j]),abs(float(r['relative_abundance'])-counts[idx,j]/denom[j]),abs(int(r['count'])-counts[idx,j])])
                gsamples.append(dict(marker=marker,display_order=selected['display_order'],feature_id=key,genus=selected['genus'],sample_id=sid,locality=loc[j],degradation=md[sid]['Degradation'],count=int(counts[idx,j]),retained_library_count=int(denom[j]),relative_abundance=float(counts[idx,j]/denom[j]),relative_abundance_percent=float(100*counts[idx,j]/denom[j]),clr_pc0p5=float(clr[idx,j]),clr_pc1=float(clr1[idx,j]),**{v:float(env[v][j]) for v in VARIABLES}))
            for var in VARIABLES:
                identity=dict(analysis_type='genus',marker=marker,response_id=key,response_label=selected['genus']); fitrow=fit(clr[idx],var,identity)
                rho1=coefficient(env[var],clr1[idx],loc)
                row=dict(marker=marker,display_order=selected['display_order'],feature_id=key,genus=selected['genus'],variable=var,variable_type='metal' if var in METALS else 'soil',
                  prevalence=int(selected['prevalence']),mean_relative_abundance=float(selected['mean_relative_abundance']),mean_relative_abundance_percent=float(selected['mean_relative_abundance_percent']),**fitrow,
                  q_BH156=0.,q_BY156=0.,primary_family_size=156,rho_pseudocount1=rho1,delta_rho_pseudocount=abs(fitrow['rho_locality_adjusted']-rho1),pseudocount_same_sign=direction(fitrow['rho_locality_adjusted'])==direction(rho1))
                tests.append(row)
        assert len(tests)==156
        for method in ['BH','BY']:
            for row,q in zip(tests,adjust([r['p_two_sided'] for r in tests],method)): row['q_'+method+'156']=float(q)
        genus_tests.extend(tests)
        ar={r['Sample']:r for r in alpharows if r['Marker']==marker}; assert set(ar)==set(SAMPLES); atests=[]
        for index in INDICES:
            vals=np.array([float(ar[s][index]) for s in SAMPLES]); assert np.isfinite(vals).all()
            for j,sid in enumerate(SAMPLES):
                assert ar[sid]['Locality']==loc[j] and ar[sid]['Category']==md[sid]['Degradation']
                asamples.append(dict(marker=marker,alpha_index=index,sample_id=sid,locality=loc[j],degradation=md[sid]['Degradation'],alpha_value=float(vals[j]),estimator=ar[sid]['Estimator'],draws=int(ar[sid]['Draws for estimator']),depth=int(ar[sid]['Depth (sequences)']),**{v:float(env[v][j]) for v in VARIABLES}))
            for var in VARIABLES:
                fitrow=fit(vals,var,dict(analysis_type='alpha',marker=marker,response_id=index,response_label=index))
                row=dict(marker=marker,alpha_index=index,variable=var,**fitrow,q_BH26=0.,q_BY26=0.,primary_family_size=26,estimator=ar[SAMPLES[0]]['Estimator'],depth=int(ar[SAMPLES[0]]['Depth (sequences)']))
                atests.append(row)
        assert len(atests)==26
        for method in ['BH','BY']:
            for row,q in zip(atests,adjust([r['p_two_sided'] for r in atests],method)): row['q_'+method+'26']=float(q)
        alpha_tests.extend(atests)
    assert max(clr_errors)<1e-10
    correlations=[dict(variable_a=a,variable_b=b,n=9,spearman_rho=float(np.corrcoef(ranks(env[a]),ranks(env[b]))[0,1]),rho_locality_adjusted=coefficient(env[a],env[b],loc),interpretation='Descriptive dependence only; no P value; does not establish laboratory derivation') for a,b in itertools.combinations(VARIABLES,2)]
    envsamples=[dict(sample_id=sid,locality=loc[j],degradation=md[sid]['Degradation'],**{v:float(env[v][j]) for v in VARIABLES}) for j,sid in enumerate(SAMPLES)]
    outputs={'GENUS_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv':genus_tests,'ALPHA_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv':alpha_tests,'GENUS_ENVIRONMENT_SAMPLE_DATA.tsv':gsamples,'ALPHA_ENVIRONMENT_SAMPLE_DATA.tsv':asamples,'ENVIRONMENT_ALL13_CORRELATIONS_DESCRIPTIVE.tsv':correlations,'ENVIRONMENT_ALL13_SAMPLE_VALUES.tsv':envsamples,'ASSOCIATION_LEAVE_ONE_OUT.tsv':loorows,'GENUS_SELECTION_CURRENT.tsv':selections,'LIBRARY_DENOMINATORS_CURRENT.tsv':denominators}
    for name,rows in outputs.items(): write(tables/name,rows)
    summary=[]
    for kind,rows,qfield,byfield in [('genus',genus_tests,'q_BH156','q_BY156'),('alpha',alpha_tests,'q_BH26','q_BY26')]:
        for marker in MARKERS:
            subset=[r for r in rows if r['marker']==marker]
            summary.append(dict(analysis_type=kind,marker=marker,tests=len(subset),min_raw_p=min(r['p_two_sided'] for r in subset),min_BH=min(r[qfield] for r in subset),min_BY=min(r[byfield] for r in subset),BH_below_0p05=sum(r[qfield]<.05 for r in subset),BY_below_0p05=sum(r[byfield]<.05 for r in subset),any_LOO_sign_change=sum(r['loo_same_sign_count']<9 for r in subset)))
    write(tables/'ASSOCIATION_FAMILY_SUMMARY.tsv',summary)
    receipt=dict(status='NUMERICAL_EXECUTION_PASS',started_utc=started,completed_utc=utc(),design_sha256=sha(planpath),script_sha256=sha(__file__),command=['python','analysis/run_associations.py','--package','<package>','--output','<new-output>'],environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform()),checks=dict(n=9,J=6,SI=3,genus_tests=312,alpha_tests=52,max_full_CLR_error=max(clr_errors),permutations=4320,identity_once=True,leave_one_out_coefficients=len(loorows)),summaries=summary,outputs=[dict(path=str((tables/name).relative_to(output)).replace('\\','/'),rows=len(rows),sha256=sha(tables/name)) for name,rows in outputs.items()])
    receipt['outputs'].append(dict(path='tables/ASSOCIATION_FAMILY_SUMMARY.tsv',rows=len(summary),sha256=sha(tables/'ASSOCIATION_FAMILY_SUMMARY.tsv')))
    dump(receiptpath,receipt); print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); package=a.package.resolve(); output=a.output.resolve()
    if output.exists(): raise FileExistsError('Choose a new output directory; no overwrites are permitted.')
    if output==package or package in output.parents: raise ValueError('Write replay outside the released package.')
    analyse(package,output)

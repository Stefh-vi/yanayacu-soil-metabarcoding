"""Independent arithmetic QA; does not import the production analysis module.

Uses pairwise-count ranks, block-mean subtraction, direct product null statistics,
and ordered definitions of adjustments. Category models use 180 distinct label
allocations and least-squares fits, unlike the production residual projectors.
"""
from pathlib import Path
import argparse, collections, csv, datetime, hashlib, itertools, json, math, sys
import numpy as np

def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_rows(path, rows):
    with path.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
def rank(v): return np.array([sum(y<x for y in v)+(sum(y==x for y in v)+1)/2 for x in v],float)
def centered(v,loc):
    a=np.array(v,float)
    for block in set(loc):
        mask=np.array(loc)==block; a[mask]-=a[mask].mean()
    return a
def corr(a,b): return float(sum(a*b)/math.sqrt(sum(a*a)*sum(b*b)))
def sign(v): return int(v>1e-12)-int(v< -1e-12)
def adjust(p,method):
    n=len(p); order=sorted(range(n),key=lambda i:p[i]); out=[0.]*n
    for pos,i in enumerate(order):
        if method=='Holm': val=max((n-j)*p[order[j]] for j in range(pos+1))
        else: val=min(n*p[order[j]]/(j+1) for j in range(pos,n))*(sum(1/j for j in range(1,n+1)) if method=='BY' else 1)
        out[i]=min(1.,val)
    return out
def verify(package, output):
    inp=package/'data/analysis_inputs'; tables=output/'tables'; results=output/'validation'
    plan=json.loads((inp/'DESIGN.json').read_text()); receipt=json.loads((results/'ASSOCIATION_EXTENSION_RECEIPT.json').read_text()); checks=[]
    def check(name,ok,detail=None): checks.append(dict(check=name,status='PASS' if bool(ok) else 'FAIL',detail=detail))
    check('Specification precedes replay',plan['created_utc']<receipt['started_utc']<receipt['completed_utc'])
    check('Packaged specification and executed code',sha(inp/'DESIGN.json')==receipt['design_sha256'] and plan['script_sha256']==receipt['script_sha256']==sha(package/'analysis/run_associations.py'))
    for src in plan['sources']:
        check('Staged input '+src['file'],sha(inp/src['file'])==src['sha256'])
    for item in receipt['outputs']: check('Output identity '+item['path'],sha(output/item['path'])==item['sha256'])
    ids=plan['samples']; variables=plan['variables']; fields=plan['metadata_fields']; md={r['Sample ID']:r for r in read(inp/'S16_Study_metadata.tsv')}; loc=[md[s]['Locality'] for s in ids]
    check('Nine matched units, J6 SI3, no P1',len(ids)==len(set(ids))==9 and set(ids)==set(md) and 'P1' not in ids and loc.count('J')==6 and loc.count('SI')==3)
    blocks=[[i for i,x in enumerate(loc) if x==block] for block in ['J','SI']]; per=[]
    for first in itertools.permutations(blocks[0]):
        for second in itertools.permutations(blocks[1]):
            row=list(range(9))
            for dest,src in zip(blocks[0]+blocks[1],first+second): row[dest]=src
            per.append(row)
    check('4320 unique blocked permutations; identity once',len(per)==len(set(map(tuple,per)))==4320 and sum(p==list(range(9)) for p in per)==1 and all(loc[j]==loc[k] for p in per for j,k in enumerate(p)))
    per=np.array(per); env={v:np.array([float(md[s][fields[v]]) for s in ids]) for v in variables}; alpha=read(inp/'Alpha_Primary_Estimates_EN.tsv')
    gt=read(tables/'GENUS_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv'); at=read(tables/'ALPHA_ENVIRONMENT_ALL13_ASSOCIATIONS.tsv'); gs=read(tables/'GENUS_ENVIRONMENT_SAMPLE_DATA.tsv'); ass=read(tables/'ALPHA_ENVIRONMENT_SAMPLE_DATA.tsv')
    loo={(r['analysis_type'],r['marker'],r['response_id'],r['variable'],r['omitted_sample']):r for r in read(tables/'ASSOCIATION_LEAVE_ONE_OUT.tsv')}
    check('Complete unique 312 genus tests',len(gt)==len({(r['marker'],r['feature_id'],r['variable']) for r in gt})==312)
    check('Complete unique 52 alpha tests',len(at)==len({(r['marker'],r['alpha_index'],r['variable']) for r in at})==52)
    check('3276 unique leave-one-out coefficients',len(loo)==3276)
    check('216 genus +36 alpha observation rows',len(gs)==216 and len(ass)==36)
    all_errors=[]; all_exact=[]; summaries=[]
    for marker in ['16S','ITS']:
        full=read(inp/(marker+'_ALL_GENUS_RANK_LINEAGE_COUNTS_DERIVED.tsv')); keys=[r['feature_id'] for r in full]; counts=np.array([[int(r[s]) for s in ids] for r in full]); totals=np.sum(counts,axis=0)
        loghalf=np.log(counts+.5); clr=loghalf-np.sum(loghalf,axis=0)/len(full); logone=np.log(counts+1); clr1=logone-np.sum(logone,axis=0)/len(full)
        source={r['Sample']:r for r in alpha if r['Marker']==marker}; alpha_values={v:np.array([float(source[s][v]) for s in ids]) for v in plan['alpha_indices']}
        errors=[]
        for r in [r for r in gs if r['marker']==marker]:
            i=keys.index(r['feature_id']); j=ids.index(r['sample_id'])
            errors.extend([abs(int(r['count'])-counts[i,j]),abs(int(r['retained_library_count'])-totals[j]),abs(float(r['relative_abundance'])-counts[i,j]/totals[j]),abs(float(r['clr_pc0p5'])-clr[i,j]),abs(float(r['clr_pc1'])-clr1[i,j])])
            errors.extend(abs(float(r[v])-env[v][j]) for v in variables)
        check(marker+' genus full-composition/sample values',max(errors)<1e-10,dict(max_abs_error=float(max(errors))))
        aerrors=[]
        for r in [r for r in ass if r['marker']==marker]:
            j=ids.index(r['sample_id']); aerrors.append(abs(float(r['alpha_value'])-alpha_values[r['alpha_index']][j])); aerrors.extend(abs(float(r[v])-env[v][j]) for v in variables)
            check(marker+' estimator '+r['sample_id']+'/'+r['alpha_index'],r['estimator']==source[r['sample_id']]['Estimator'] and int(r['draws'])==int(source[r['sample_id']]['Draws for estimator']) and int(r['depth'])==int(source[r['sample_id']]['Depth (sequences)']))
        check(marker+' alpha sample values',max(aerrors)<1e-10,dict(max_abs_error=float(max(aerrors))))
        for kind,tests,idfield,qsize in [('genus',gt,'feature_id',156),('alpha',at,'alpha_index',26)]:
            rows=[r for r in tests if r['marker']==marker]; ps=[]; errors=[]; looerrors=[]; exact=[]; signchecks=[]
            check(marker+' '+kind+' complete fixed family',len(rows)==qsize and all(int(r['primary_family_size'])==qsize and int(r['n'])==9 and int(r['total_permutations'])==4320 for r in rows) and all(sum(r['variable']==v for r in rows)==(12 if kind=='genus' else 2) for v in variables))
            for r in rows:
                response=r[idfield]; v=r['variable']; raw=clr[keys.index(response)] if kind=='genus' else alpha_values[response]
                x=centered(rank(env[v]),loc); y=centered(rank(raw),loc); rho=corr(x,y)
                null=np.sum(x[per]*y,axis=1)/math.sqrt(sum(x*x)*sum(y*y)); count=sum(abs(z)>=abs(rho)-1e-12 for z in null); p=count/4320; ps.append(p)
                errors.extend([abs(rho-float(r['rho_locality_adjusted'])),abs(p-float(r['p_two_sided']))]); exact.append(count==int(r['extreme_permutations']))
                loo_values=[]
                for j,sid in enumerate(ids):
                    keep=[i for i in range(9) if i!=j]; subloc=[loc[i] for i in keep]; rr=corr(centered(rank(env[v][keep]),subloc),centered(rank(raw[keep]),subloc)); loo_values.append(rr)
                    looerrors.append(abs(rr-float(loo[kind,marker,response,v,sid]['rho_locality_adjusted'])))
                looerrors.extend([abs(min(loo_values)-float(r['loo_rho_min'])),abs(max(loo_values)-float(r['loo_rho_max']))]); signchecks.append(sum(sign(rr)==sign(rho) for rr in loo_values)==int(r['loo_same_sign_count']))
                if kind=='genus':
                    rho1=corr(x,centered(rank(clr1[keys.index(response)]),loc)); errors.extend([abs(rho1-float(r['rho_pseudocount1'])),abs(abs(rho-rho1)-float(r['delta_rho_pseudocount']))]); signchecks.append((sign(rho)==sign(rho1))==(r['pseudocount_same_sign']=='True'))
            check(marker+' '+kind+' independent rho and exact P',all(exact) and max(errors)<1e-10,dict(max_abs_error=float(max(errors)),extreme_integer_matches=int(sum(exact))))
            check(marker+' '+kind+' independent LOO and signs',max(looerrors)<1e-10 and all(signchecks),dict(max_abs_error=float(max(looerrors))))
            for method in ['BH','BY']:
                q=adjust(ps,method); err=max(abs(a-float(r['q_'+method+str(qsize)])) for a,r in zip(q,rows)); check(marker+' '+kind+' '+method+str(qsize),err<1e-10,dict(max_abs_error=float(err)))
            all_errors.extend(errors+looerrors); all_exact.extend(exact); summaries.append(dict(marker=marker,analysis_type=kind,tests=len(rows),max_abs_error=float(max(errors+looerrors))))
    # Independently refit four preserved alpha category models and 28 environmental
    # models over the 180 distinct category allocations (P denominator 180).
    categories=['Low degradation','Moderate degradation','High degradation']; g=np.array([categories.index(md[s]['Degradation']) for s in ids]); label_alloc=np.unique(g[per],axis=0); check('180 unique category allocations',len(label_alloc)==180)
    x0=np.array([[1.,float(x=='SI')] for x in loc])
    def category_fit(y,labels):
        x1=np.column_stack([x0,labels==1,labels==2]); b0=np.linalg.lstsq(x0,y,rcond=None)[0]; b1=np.linalg.lstsq(x1,y,rcond=None)[0]
        s0=sum((y-x0@b0)**2); s1=sum((y-x1@b1)**2); return max(0.,s0-s1)/2/(s1/5),(s0-s1)/s0
    category=read(inp/'Table_2_Category_Associations_EN.tsv'); alpha_errors=[]; alpha_exact=[]; refitted_alpha=[]
    for r in category:
        if r['Endpoint'] not in ['Observed ASVs','Shannon diversity']: continue
        field='Observed ASVs' if r['Endpoint']=='Observed ASVs' else 'Shannon (ln, nats)'; lookup={a['Sample']:a for a in alpha if a['Marker']==r['Marker']}; y=np.array([float(lookup[s][field]) for s in ids]); f,r2=category_fit(y,g)
        p=sum(category_fit(y,a)[0]>=f-1e-10*max(1,abs(f)) for a in label_alloc)/180
        alpha_errors.extend([abs(f-float(r['F statistic'])),abs(r2-float(r['Partial R2 after locality'])),abs(p-float(r['P raw']))]); alpha_exact.append(abs(p-float(r['P raw']))<1e-10)
        refitted_alpha.append({'Marker':r['Marker'],'Endpoint':r['Endpoint'],'n':9,'df_effect':2,'df_residual':5,'F_statistic_refitted':float(f),'partial_R2_refitted':float(r2),'P_raw_refitted':float(p),'P_Holm_current_family':float(r['P Holm']),'family_size':int(r['Family size']),'scope':'Alpha category refit; family adjustment includes retained beta endpoints'})
    check('Four preserved alpha category fits',len(alpha_exact)==4 and all(alpha_exact) and max(alpha_errors)<1e-10,dict(max_abs_error=float(max(alpha_errors))))
    for marker,n in [('16S',3),('ITS',4)]:
        rows=[r for r in category if r['Marker']==marker]; q=adjust([float(r['P raw']) for r in rows],'Holm'); check(marker+' preserved category Holm'+str(n),len(rows)==n and all(int(r['Family size'])==n for r in rows) and max(abs(v-float(r['P Holm'])) for v,r in zip(q,rows))<1e-10)
    environmental=read(inp/'ENVIRONMENT_CATEGORY_28.tsv'); ee=[]; ep=[]; ek=[]; refitted_environment=[]
    allocation=[]
    for first in itertools.combinations(range(9),3):
        remaining=[i for i in range(9) if i not in first]
        for second in itertools.combinations(remaining,3): allocation.append([first,second,[i for i in remaining if i not in second]])
    allocation=np.array(allocation); check('1680 KW category allocations',len(allocation)==1680)
    for r in environmental:
        vals=np.array([float(md[s][r['variable']]) for s in ids]); y=rank(vals); f,r2=category_fit(y,g); p=sum(category_fit(y,a)[0]>=f-1e-10*max(1,abs(f)) for a in label_alloc)/180; ep.append(p)
        ee.extend([abs(f-float(r['F_partial'])),abs(r2-float(r['partial_R2_rank'])),abs(p-float(r['P_locality']))])
        counts=collections.Counter(vals); tie=1-sum(c**3-c for c in counts.values())/(9**3-9); h=(12/(9*10)*sum(sum(y[g==j])**2/3 for j in range(3))-30)/tie
        hs=(12/(9*10)*np.sum(np.sum(y[allocation],axis=2)**2/3,axis=1)-30)/tie; kp=sum(hs>=h-1e-10*max(1,abs(h)))/1680; ek.append(kp); ee.extend([abs(h-float(r['KW_H'])),abs(kp-float(r['KW_P_exact']))])
        refitted_environment.append({'variable':r['variable'],'n':9,'F_partial_refitted':float(f),'partial_R2_rank_refitted':float(r2),'P_locality_refitted':float(p),'P_BH28_refitted':0.,'P_BY28_refitted':0.,'P_Holm28_refitted':0.,'KW_H_refitted':float(h),'KW_P_exact_refitted':float(kp),'KW_P_BH28_refitted':0.})
    check('All28 environmental partial F/P/R2 and exact KW',len(environmental)==28 and max(ee)<1e-10,dict(max_abs_error=float(max(ee))))
    for method in ['BH','BY','Holm']:
        q=adjust(ep,method); err=max(abs(v-float(r['P_'+method+'28'])) for v,r in zip(q,environmental)); check('Environment preserved '+method+'28',err<1e-10,dict(max_abs_error=float(err)))
    q=adjust(ek,'BH'); check('Environment preserved KW BH28',max(abs(v-float(r['KW_P_BH28'])) for v,r in zip(q,environmental))<1e-10)
    for method in ['BH','BY','Holm']:
        for r,value in zip(refitted_environment,adjust(ep,method)): r['P_'+method+'28_refitted']=float(value)
    for r,value in zip(refitted_environment,adjust(ek,'BH')): r['KW_P_BH28_refitted']=float(value)
    write_rows(results/'REFITTED_ALPHA_CATEGORY.tsv',refitted_alpha)
    write_rows(results/'REFITTED_ENVIRONMENT_CATEGORY.tsv',refitted_environment)
    corrows=read(tables/'ENVIRONMENT_ALL13_CORRELATIONS_DESCRIPTIVE.tsv'); ce=[]
    for r in corrows:
        a,b=r['variable_a'],r['variable_b']; ra,rb=rank(env[a]),rank(env[b]); ce.extend([abs(corr(ra-ra.mean(),rb-rb.mean())-float(r['spearman_rho'])),abs(corr(centered(ra,loc),centered(rb,loc))-float(r['rho_locality_adjusted']))])
    check('78 descriptive dependence pairs',len(corrows)==78 and max(ce)<1e-10)
    out=dict(status='PASS' if all(r['status']=='PASS' for r in checks) else 'FAIL',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),script_sha256=sha(__file__),approach=__doc__,python=sys.version,numpy=np.__version__,checks_count=len(checks),association_exact_integer_matches=int(sum(all_exact)),association_max_abs_error=float(max(all_errors)),summaries=summaries,checks=checks,limits='Independent numerical verification; does not validate assay provenance, laboratory independence, biological identity, exchangeability or causality. The three existing beta raw tests were not refit from distances; their provenance and multiplicity are retained.')
    (results/'ASSOCIATION_INDEPENDENT_QA.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8'); print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
    if out['status']!='PASS': raise AssertionError([r for r in checks if r['status']=='FAIL'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); verify(a.package.resolve(),a.output.resolve())

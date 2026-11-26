"""Read-only arithmetic/provenance audit of retained ANCOM and dispersion tables.

This does not refit ANCOM-BC2, PERMANOVA or PERMDISP models.
"""
from pathlib import Path
import argparse, csv, datetime, hashlib, json, shutil

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def holm(p):
    n=len(p); order=sorted(range(n),key=lambda j:p[j]); result=[0.]*n
    for pos,j in enumerate(order): result[j]=min(1.,max((n-k)*p[order[k]] for k in range(pos+1)))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); package=a.package.resolve(); output=a.output.resolve()
    inp=package/'data/analysis_inputs/retained_statistics'
    names=['S17_ANCOM_global.tsv','S18_ANCOM_contrasts.tsv','S06_Dispersion_diagnostics.tsv']
    sources=json.loads((package/'data/analysis_inputs/DESIGN.json').read_text())['retained_statistics_sources']
    checks=[]
    def check(name,ok,details=None): checks.append(dict(check=name,status='PASS' if bool(ok) else 'FAIL',details=details))
    for src in sources: check('source hash '+src['file'],sha(inp/src['file'])==src['sha256'])
    globalrows=read(inp/names[0]); contrasts=read(inp/names[1]); disp=[r for r in read(inp/names[2]) if r['Test']=='Dispersion diagnostic']; details=[]
    for marker,n in [('16S',356),('ITS',179)]:
        rows=[r for r in globalrows if r['Marker']==marker]; vals=[float(r['Global P raw']) for r in rows]; adjusted=holm(vals)
        err=max(abs(q-float(r['Global P Holm'])) for q,r in zip(adjusted,rows)); check(marker+' complete ANCOM family/Holm',len(rows)==n and all(int(r['Multiplicity family size'])==n for r in rows) and err<1e-10,dict(max_abs_error=err))
        check(marker+' unique global lineage keys',len({r['Full-lineage feature ID'] for r in rows})==n)
        cr=[r for r in contrasts if r['Marker']==marker]; lookup={r['Full-lineage feature ID']:r for r in rows}; errs=[]
        check(marker+' two contrast rows per lineage',len(cr)==2*n and all(sum(r['Full-lineage feature ID']==k for r in cr)==2 for k in lookup))
        for r in cr:
            g=lookup[r['Full-lineage feature ID']]; coef=float(r['Coefficient (natural log)']); se=float(r['Standard error'])
            errs.extend([abs((coef-se)-float(r['Coefficient - 1 SE'])),abs((coef+se)-float(r['Coefficient + 1 SE'])),abs(float(r['Global P raw'])-float(g['Global P raw'])),abs(float(r['Global P Holm'])-float(g['Global P Holm']))])
        check(marker+' contrast coefficients +/-1SE and global P consistency',max(errs)<1e-10,dict(max_abs_error=max(errs)))
        details.append(dict(marker=marker,global_rows=n,contrast_rows=len(cr),global_P_Holm_below_0p05=sum(float(r['Global P Holm'])<.05 for r in rows),global_robust_discoveries=sum(r['Global robust discovery flag']=='True' for r in rows),valid_primary_degradation_endpoints=sum(r['Primary degradation endpoint valid']=='True' for r in rows),auxiliary_invalid_rows=sum(r['All primary coefficients valid']!='True' for r in rows)))
    its=[r for r in disp if r['Marker']=='ITS']; iq=holm([float(r['P raw']) for r in its]); check('Retained ITS dispersion Holm2',len(its)==2 and max(abs(q-float(r['P adjusted'])) for q,r in zip(iq,its))<1e-10)
    check('Retained dispersion 3 diagnostics',len(disp)==3 and sum(r['Marker']=='16S' for r in disp)==1)
    out=dict(status='PASS' if all(r['status']=='PASS' for r in checks) else 'FAIL',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,source_bindings=sources,summaries=details,checks=checks,category_scope='Current dispersion-only input; no obsolete category records are included.',script_sha256=sha(__file__))
    result=output/'validation/RETAINED_STATISTICS_QA.json'; result.parent.mkdir(parents=True,exist_ok=True); result.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8'); print(json.dumps(out,indent=2))
    if out['status']!='PASS': raise AssertionError([r for r in checks if r['status']=='FAIL'])

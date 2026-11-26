"""Verify current sample metadata against the preserved preparation record."""
from pathlib import Path
import argparse, csv, json

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args(); inp=a.package/'data/analysis_inputs'
    current={r['Sample ID']:r for r in read(inp/'S16_Study_metadata.tsv')}
    source={r['SampleID']:r for r in read(inp/'metadata_preparation_record.tsv')}
    fields={'Degradation':'zone','Locality':'region','Vegetation':'vegetation','Soil texture':'soil_type','Sand (%)':'sand','Silt (%)':'silt','Clay (%)':'clay','pH (dimensionless)':'pH','Conductivity (unit unconfirmed)':'conductivity','Calcium carbonate (unit unconfirmed)':'Calcium carbonate','Total nitrogen (unit unconfirmed)':'Total Nitrogen','Available phosphorus (unit unconfirmed)':'available phosphorus','Available potassium (unit unconfirmed)':'available potassium','Cation exchange capacity (unit unconfirmed)':'CIC','Exchangeable acidity (unit unconfirmed)':'Interchangeable acidity','Exchangeable aluminum (unit unconfirmed)':'Interchangeable aluminum','Organic matter (unit unconfirmed)':'Organic matter','Total carbon (unit unconfirmed)':'Total carbon'}
    for f in ['Calcium','Magnesium','Sodium','Potassium','Cu','Fe','Pb','Zn','Ba','Al','Mn','Co','Sr','Cr']:
        fields[f+' (unit unconfirmed)']=f
    checks=[]
    def check(name,condition): checks.append({'check':name,'status':'PASS' if bool(condition) else 'FAIL'})
    check('Nine identical sample identifiers',set(current)==set(source)==set(['A1','A2','A3','B1','B2','B3','M1','M2','M3']))
    for sid,r in current.items():
        for k,f in fields.items(): check(sid+'/'+k,r[k]==source[sid][f] or (k not in ['Degradation','Locality','Vegetation','Soil texture'] and float(r[k])==float(source[sid][f])))
        original=[float(source[sid][f+'_original']) for f in ['sand','silt','clay']]
        used=[float(r[f+' (%)']) for f in ['Sand','Silt','Clay']]
        check(sid+' analytical texture sum',abs(sum(used)-100)<1e-9)
        if sid=='A1':
            check('A1 recorded original sum 102',abs(sum(original)-102)<1e-9)
            check('A1 normalized rounded values',all(round(x*100/102,2)==y for x,y in zip(original,used)))
        else: check(sid+' texture unchanged',original==used)
    result={'status':'PASS' if all(x['status']=='PASS' for x in checks) else 'FAIL','checks':checks,'scope':'Record consistency and computational preparation only; not a laboratory correction or confirmation of assay units/methods.'}
    a.report.parent.mkdir(parents=True,exist_ok=True)
    a.report.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':len(checks)}))
    if result['status']!='PASS': raise AssertionError([x for x in checks if x['status']!='PASS'])

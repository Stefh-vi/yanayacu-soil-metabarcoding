"""Check all released bytes against SHA256_MANIFEST.tsv using only Python's standard library."""
from pathlib import Path
import argparse, csv, hashlib, json, sys

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--report',type=Path)
    a=ap.parse_args(); root=a.root.resolve()
    if a.report and (a.report.resolve()==root or root in a.report.resolve().parents):
        ap.error('Write the verification report outside the immutable release.')
    manifest=root/'SHA256_MANIFEST.tsv'
    with manifest.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    problems=[]; seen=set()
    for row in rows:
        rel=row['path']; p=(root/rel).resolve()
        if rel in seen or root not in p.parents:
            problems.append({'path':rel,'error':'Duplicate or unsafe path'}); continue
        seen.add(rel)
        if not p.is_file(): problems.append({'path':rel,'error':'Missing'}); continue
        if p.stat().st_size!=int(row['bytes']) or digest(p)!=row['sha256']:
            problems.append({'path':rel,'error':'Bytes differ'})
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and '.git' not in p.relative_to(root).parts and '__pycache__' not in p.relative_to(root).parts and p!=manifest}
    for rel in sorted(actual-seen): problems.append({'path':rel,'error':'Unlisted file'})
    report={'status':'PASS' if not problems else 'FAIL','manifest_entries':len(rows),'problems':problems,'scope':'File identity only. This does not establish biological validity or repeat analyses.'}
    if a.report: a.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    return 0 if not problems else 1

if __name__=='__main__': sys.exit(main())

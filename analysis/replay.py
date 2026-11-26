"""Reproduce the current numerical analyses and compare released table hashes."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--output', type=Path, required=True)
    a=p.parse_args(); package=a.package.resolve(); output=a.output.resolve()
    if output.exists(): raise FileExistsError('Choose a new output directory; prior outputs are preserved.')
    if output==package or package in output.parents: raise ValueError('Replay must write outside the released package.')
    for script in ['run_associations.py','verify_independently.py','verify_retained_statistics.py']:
        subprocess.run([sys.executable,str(package/'analysis'/script),'--package',str(package),'--output',str(output)],check=True)
    expected=json.loads((package/'analysis/EXPECTED_TABLES.json').read_text())
    rows=[]
    for item in expected:
        digest=sha(output/item['path'])
        rows.append({**item,'replayed_sha256':digest,'equal':digest==item['sha256']})
    result={'status':'PASS' if all(r['equal'] for r in rows) else 'FAIL','tables_compared':len(rows),'results':rows,
            'scope':'Current association output identity and independent arithmetic checks. Sequence processing, taxonomy, beta-distance models, rarefaction, ANCOM-BC2 and PERMDISP are not refitted by this command.'}
    (output/'validation/REPLAY_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if result['status']!='PASS': raise AssertionError('Current replay differs from released output; investigate without replacing reference tables.')

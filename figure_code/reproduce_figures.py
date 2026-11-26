"""Reproduce current figure surfaces in a new output directory.

Thirteen figures are rebuilt from frozen numeric tables. Two native vector figures
receive the exact source-derived text revisions used by the current edition.
Sixteen other current figures are retained byte-for-byte. No model is fitted.
"""
from pathlib import Path
import argparse,csv,hashlib,html,json,os,re,shutil,subprocess,sys,platform
import numpy as np
import pandas as pd
import matplotlib
import PIL
import pymupdf

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/'figure_code'
NUMERIC_RENDER={'02','03','05','06','S03','S04','S05','S06','S07','S16','S17','S21','S23'}
VECTOR_TEXT={'S09','S11'}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def gallery(out,reg):
    cards=[]
    for row in reg:
        base=Path(row['scope'])/row['figure']
        cap=(ROOT/row['caption_path']).read_text(encoding='utf-8')
        links=' | '.join('<a href="'+base.with_suffix('.'+ext).as_posix()+'">'+ext.upper()+'</a>' for ext in ['pdf','svg','png'])
        cards.append('<article><h2>'+html.escape(row['heading'])+'</h2><img loading="lazy" src="'+base.with_suffix('.png').as_posix()+'" alt="'+html.escape(row['heading'])+'"><p>'+links+'</p><details><summary>Legend</summary><pre>'+html.escape(cap)+'</pre></details></article>')
    page='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Yanayacu current figures</title><style>body{font:16px Arial;max-width:1500px;margin:auto;padding:24px;color:#1c303a}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(510px,1fr));gap:20px}article{border:1px solid #ccc;padding:18px}img{width:100%;height:440px;object-fit:contain}pre{white-space:pre-wrap;font:15px Arial;line-height:1.5}h2{font-size:19px}</style><h1>Yanayacu current scientific figures</h1><p>Six main figures and 25 supplementary figures. Statistical families and interpretation are stated in each legend.</p><main>'+''.join(cards)+'</main></html>'
    (out/'FIGURE_GALLERY.html').write_text(page,encoding='utf-8')

def verify(out,reg,qa):
    records=[]
    for row in reg:
        base=Path(row['scope'])/row['figure'];target=out/base.with_suffix('.pdf');source=ROOT/'figures'/base.with_suffix('.pdf')
        with pymupdf.open(source) as expected,pymupdf.open(target) as actual:
            assert len(expected)==len(actual)==1,row['figure_id']
            # PDF text strips contain each rendered numeric annotation, axis tick and sample ID.
            et=expected[0].get_text();at=actual[0].get_text()
            tokens=lambda t:re.findall(r'[-+]?\d+(?:[,.]\d+)*',t)
            en,an=tokens(et),tokens(at)
            numeric_equal=en==an
            ep=expected[0].get_pixmap(dpi=100,alpha=False);ap=actual[0].get_pixmap(dpi=100,alpha=False)
            shape_equal=(ep.width,ep.height)==(ap.width,ap.height)
            pixel_equal=shape_equal and ep.samples==ap.samples
            mae=None
            if shape_equal:mae=float(np.abs(np.frombuffer(ep.samples,dtype=np.uint8).astype(float)-np.frombuffer(ap.samples,dtype=np.uint8).astype(float)).mean())
            rendered=qa/'renders';rendered.mkdir(exist_ok=True);ap.save(rendered/(row['figure_id']+'.png'))
        untouched=row['figure_id'] not in NUMERIC_RENDER|VECTOR_TEXT
        files_equal=all(digest(ROOT/f['path'])==digest(out/Path(f['path']).relative_to('figures')) for f in row['files'])
        records.append({'figure_id':row['figure_id'],'scope':row['reproduction_scope'],'numeric_tokens_match':numeric_equal,'expected_numeric_token_count':len(en),'rendered_numeric_token_count':len(an),'page_dimensions_match':shape_equal,'rendered_pixels_equal_at_100dpi':pixel_equal,'mean_absolute_pixel_difference_0_255':mae,'all_native_files_byte_identical':files_equal,'native_preservation_required':untouched,'status':'PASS' if numeric_equal and shape_equal and (not untouched or files_equal) else 'FAIL'})
        if not numeric_equal:
            (qa/(row['figure_id']+'_NUMERIC_TOKEN_DIFF.json')).write_text(json.dumps({'expected':en,'actual':an},indent=2),encoding='utf-8')
    result={'status':'PASS' if all(x['status']=='PASS' for x in records) else 'FAIL','figure_sets':len(reg),'rendered_from_numeric_tables':len(NUMERIC_RENDER),'native_vector_text_updates':len(VECTOR_TEXT),'native_sets_retained_byte_identically':len(reg)-len(NUMERIC_RENDER)-len(VECTOR_TEXT),'scope':'Rendering/numeric annotation and native-artifact preservation check. No statistical or biological validation; font/render versions may change PDF bytes and pixels.','records':records}
    (qa/'FIGURE_REPRODUCTION_CHECK.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if result['status']!='PASS':raise RuntimeError('Figure readback failed; inspect FIGURE_REPRODUCTION_CHECK.json')
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    out=args.output.resolve()
    if out.exists():raise SystemExit('Choose a new, nonexistent output directory. Current publication files are never overwritten.')
    checks=[]
    for r in csv.DictReader((ROOT/'validation/figures/SOURCE_INVENTORY.tsv').open(encoding='utf-8'),delimiter='\t'):
        p=ROOT/r['release_path'];ok=digest(p)==r['release_sha256'];checks.append({'path':r['release_path'],'hash_matches':ok})
    if not all(r['hash_matches'] for r in checks):raise RuntimeError('Input/source identity mismatch')
    out.mkdir(parents=True);qa=out/'validation';qa.mkdir()
    (qa/'INPUT_CHECK.json').write_text(json.dumps({'status':'PASS','records':checks},indent=2),encoding='utf-8')
    reg=json.loads((ROOT/'validation/figures/FIGURE_REGISTER.json').read_text(encoding='utf-8'))
    for scope in ['MAIN','SUPPLEMENTARY']:shutil.copytree(ROOT/'figures'/scope,out/scope)
    env=os.environ.copy();env['YANAYACU_FIGURE_OUTPUT']=str(out);env['YANAYACU_FIGURE_WORK']=str(qa)
    def run(script,*args):
        result=subprocess.run([sys.executable,str(CODE/script),*map(str,args)],env=env,text=True,capture_output=True,encoding='utf-8',errors='replace')
        (qa/(script+'.log')).write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
        if result.returncode:raise RuntimeError('Rendering failed: '+script+'; inspect its log')
    run('build_associations.py');run('build_alpha.py','--package',ROOT);run('build_beta.py')
    for ext in ['pdf','svg','png']:
        suffix='_600dpi.png' if ext=='png' else '.'+ext
        shutil.copy2(qa/'alpha_render'/('Figure_02_Alpha_Diversity_Revised'+suffix),out/'MAIN'/('Figure_02_Alpha_Diversity_Revised.'+ext))
    run('annotate_s09.py');run('clarify_s11_title.py')
    caps=json.loads((ROOT/'data/figure_inputs/FIGURE_CAPTIONS.json').read_text(encoding='utf-8'))
    for row in reg:
        base=out/row['scope']/row['figure'];cap=caps[row['figure_id']]
        base.with_suffix('.caption.md').write_text('# '+cap['heading']+'\n\n'+cap['body']+'\n',encoding='utf-8')
        if row['figure_id'] in {'05','06','S25'}:
            p=base.with_suffix('.pdf');d=pymupdf.open(p);meta=d.metadata
            if meta.get('title')!=row['heading']:
                meta['title']=row['heading'];d.set_metadata(meta);tmp=p.with_suffix('.metadata.pdf');d.save(tmp);d.close();tmp.replace(p)
            else:d.close()
    environment={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'Pillow':PIL.__version__,'PyMuPDF':pymupdf.__version__}
    (qa/'RENDER_ENVIRONMENT.json').write_text(json.dumps(environment,indent=2)+'\n',encoding='utf-8')
    gallery(out,reg);result=verify(out,reg,qa)
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))

if __name__=='__main__':main()

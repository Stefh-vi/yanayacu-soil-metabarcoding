"""Add source-derived test probabilities to the native S09 footer only."""
from pathlib import Path
import csv,json,hashlib
import pymupdf
from matplotlib.font_manager import findfont

from paths import ROOT, INP, OUT, QA
src=INP/'native/S09'
out=OUT/'SUPPLEMENTARY/Figure_S09_ITS_Alternative_Ordination_REVISED'
rows=list(csv.DictReader((INP/'tables/T02_Category_associations.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
beta=[r for r in rows if r['Marker']=='ITS' and ('Bray' in r['Endpoint'] or 'Hellinger' in r['Endpoint'])]
assert len(beta)==2
assert len({r['P raw'] for r in beta})==len({r['P Holm'] for r in beta})==1
assert {int(r['Family size']) for r in beta}=={4}
old='Axis percentages use positive inertia; existing conditional tests are in the caption.'
new=f"Both distances: P = {float(beta[0]['P raw']):.4f}; Holm-adjusted P = {float(beta[0]['P Holm']):.4f} (4 endpoints)."
svg=(src/'Figure_S09_ITS_Alternative_Ordination.svg').read_text(encoding='utf-8')
assert svg.count(old)==1
out.with_suffix('.svg').write_text(svg.replace(old,new),encoding='utf-8')
doc=pymupdf.open(src/'Figure_S09_ITS_Alternative_Ordination.pdf');page=doc[0]
hits=page.search_for(old);assert len(hits)==1
spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
line_spans=[s for s in spans if pymupdf.Rect(s['bbox']).intersects(hits[0])]
span=min(line_spans,key=lambda s:s['origin'][0])
fontpath=findfont('DejaVu Sans');font=pymupdf.Font(fontfile=fontpath)
assert font.text_length(new,fontsize=span['size']) < page.rect.width-span['origin'][0]-20
page.add_redact_annot(hits[0],fill=(1,1,1))
page.apply_redactions(images=0,graphics=0,text=0)
page.insert_font(fontname='V3DejaVu',fontbuffer=Path(fontpath).read_bytes())
page.insert_text(span['origin'],new,fontname='V3DejaVu',fontsize=span['size'],color=(83/255,107/255,117/255))
doc.save(out.with_suffix('.pdf'),garbage=4,deflate=True)
doc.close()
with pymupdf.open(out.with_suffix('.pdf')) as d:d[0].get_pixmap(dpi=400).save(out.with_suffix('.png'))
receipt={'scope':'Native footer text only; ordination points and axes unchanged','source_pdf_sha256':hashlib.sha256((src/'Figure_S09_ITS_Alternative_Ordination.pdf').read_bytes()).hexdigest(),'old_text':old,'new_text':new,'tests':beta}
(QA/'S09_ANNOTATION_CHANGE.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(new)

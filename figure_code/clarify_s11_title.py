"""Clarify the technical sensitivity title; retain all quantitative artwork."""
from pathlib import Path
import shutil,json,hashlib
import pymupdf
from matplotlib.font_manager import findfont
from paths import ROOT, INP, OUT, QA
stem='Figure_S11_ITS_Error_Model_Sensitivity'
out=OUT/'SUPPLEMENTARY'/stem
inputs=INP/'native/S11'
for ext in ['pdf','svg']:
    p=inputs/(stem+'.'+ext)
    assert p.is_file(), p.name
old='Fungal community response to balanced error learning'
new='Fungal community sensitivity to balanced error learning'
svg=(inputs/(stem+'.svg')).read_text(encoding='utf-8');assert svg.count(old)==1
out.with_suffix('.svg').write_text(svg.replace(old,new),encoding='utf-8')
doc=pymupdf.open(inputs/(stem+'.pdf'));p=doc[0];hits=p.search_for(old);assert len(hits)==1
spans=[s for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if pymupdf.Rect(s['bbox']).intersects(hits[0])]
s=min(spans,key=lambda v:v['origin'][0]);fp=Path(findfont('DejaVu Sans',fallback_to_default=True));
# The native heading uses a bold font.
from matplotlib.font_manager import FontProperties
fp=Path(findfont(FontProperties(family='DejaVu Sans',weight='bold')))
font=pymupdf.Font(fontbuffer=fp.read_bytes());assert font.text_length(new,fontsize=s['size'])<p.rect.width-s['origin'][0]-10
p.add_redact_annot(hits[0],fill=(1,1,1));p.apply_redactions(images=0,graphics=0,text=0)
p.insert_font(fontname='V3Heading',fontbuffer=fp.read_bytes())
color=tuple(((s['color']>>shift)&255)/255 for shift in [16,8,0])
p.insert_text(s['origin'],new,fontname='V3Heading',fontsize=s['size'],color=color)
doc.save(out.with_suffix('.pdf'),garbage=4,deflate=True);doc.close()
with pymupdf.open(out.with_suffix('.pdf')) as d:d[0].get_pixmap(dpi=400).save(out.with_suffix('.png'))
(QA/'S11_TITLE_CHANGE.json').write_text(json.dumps({'old':old,'new':new,'scope':'Title only; technical quantities and all plotted geometry unchanged','source_sha256':hashlib.sha256((inputs/(stem+'.pdf')).read_bytes()).hexdigest()},indent=2),encoding='utf-8')
print(new)

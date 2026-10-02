from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
assert not (root/'manuscript').exists()
assert not (root/'publications/Four_Dimensional_Connectedness_Loci.pdf').exists()
assert hashlib.sha256((root/'publications/BMD2026_Poster.pdf').read_bytes()).hexdigest()=='7fcba8b565c848c18a668c4b46f368079477179bef952ce578dafd01fdc8adc2'
for name in ['citation.bib','bibitems.tex']:
 text=(root/name).read_text()
 assert text.count('provisional title')==2
 assert 'in preparation' in text
print('Submitted poster verified; manuscript files withheld; provisional citations present.')

#!/usr/bin/env python3
"""Extract a hero page, dropping unused resources without altering its artwork.

Matplotlib shares an XObject dictionary across pages. Keeping that full dictionary
when extracting one page copies images used only by other pages. Pruning by the
page's actual Do operators is lossless and retains all original image resolution.
"""
import argparse, hashlib, json
from pathlib import Path
from pypdf import PdfReader, PdfWriter
from pypdf.generic import ContentStream, DictionaryObject, NameObject
import pypdf.generic._base as pdf_numbers
import fitz
import numpy as np

# Preserve the source's decimal page boxes; the library's default rounding can
# change only the antialiased white page-edge row by one RGB level.
pdf_numbers.FLOAT_WRITE_PRECISION = 14

parser=argparse.ArgumentParser()
parser.add_argument('input',type=Path)
parser.add_argument('output',type=Path)
parser.add_argument('--page',type=int,default=1)
args=parser.parse_args()
reader=PdfReader(args.input);page=reader.pages[args.page-1]
resources=DictionaryObject(page['/Resources'])
objects=resources.get('/XObject',DictionaryObject()).get_object()
used={str(operands[0]) for operands,operator in ContentStream(page.get_contents(),reader).operations if operator==b'Do'}
def write_with_objects(retained):
    active=DictionaryObject(resources)
    active[NameObject('/XObject')]=DictionaryObject({key:value for key,value in objects.items() if str(key) in retained})
    page[NameObject('/Resources')]=active
    writer=PdfWriter();writer.pdf_header=reader.pdf_header;writer.add_page(page)
    writer.add_metadata({'/Title':'Attractors and contact differences in a four-dimensional atlas','/Author':'Bernat Espigulé'})
    with args.output.open('wb') as stream:writer.write(stream)

write_with_objects(used)
original=fitz.open(args.input)[args.page-1].get_pixmap(matrix=fitz.Matrix(.8,.8),alpha=False)
extracted=fitz.open(args.output)[0].get_pixmap(matrix=fitz.Matrix(.8,.8),alpha=False)
assert original.width==extracted.width and original.height==extracted.height
retained_for_transparency=[]
if original.samples!=extracted.samples:
    # MuPDF scans unused alpha-image resources to choose its page compositing
    # path. Removing every alpha image from an otherwise opaque page can change
    # only a white page-edge row by one RGB level. Retain the smallest such
    # resource to preserve even this renderer behaviour exactly.
    a=np.frombuffer(original.samples,np.uint8).reshape(original.height,original.width,3)
    b=np.frombuffer(extracted.samples,np.uint8).reshape(extracted.height,extracted.width,3)
    delta=np.abs(a.astype(int)-b.astype(int))
    assert delta.max()<=1 and not np.any(delta[1:-1,1:-1]),'Unexpected artwork difference.'
    candidates=[]
    for key,value in objects.items():
        obj=value.get_object()
        if str(key) not in used and '/SMask' in obj:
            size=len(obj._data)+len(obj['/SMask'].get_object()._data)
            candidates.append((size,str(key)))
    for _,key in sorted(candidates):
        write_with_objects(used|{key})
        extracted=fitz.open(args.output)[0].get_pixmap(matrix=fitz.Matrix(.8,.8),alpha=False)
        if original.samples==extracted.samples:
            retained_for_transparency=[key]
            break
assert original.samples==extracted.samples,'The extracted figure must be pixel-identical.'
record={'input':args.input.name,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
        'page':args.page,'xobjects_before':len(objects),'xobjects_used':len(used),
        'retained_for_transparency':retained_for_transparency,
        'pixel_identical':True,'render_dimensions':[original.width,original.height],
        'output':args.output.name,'output_bytes':args.output.stat().st_size,
        'output_sha256':hashlib.sha256(args.output.read_bytes()).hexdigest()}
args.output.with_suffix('.extraction.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))

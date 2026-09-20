#!/usr/bin/env python3
"""Check the self-contained page without network access.
Requires Playwright and Chromium. CHROMIUM_PATH optionally selects a binary.
"""
from pathlib import Path
import json,os,shutil
from playwright.sync_api import sync_playwright
root=Path(__file__).resolve().parents[1]
errors=[];requests=[]
with sync_playwright() as p:
 args={'headless':True,'args':['--no-sandbox','--enable-unsafe-swiftshader']}
 exe=os.environ.get('CHROMIUM_PATH') or shutil.which('chromium') or shutil.which('google-chrome')
 if exe:args['executable_path']=exe
 browser=p.chromium.launch(**args)
 page=browser.new_page(viewport={'width':1440,'height':1100})
 page.on('pageerror',lambda e:errors.append(str(e)))
 page.on('request',lambda r:requests.append(r.url) if not r.url.startswith(('data:','file:')) else None)
 page.set_content((root/'interactive/Four_Family_Stable_Atlas.html').read_text(),wait_until='load')
 page.wait_for_function('window.ATLAS_QA && ATLAS_QA.renderCount>0')
 data=page.evaluate('({families:ATLAS_QA.familyNames,curves:Object.keys(ATLAS_QA.curves).length,overflow:document.documentElement.scrollWidth>innerWidth})')
 assert data['families']==['DD10','DD11','DO10','DO11'] and data['curves']==12 and not data['overflow']
 page.locator('#allWeights').uncheck()
 page.locator('#weight').evaluate("el=>{el.value='0.35';el.dispatchEvent(new Event('input',{bubbles:true}));}")
 assert page.evaluate('ATLAS_QA.state.weight')==.35
 for family in data['families']:
  for probe in 'ABC':page.locator(f'.curve-button[data-family="{family}"][data-probe="{probe}"]').click()
 page.set_viewport_size({'width':390,'height':844});page.wait_for_timeout(120)
 assert page.evaluate('document.documentElement.scrollWidth<=innerWidth')
 assert not errors and not requests
 browser.close()
print(json.dumps({'status':'passed','mode':'inline self-contained HTML','desktop_width':1440,'mobile_width':390,'curve_examples_checked':12,'runtime_errors':errors,'external_requests':requests},indent=2))

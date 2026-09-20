#!/usr/bin/env python3
"""Bounded offline verification of the published exact records."""
from pathlib import Path
import importlib.util, json, subprocess, sys
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("independent",root/"computation/verify_certificate.py")
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
results=[]
for p in sorted((root/"computation/certificates").glob("*.json")):
 d=json.loads(p.read_text())
 if not isinstance(d,dict) or d.get("schema")!="rational-planar-ifs-exclusion-1":continue
 results.append({"file":p.name,**mod.verify(d)})
assert len(results)==7
for name in ("verify_contact_identities.py","verify_fixed_address_independently.py","check_exact_families.py"):
 subprocess.run([sys.executable,str(root/"verification"/name)],cwd=root,check=True,timeout=150)
print(json.dumps({"status":"passed","certificates":results},indent=2))

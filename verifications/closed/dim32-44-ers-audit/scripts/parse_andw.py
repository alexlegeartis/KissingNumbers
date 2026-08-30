#!/usr/bin/env python3
"""Parse Brouwer's Andw.html constant-weight-code tables into (n,d,w) -> (lo,hi,ref)."""
import re, sys, json, os
HERE=os.path.dirname(os.path.abspath(__file__))
path=os.path.join(HERE,'..','data','Andw.html')
html=open(path,encoding='utf-8').read()

# split into sections by <h1><a name="dK">
secs=re.split(r'<h1><a name="d(\d+)">',html)
out={}
for i in range(1,len(secs),2):
    d=int(secs[i]); body=secs[i+1]
    # take the first <table ...> ... </table>
    m=re.search(r'<table[^>]*>(.*?)</table>',body,re.S)
    if not m: continue
    tbl=m.group(1)
    rows=re.findall(r'<tr>(.*?)</tr>',tbl,re.S)
    header=None
    for r in rows:
        cells=re.findall(r'<t([hd])[^>]*>(.*?)</t[hd]>',r,re.S)
        if not cells: continue
        if cells[0][0]=='h' and header is None and len(cells)>2 and 'n\w' in cells[0][1]:
            header=[c[1] for c in cells[1:]]
            header=[int(re.sub(r'<[^>]+>','',h).strip()) for h in header]
            continue
        if cells[0][0]!='h': continue
        try: n=int(re.sub(r'<[^>]+>','',cells[0][1]).strip())
        except ValueError: continue
        for j,(t,c) in enumerate(cells[1:]):
            if header is None or j>=len(header): continue
            w=header[j]
            raw=c
            txt=re.sub(r'<sup>.*?</sup>','',raw,flags=re.S)
            link=re.search(r'href="([^"]+)"',raw)
            txt=re.sub(r'<[^>]+>','',txt)
            txt=txt.replace('&nbsp;',' ').strip()
            if not txt or txt=='-': continue
            # forms: "1667", "46-49", "48-58"
            mm=re.match(r'^\s*([0-9]+)\s*(?:[-\u2013]\s*([0-9]+))?\s*$',txt)
            if not mm: 
                out.setdefault('UNPARSED',[]).append((d,n,w,txt))
                continue
            lo=int(mm.group(1)); hi=int(mm.group(2)) if mm.group(2) else lo
            sups=re.findall(r'<sup>(.*?)</sup>',raw,re.S)
            out[(n,d,w)]=(lo,hi,link.group(1) if link else None,''.join(re.sub(r'<[^>]+>','',s) for s in sups))
if 'UNPARSED' in out:
    print("UNPARSED:",out.pop('UNPARSED'),file=sys.stderr)
ks=sorted(k for k in out if isinstance(k,tuple))
print("total entries:",len(ks))
print("d values:",sorted(set(k[1] for k in ks)))
print("n range:",min(k[0] for k in ks),max(k[0] for k in ks))
json.dump({f"{k[0]},{k[1]},{k[2]}":out[k] for k in ks},open(os.path.join(HERE,'..','data','andw_parsed.json'),'w'),indent=0)

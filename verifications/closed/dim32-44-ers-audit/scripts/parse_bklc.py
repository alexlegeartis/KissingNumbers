#!/usr/bin/env python3
"""Parse a codetables.de BKLC grid (n rows, k cols) -> best known linear [n,k,d]."""
import re,os,sys,json
HERE=os.path.dirname(os.path.abspath(__file__))
h=open(os.path.join(HERE,'..','data',sys.argv[1] if len(sys.argv)>1 else 'bklc_30_50.html'),encoding='latin-1').read()
m=re.search(r'<TABLE BORDER=1 CELLPADDING=2.*?</TABLE>',h,re.S|re.I)
tbl=m.group(0)
rows=re.findall(r'<TR>(.*?)</TR>',tbl,re.S|re.I)
hdr=[re.sub(r'<[^>]+>','',c).strip() for c in re.findall(r'<TD[^>]*>(.*?)</TD>',rows[0],re.S|re.I)]
ks=[int(x) for x in hdr[1:]]
D={}
for r in rows[1:]:
    cells=[re.sub(r'<[^>]+>','',c).strip().replace('&nbsp;',' ') for c in re.findall(r'<TD[^>]*>(.*?)</TD>',r,re.S|re.I)]
    if not cells: continue
    try: n=int(cells[0])
    except ValueError: continue
    for k,c in zip(ks,cells[1:]):
        c=c.strip()
        mm=re.match(r'^(\d+)(?:\s*[-\u2013]\s*(\d+))?$',c)
        if not mm: continue
        D[(n,k)]=int(mm.group(1))   # LOWER end = achievable
json.dump({f"{a},{b}":D[(a,b)] for a,b in sorted(D)},open(os.path.join(HERE,'..','data','bklc_parsed.json'),'w'),indent=0)
# best k for each (n,d): largest k with d(n,k) >= d
print("Best known LINEAR lower bounds A(n,d) >= 2^k  (codetables.de, Grassl)")
print("%4s %4s %5s %14s"%("n","d","k","2^k"))
for d in range(6,14):
    print("-"*36)
    for n in range(30,51):
        best=max([k for (nn,k) in D if nn==n and D[(nn,k)]>=d]+[0])
        if best: print("%4d %4d %5d %14d"%(n,d,best,1<<best))

#!/usr/bin/env python3
"""Parse the Litsyn-Rains-Sloane 'Table of Nonlinear Binary Codes' snapshot into
lower bounds A(n,d) >= N*2^k for odd d, plus the derived even-d values."""
import re, html, os, sys, json
HERE=os.path.dirname(os.path.abspath(__file__))
s=open(os.path.join(HERE,'..','data','LRS_tableand.html'),encoding='latin-1').read()
KEY={}
kidx=s.lower().find('<a name="key"')
for line in html.unescape(re.sub(r'<[^>]+>','',s[kidx:kidx+4000])).splitlines():
    m=re.match(r'\s*([A-Z0-9]{1,2})\s*=\s*(.+?)\s*$',line)
    if m: KEY[m.group(1)]=m.group(2)
def table(d):
    i=s.find('NAME="dist%d"'%d); j=s.find('NAME="dist%d"'%(d+2))
    seg=s[i:j if j>0 else len(s)]; rows=[]
    for m in re.finditer(r'<tr>(.*?)(?=<tr>|</table>|$)',seg,re.S|re.I):
        c=re.findall(r'<t[dh]>(.*?)</t[dh]>',m.group(1),re.S|re.I)
        if len(c)<6: continue
        c=[html.unescape(re.sub(r'<[^>]+>','',x)).strip() for x in c]
        if not c[0].isdigit(): continue
        rows.append((int(c[0]),int(c[1]),int(c[2]),int(c[3]),c[4],c[5]))
    return rows
DB={}
for d in (3,5,7,9,11,13):
    rows=table(d)
    listed={r[0]:r for r in rows}
    ns=sorted(listed)
    for n in range(4,60):
        if n in listed:
            n0,dd,k,N,ty,ref=listed[n]; DB[(n,d)]=(N<<k,'listed',ty,ref,n0)
        else:
            nxt=[x for x in ns if x>n]
            if not nxt: continue
            n0=nxt[0]; _,dd,k,N,ty,ref=listed[n0]; kk=k-(n0-n)
            if kk<0: continue
            DB[(n,d)]=(N<<kk,'shortened from n=%d'%n0,ty,ref,n0)
print("LRS lower bounds A(n,d) (snapshot 'last updated November 24, 1999')")
print("%4s %10s %14s  %-22s %-4s %s"%("n","d","A(n,d)>=","provenance","type","ref"))
for d in (7,9,11,13):
    print("-"*92)
    for n in range(26,50):
        if (n,d) in DB:
            v,prov,ty,ref,n0=DB[(n,d)]
            print("%4d %10d %14d  %-22s %-4s %-10s %s"%(n,d,v,prov,ty,ref,KEY.get(ty,'')))
json.dump({f"{k[0]},{k[1]}":DB[k] for k in DB},open(os.path.join(HERE,'..','data','lrs_parsed.json'),'w'),indent=0)

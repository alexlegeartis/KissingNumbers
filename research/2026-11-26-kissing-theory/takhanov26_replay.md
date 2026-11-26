# Reproducing the local 26D obstruction from the upstream certificate

Run with Takhanov's pinned release 77f0860bf93a745140e6e25481d3f60773b11441, in version2/certificates. The last step can be run outside that directory if the two JSON file locations below are updated. This is a finite replay of the two reported obstructions; the sharp 1/9 proof is analytic in the accompanying long-form note.

```python
from common import leech_shell, I, J
import numpy as np, json, hashlib
assert hashlib.sha256(open('data/d26_certificate.npz','rb').read()).hexdigest()=='35fd8143cf0c7e6b0146719485e00c86952cbd5e404c7881c5f6b92fa25b21ec'
w=json.load(open('data/d26_blocks.json'))
z=np.load('data/d26_certificate.npz',allow_pickle=False)
L=np.delete(leech_shell(),z['deleted_leech_indices'],axis=0)
A=np.vstack([np.pad(L,((0,0),(0,2))),z['added_integer']])
X=A[np.any(A[:,J],axis=1)]
assert len(X)==200864
D=np.unique(X[:,I+[24,25]] @ np.array(w['frame_rows']).T,axis=0)
assert len(D)==2601 and np.all(D%32==0)
json.dump((D//32).tolist(),open('replay42-P.json','w'))
```

```js
function verify(P,Q,B,T) {
 const req=(v,n)=>{if(!v)throw Error(n)};const range=n=>Array.from({length:n},(_,i)=>i); const pop=x=>{x>>>=0;x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
 let chambers=0,withM=0,K4=0,contactChecks=0;
 for(let b of B){let R=b.slice(0,4),a=[...b,17-b[4]];
  for(let i of R)req(Q.some(q=>q.includes(8)&&q.includes(9)&&q.includes(i)),1);
  for(let s=0;s<64;s++){const q=[1,1,1,1,0,2], seen=new Map;
   for(let p of P){let u=a.map((v,i)=>p[v]*((s>>i&1)?-1:1));if(u.reduce((t,v)=>t+Math.max(v,0),0)>4)seen.set(u.join(','),u);}
   req(seen.has(q.join(',')),2);req([...seen.values()].every(u=>u.every((v,i)=>v<=q[i])),3);chambers++;
  }
 }
 for(let U of Q){let O=range(10).filter(i=>!U.includes(i));
  let E=Q.filter(V=>V.filter(i=>U.includes(i)).length===2).map(V=>V.filter(i=>!U.includes(i)));
  if(E.some(a=>E.some(b=>E.some(c=>new Set([...a,...b,...c]).size===6)))){withM++;continue;}
  let J=O.filter(i=>i<8);req(J.length===4,4);for(let i of J)for(let j of J)if(i<j)req(E.some(e=>e.includes(i)&&e.includes(j)),5);
  for(let s=0;s<16;s++)for(let j of [8,9])for(let h of [-2,2]){req(P.some(p=>U.reduce((t,i,k)=>t+p[i]*((s>>k&1)?1:-1),0)===4&&O.every(i=>p[i]===(i===j?h:0))),6);contactChecks++;}
  K4++;
 }
 let masks=range(32).map(s=>T.reduce((a,t,i)=>a+((8*t.slice(0,5).reduce((d,v,k)=>d+((s>>k&1)?1:-1)*v,0)>9*t[5])?2**i:0),0)), adj=range(32).map(s=>range(5).reduce((a,i)=>a+2**(s^(1<<i)),2**s)>>>0);
 let n=0,best=-50,unique=0,nonempty=-99,hist={};
 function rec(c,k,ex){let val=k-pop(ex);hist[val]=(hist[val]||0)+1;n++;if(k)nonempty=Math.max(nonempty,val);if(val>best){best=val;unique=0;}if(val===best)unique++;while(c){let b=c&-c,s=31-Math.clz32(b);c=(c^b)>>>0;rec((c&~adj[s])>>>0,k+1,ex|masks[s]);}}
 rec(0xffffffff,0,0);
 req(chambers===256&&withM===10&&K4===4&&contactChecks===256&&n===254475&&best===0&&unique===1&&nonempty===-1,7);return {chambers,withM,K4,contactChecks,n,best,unique,nonempty,near:[-3,-2,-1,0].map(i=>hist[i])};
}
const fs=require('fs');
const P=JSON.parse(fs.readFileSync('replay42-P.json'));
const w=JSON.parse(fs.readFileSync('data/d26_blocks.json'));
console.log(verify(P,w.four_coordinate_supports,w.five_coordinate_supports,w.template_numerators_and_denominators));
```

Actual independent output: chambers 256, withM 10, K4 4, contactChecks 256, n 254475, best 0, unique 1, nonempty −1, near [12,4,1,1].
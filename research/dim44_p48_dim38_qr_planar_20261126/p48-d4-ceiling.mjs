
/* P48 D4 full-fibre ceiling, Lower57. Run under standard Node. 
   Pinned input is KissingNumbers main 26998bee3ff7338c043a885913d6e9a9cfacf4a0,
   verifications/improved/dim44-d4-section/families.npz. Inherit moment proof in moments.py. */
"use strict";
const loadNPZ=(function(){
const ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
function unb64(s){const v=[];let word=0,bits=0;for(let c of s){if(c==="=")break;let a=ABC.indexOf(c);if(a<0)continue;word=(word<<6)|a;bits+=6;if(bits>=8){bits-=8;v.push((word>>bits)&255);}}return Uint8Array.from(v);}
function inflate(a,n){
 let i=0,b=0,z=0,out=new Uint8Array(n),q=0;
 const bits=k=>{while(z<k){if(i>=a.length)throw Error("read early");b|=a[i++]<<z;z+=8;}let r=b&((1<<k)-1);b>>>=k;z-=k;return r;};
 function make(l){
  let count=new Int32Array(16),next=new Int32Array(16),arr=Array.from({length:16},()=>new Map());
  for(let x of l)if(x)count[x]++;
  for(let k=1;k<=15;k++)next[k]=(next[k-1]+count[k-1])*2;
  l.forEach((x,s)=>{if(x){let c=next[x]++,r=0;for(let k=0;k<x;k++,c>>=1)r=(r<<1)|(c&1);arr[x].set(r,s);}});
  return arr;
 }
 const read=T=>{let c=0;for(let k=1;k<=15;k++){c|=bits(1)<<(k-1);if(T[k].has(c))return T[k].get(c);}throw Error("bad code "+i);};
 const L0=[3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258];
 const D0=[1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577];
 let last=0;while(!last){
  last=bits(1);const ty=bits(2);let A,D;if(ty===0){bits(z%8);let s=bits(16),t=bits(16);if(s+t!==65535)throw Error("plain");for(let j=0;j<s;j++)out[q++]=bits(8);continue;}
  if(ty===1){A=make(Array.from({length:288},(_,i)=>i<144?8:i<256?9:i<280?7:8));D=make(Array(32).fill(5));}
  else if(ty===2){let hl=bits(5)+257,hd=bits(5)+1,hc=bits(4)+4,c=Array(19).fill(0);
   let ord=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15];for(let j=0;j<hc;j++)c[ord[j]]=bits(3);const C=make(c),l=[];
   while(l.length<hl+hd){let u=read(C);if(u<16){l.push(u);continue;}let rep=u===16?bits(2)+3:u===17?bits(3)+3:bits(7)+11, x=u===16?l.at(-1):0;for(let j=0;j<rep;j++)l.push(x);}
   A=make(l.slice(0,hl));D=make(l.slice(hl)); 
  }else throw Error("eof");
  while(true){let s=read(A);if(s===256)break;if(s<256){out[q++]=s;continue;}
   let l=s-257,e=l===28?0:Math.max(0,(l>>2)-1),k=L0[l]+bits(e),t=read(D),d=D0[t]+bits(Math.max(0,(t>>1)-1));
   for(let j=0;j<k;j++,q++)out[q]=out[q-d];
  }
 }
 if(q!==n)throw Error("length "+q+" "+n);return out;
}
const ascii=a=>Array.from(a,x=>String.fromCharCode(x)).join("");
const u16=(a,i)=>a[i]+256*a[i+1],u32=(a,i)=>(a[i]+256*a[i+1]+65536*a[i+2]+16777216*a[i+3])>>>0;
function zip(a){
 let Z={};for(let p=0;p<a.length;){
  let m=u32(a,p);if(m!==0x04034b50)break;
  let t=u16(a,p+8),dry=u32(a,p+18),wet=u32(a,p+22),ns=u16(a,p+26),ex=u16(a,p+28);
  let name=ascii(a.slice(p+30,p+30+ns)).replace(/\.npy$/,"");
  for(let j=p+30+ns;j<p+30+ns+ex;){let id=u16(a,j),len=u16(a,j+2);if(id===1){wet=u32(a,j+4);dry=u32(a,j+12);}j+=4+len;}
  p+=30+ns+ex;Z[name]=t===8?inflate(a.slice(p,p+dry),wet):a.slice(p,p+dry);p+=dry;
 }
 return Z;
}
function npy(a){
 let k=a[6]===1?10:12,n=a[6]===1?u16(a,8):u32(a,8),s=ascii(a.slice(k,k+n));
 let shape=s.match(/'shape': *\(([^)]*)\)/)[1].split(",").map(x=>x.trim()).filter(x=>x).map(Number),D=a.slice(k+n),u=s.match(/'descr': *'([^']+)'/)[1];
 let d=new DataView(D.buffer,D.byteOffset,D.byteLength),nb=Number(u.slice(2)),N=shape.reduce((x,y)=>x*y,1), vals=Array.from({length:N},(_,i)=>nb===8?Number(d.getBigInt64(i*nb,true)):nb===4?d.getInt32(i*nb,true):nb===2?d.getInt16(i*nb,true):d.getInt8(i*nb));
 return {shape,dtype:u,vals};
}
return s=>Object.fromEntries(Object.entries(zip(unb64(s))).map(([k,a])=>[k,npy(a)]));
})();
const F="UEsDBC0AAAAIAAAAIQCYvkkF//////////8RABQAYWZmaW5lX2NvdW50cy5u\ncHkBABAAwDMAAAAAAABTAwAAAAAAAO1YPWtUQRR90USCENBoYblgsUa2CIqw\nihgL0cIoIkkhiGYx6weIkV0RRIQUIvY2linszD+wsfAHRLDYysrCUmwU2eya\n996cgMd9zJlleLd5t7nMzOHOnDtvztx5b64uX7l2YyJ5mjyvr7a7dzr1M7X6\n2QfNeqNWv7vWedJpPbq91lltp/0XWw+77Z3+7v3W4/ZO+9jJE81G7dRco/ai\nNqbtTy4vJKmdT5KRHvZhaSrzR1ybPWxr5lDm51y7N8wNHrjXB3M/5dp9N94n\nnBpPXd/69dxPJ6N9KE7lAZvgDjLwBT/mzfNiHp4/NB7jGA9rrOQMfrpx+JlI\nOIyj7UnXriGvajzOH+cRdtgF3i5Yj+seqjgrvlb7a3UusT7OB+cF+zb0mBpP\nXd/pvbnfk4z2obiE+nn9+F5C9coXz4qv1f6qustWpPuqDnGcAe0Dx/fhkBce\nRxum8lXjFdD6b/1qXlR9UXHq96fWJep3peKs+IaeD5++qLjYuhH7/ojNV63X\n1Hs/dn2g4tTza8U3dp2j1n/j3h9svncE46z4VvpcDt9Kn8vhO+77suicq7jY\n5xfmqxOt+Kr3h4oL1TVfHRv7PrLia3V+1fOmvt9UnBVfmO89Hfs+iv1fRb23\nrPhW+lwO30qfy+H79tKze6n/+CsfZw/c+va7A6mfdwU5e+C2BkvHU9/7k/f/\nnsz74WHTs58y/8qta9H1L/4Lk+Op69vsLWSS0XTzsg/FqTw2buX78c3hfpBn\nvuDHvHlezMPzh8ZjHOOBax69sC/1L93BgH8/iIPDONq76/O8479gf8R4nD/O\nI3ArX2czf5Pmw/yNQJwVX6v9tTqXWB/ng/OCfVt27UnyvG++eOr6NvrnMv/d\n9bMPxeF7+ez6mQe+F1l3xXhWfK32V9Vd6D30nXW/KF6RDmE/Nl3/vOvH/zHs\nh4pDXngcbeRF5avG47xwfkLzouqLipPPm1iXqN+VirPiG3o+fPoi102RdSP2\n/RGbr1qvqfd+7PpAxcm6ZsQ3dp2j1n9y3S7Wf+PWiWXxrfS5HL6VPpfDd9z3\nZdE5V3Gxz69aJ1rxVe8PFReqa746NvZ9ZMXX6vyq5019v6k4K77qezr2fRT7\nv4r8H9iIb6XP5fCt9Lkcvn8BUEsDBC0AAAAIAAAAIQBnSYCj//////////8K\nABQAbGFiZWxzLm5weQEAEACgAgAAAAAAAJcAAAAAAAAAm+wX6hsQychQxlCt\nnpJanFykbqWgbpNpoa6joJ6WX1RSlJgXn1+UkgoSd0vMKU4FihdnJBakAvka\nhuY6CiaaOgq1CuQCLgYoyMFBw8AiNHEzNHl08UVo8ujiuMw3wyGfg0Mduvno\n6tD1o9M30PhmaPQNND66OnR5dHXo9ixCo9Hdi64OXR5dHQOaOnR3MKDxccUf\nLvcDAFBLAwQtAAAACAAAACEAKrpPtf//////////CQAUAHJhZGlpLm5weQEA\nEAAIAQAAAAAAAFYAAAAAAAAAm+wX6hsQychQxlCtnpJanFykbqWgbpNpoa6j\noJ6WX1RSlJgXn1+UkgoSd0vMKU4FihdnJBakAvkahuY6mjoKtQrkA64Z7Axg\nsISFNFqAFUKHsFGHBgBQSwMELQAAAAgAAAAhAPHydaH//////////xAAFABj\nYWxpYnJhdGlvbnMubnB5AQAQAAABAAAAAAAAXAAAAAAAAACb7BfqGxDJyFDG\nUK2eklqcXKRupaBuk2mhrqOgnpZfVFKUmBefX5SSChJ3S8wpTgWKF2ckFqQC\n+RomOgommjoKtQpkA66//yGAAQqY0Whc8gw45HGZg67+H5QGAFBLAwQtAAAA\nCAAAACEAVYR63P//////////DwAUAHBvcHVsYXRpb25zLm5weQEAEAAIAQAA\nAAAAAFsAAAAAAAAAm+wX6hsQychQxlCtnpJanFykbqWgbpNpoa6joJ6WX1RS\nlJgXn1+UkgoSd0vMKU4FihdnJBakAvkahuY6mjoKtQrkA64Ff+wYQKCPgQEr\n3YaDtmSG0DoJjFShAVBLAQItAy0AAAAIAAAAIQCYvkkFUwMAAMAzAAARAAAA\nAAAAAAAAAACAAQAAAABhZmZpbmVfY291bnRzLm5weVBLAQItAy0AAAAIAAAA\nIQBnSYCjlwAAAKACAAAKAAAAAAAAAAAAAACAAZYDAABsYWJlbHMubnB5UEsB\nAi0DLQAAAAgAAAAhACq6T7VWAAAACAEAAAkAAAAAAAAAAAAAAIABaQQAAHJh\nZGlpLm5weVBLAQItAy0AAAAIAAAAIQDx8nWhXAAAAAABAAAQAAAAAAAAAAAA\nAACAAfoEAABjYWxpYnJhdGlvbnMubnB5UEsBAi0DLQAAAAgAAAAhAFWEetxb\nAAAACAEAAA8AAAAAAAAAAAAAAIABmAUAAHBvcHVsYXRpb25zLm5weVBLBQYA\nAAAABQAFACkBAAA0BgAAAAA=\n",data=loadNPZ(F);
const model=(function(data,cal){
const R=324,Adj=[[216,108,108,108],[108,108,54,54],[108,54,108,54],[108,54,54,108]];
const mod=x=>(x%R+R)%R;
function row(v){return Adj.map(b=>v.reduce((s,x,i)=>s+x*b[i],0));}
function q(v,a){return v.reduce((s,x,i)=>s+x*a[i],0);}
const idx=new Map(), dual=[];
for(let a=-6;a<=6;a++)for(let b=-6;b<=6;b++)for(let c=-6;c<=6;c++)for(let d=-6;d<=6;d++){
 let v=[a,b,c,d],w=row(v),s=q(v,w),L=w.map(mod);if(s>8*R)continue;
 let key=L.join();if(!idx.has(key)||idx.get(key).s>s)idx.set(key,{v,L,s,w});
 dual.push({v,key,s});
}
const vals=[...idx.values()].sort((a,b)=>{for(let i=0;i<4;i++)if(a.L[i]!==b.L[i])return a.L[i]-b.L[i];return 0});
const D=new Map(vals.map((v,i)=>[v.L.join(),i]));const bins=[];let orb=0;
for(let i=0;i<vals.length;i++){
 let u=vals[i],j=D.get(u.L.map(v=>mod(-v)).join());
 u.orb=j<i?vals[j].orb:orb++;
}
const forms=Array.from({length:328},(_,i)=>data.affine_counts.vals.slice(i*5,i*5+5));
function count(v,s,c=cal){const r=forms[v.orb+(s===6?0:164)];return r[0]+c.reduce((s,v,j)=>s+v*r[j+1],0);}
const V=[];
for(let u of vals)for(let s of [6,8]){let r=s*R-u.s;if(count(u,s)>0&&r<=6*R)V.push({...u,r,k:s,c:count(u,s)});}
const origin=data.labels.vals.length/4;
const prev=Array.from({length:origin},(_,i)=>{
 let L=data.labels.vals.slice(i*4,i*4+4), r=data.radii.vals[i];
 let j=V.findIndex(v=>v.L.join()===L.join()&&v.r===r);
 
 return j;
}).filter(j=>j>=0);
function edges(v,mode=2, calibrations){
  const w=V[v];return V.map((u,j)=>{
   if(j===v)return 1;
   const L=u.L.map((x,i)=>mod(x-w.L[i])),z=vals[D.get(L.join())];let min=6;
   if(mode>=1&&count(z,6,calibrations)===0){min=8;if(mode>=2&&count(z,8,calibrations)===0)min=10;}
   let r=u.r+w.r-min*R+z.s;
   return +(r<=0||r*r<=u.r*w.r);
  });
}
const E=[null,null,V.map((_,i)=>edges(i,2))];
const hist={};
for(let u of vals){let k=String(u.s/R);hist[k]??=[0,0,0];hist[k][0]++;if(count(u,6)===0)hist[k][1]++;if(count(u,6)===0&&count(u,8)===0)hist[k][2]++;}
return {V,vals,forms,E,hist,prev,dual,Dkeys:[...D.keys()]};
});
const poly=(function(full){
const gcd=(a,b)=>{if(a<0n)a=-a;if(b<0n)b=-b;while(b){let r=a%b;a=b;b=r;}return a;};
function det(a){if(a.length===0)return 1n;let s=0n;for(let j=0;j<a.length;j++)s+=(j%2?-1n:1n)*a[0][j]*det(a.slice(1).map(r=>r.filter((_,k)=>k!==j)));return s;}
const unique=new Map(),A=full.filter(r=>r[0]+276*Math.min(0,...r.slice(1))<=0).map(r=>r.map(BigInt));
function dfs(i,P){if(P.length<4){for(let j=i;j<=A.length-(4-P.length);j++)dfs(j+1,P.concat([A[j]]));return;}
const mat=P.map(r=>r.slice(1)), rhs=P.map(r=>-r[0]);let d=det(mat);if(d===0n)return;
let ns=mat[0].map((_,i)=>det(mat.map((r,j)=>r.map((x,k)=>i===k?rhs[j]:x))));
if(d<0n){d=-d;ns=ns.map(x=>-x);}
const g=ns.reduce(gcd,d);d/=g;ns=ns.map(x=>x/g);
if(A.some(row=>row.slice(1).reduce((s,a,j)=>s+a*ns[j],row[0]*d)<0n))return;
unique.set(d+":"+ns.join(),[d,...ns].map(Number));
}
dfs(0,[]);
return [...unique.values()];
});
const run=(function(m,mode){
const V=m.V, E=m.E[mode],n=V.length;
const w=V.map(v=>v.c);
const degree=E.map(r=>r.reduce((x,y)=>x+y,0));
const desc=V.map((_,i)=>i).sort((i,j)=>w[j]-w[i]||degree[j]-degree[i]||i-j);
let b=m.prev, best=b.reduce((s,i)=>s+w[i],0),nodes=0,tests=0, cmp=0, done=true;
function pack(P){
 const C=[], fresh=P.slice().sort((i,j)=>w[j]-w[i]||degree[j]-degree[i]||i-j);
 for(let i of fresh){let k=C.findIndex(a=>a.every(j=>!E[i][j]));if(k<0)C.push([i]);else C[k].push(i);}
 const cc=[], ub=[];let mass=0;
 for(let K of C){let cur=0;
  K=K.slice().sort((i,j)=>w[i]-w[j]||i-j);
  for(let i of K){cur=Math.max(cur,w[i]);ub.push(mass+cur);cc.push(i);}
  mass+=cur;
 }
 return [cc,ub];
}
function ac(P,C,N){
 if(nodes++>=1e6){done=false;return;}
 if(!P.length){if(N>best){best=N;b=C;}return;}
 let [ord,ub]=pack(P);
 for(let i=ord.length-1;i>=0;i--){
  if(N+ub[i]<=best)return;
  let v=ord[i],newN=N+w[v]; 
  let Q=[];for(let j=0;j<i;j++)if(E[v][ord[j]])Q.push(ord[j]);
  ac(Q,[...C,v],newN);
  if(!done)return;
 }
}
for(let i of desc){let chosen=[i],n=w[i];
for(let j of desc)if(j!==i&&chosen.every(z=>E[z][j])){chosen.push(j);n+=w[j];}
if(n>best){best=n;b=chosen;}}
const greedy=best;
let anchors=V.map((_,i)=>i).filter(i=>V[i].r===1944);
const anchored=[];
for(let a of anchors){
 let P=desc.filter(i=>i!==a&&E[a][i]);let prior=best;
 ac(P,[a],w[a]);anchored.push([a,best,best-prior,nodes]);
}
ac(desc.filter(i=>V[i].r!==1944),[],0);
return {best,b,greedy,nodes,done,anchored,scope:"Exact finite universal P48 D4 projected shells from 6/8, positive nonzero radius<=6; not all actual endpoint conflict possibilities or other embeddings."};
});
const scope="Complete norm6/8 projected fibres, radius squared≤6, and shortest-distance certificate using exact empty source-shell bins; not an upper bound for unrestricted 44D configurations or endpoint-specific subfamilies.";
const demand=(t,m)=>{if(!t)throw Error(m);};
const orig=[0,0,0,134],M=model(data,orig),F0=[...new Map(M.forms.map(r=>[r.join(),r])).values()];
demand(F0.some(r=>r.join()==="276,-1,-1,-1,-1"),"simplex cap");
for(let i=1;i<=4;i++)demand(F0.some(r=>r[0]===0&&r.every((v,j)=>v===+(j===i))),"positive calibration coords");
const P=poly(F0),out=[];
demand(P.length===12&&P.every(x=>x[0]===1),"12 integral vertices");
demand(M.vals.length===324&&M.V.length===404,"complete source classes");
const o=run(M,2);demand(o.done&&o.best===4940961,"fixed actual frame");
for(let [d,...v] of P){
 const m=model(data,v);
 demand(m.prev.every(i=>m.prev.every(j=>m.E[2][i][j])),"retained old construction remains a legal lower incumbent");
 const r=run(m,2);demand(r.done,"incomplete exact branch-and-bound");
 out.push([v,r.best,r.nodes]);
}
demand(Math.max(...out.map(x=>x[1]))===4940961,"global envelope");
const ret={bound:4940961,vertices:out,fixed_frame:[orig,o.best,o.nodes],scope};
if(typeof report==="function")report(ret);else console.log(JSON.stringify(ret,null,2));


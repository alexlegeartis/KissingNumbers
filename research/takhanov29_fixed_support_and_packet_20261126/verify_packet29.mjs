import {readFileSync} from "node:fs";
import {createHash} from "node:crypto";
import assert from "node:assert/strict";
const input=readFileSync(process.argv[2]);
assert.equal(createHash("sha256").update(input).digest("hex"),"b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824");
function uz(buf) {let s=0,o=[]; const bits=n=>{let a=0;for(let i=0;i<n;i++)a|=((buf[s>>3]>>(s++&7))&1)<<i;return a};function tree(ls){let c=new Array(16).fill(0);for(let l of ls)if(l)c[l]++;let n=new Array(16).fill(0);for(let i=1;i<16;i++)n[i]=(n[i-1]+c[i-1])<<1;let m=new Map();ls.forEach((l,i)=>{if(l)m.set((1<<l)|n[l]++,i)});return ()=>{let h=1;for(let i=1;i<16;i++){h=(h<<1)|bits(1);if(m.has(h))return m.get(h)}throw Error('huff')}}
const LB=[3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258],DB=[1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577];
let go=true;while(go){go=!bits(1);let tp=bits(2);if(!tp){s=(s+7)&~7;let l=bits(16);bits(16);while(l--)o.push(bits(8));continue} let L,D;
if(tp===1){L=tree(Array.from({length:288},(_,i)=>i<144?8:i<256?9:i<280?7:8));D=tree(Array(32).fill(5))}
else if(tp===2){let hl=bits(5)+257,hd=bits(5)+1,hc=bits(4)+4, cl=Array(19).fill(0), pl=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15];for(let j=0;j<hc;j++)cl[pl[j]]=bits(3);let de=tree(cl),ll=[];
while(ll.length<hl+hd){let a=de();if(a<16)ll.push(a);else if(a===16)ll.push(...Array(3+bits(2)).fill(ll.at(-1)));else if(a===17)ll.push(...Array(3+bits(3)).fill(0));else ll.push(...Array(11+bits(7)).fill(0));}
L=tree(ll.slice(0,hl));D=tree(ll.slice(hl));}else throw Error('block type');
for(;;){let a=L();if(a===256)break;if(a<256){o.push(a);continue}let k=a-257,len=LB[k]+bits(k<8||k===28?0:(k-4)>>2),d=D(),dist=DB[d]+bits(d<4?0:(d-2)>>1);while(len--)o.push(o[o.length-dist]);if(o.length>5000000)throw Error('len');}}
return Uint8Array.from(o)}
function npz64(z){let bytes=Buffer.from(z,'base64');let v=new DataView(bytes.buffer,bytes.byteOffset,bytes.byteLength),p=0,o={};while(v.getUint32(p,true)===0x04034b50){let nl=v.getUint16(p+26,true),el=v.getUint16(p+28,true),n=String.fromCharCode(...bytes.slice(p+30,p+30+nl));let c=v.getUint32(p+18,true);if(c===0xffffffff)c=Number(v.getBigUint64(p+30+nl+12,true));let t=p+30+nl+el, d=uz(bytes.slice(t,t+c)), h=d[8]|d[9]<<8, head=String.fromCharCode(...d.slice(10,h+10));o[n]=[head,d.slice(10+h)];p=t+c;}return o}

const NP=npz64(input.toString('base64'));
const view=k=>new DataView(NP[k+".npy"][1].buffer), aj=view("added_integer"), ix=view("midpoint_rows"), di=view("deleted_leech_indices");
const ADA=Array.from({length:67868},(_,i)=>Int16Array.from({length:29},(_,j)=>aj.getInt16(2*(29*i+j),true)));
const MID=new Set(Array.from({length:64905},(_,i)=>ix.getInt32(4*i,true))), BODY=new Map(), OTHER=[];
const bdot=(a,b,at=0,t=24)=>{let d=0;for(let j=at;j<t;j++)d+=a[j]*b[j];return d;}
const kt=(t)=>Array.from(t).join(',');
for(let i=0;i<ADA.length;i++){const r=ADA[i];if(!MID.has(i)){OTHER.push(r);continue;}const b=r.slice(0,24),t=r.slice(24),s=kt(b);let o=BODY.get(s);if(o)o.ts.push(t);else BODY.set(s,{b,ts:[t]});}
const DB=[...BODY.values()].filter(x=>x.ts.length==2), SB=[...BODY.values()].filter(x=>x.ts.length==1);
const rootpal=[];for(let i=0;i<5;i++)for(let j=i+1;j<5;j++)for(let s of [1,-1]){let v=new Int16Array(5);v[i]=4;v[j]=4*s;rootpal.push(v);}
function canon(x){let i=x.findIndex(z=>z);return x[i]>=0?kt(x):kt(x.map(z=>-z));}
const pm=new Map(rootpal.map((x,i)=>[kt(x),i]));for(let o of DB.concat(SB))o.k=pm.get(canon(o.ts[0]));


const I6=[31061,32169,31680,31719,30909,32319], S=[...BODY.values()], six=I6.map(i=>S[i]), sc=new Set(six.map(x=>kt(x.b)));
const phases=new Map();for(let r of ADA){if(sc.has(kt(r.slice(0,24))))continue;let t=r.slice(24),k=kt(t);let o=phases.get(k);
 if(!o){o={t,v:Int16Array.from({length:6},()=>300)};phases.set(k,o);}
 for(let i=0;i<6;i++){let z=64-bdot(r,six[i].b);if(z<o.v[i])o.v[i]=z}
}
const W6=[...phases.values()], nA=W6.length;
const T6=[...W6].filter(q=>q.v.some(x=>x<32));
function allowed(v,n,i){for(let q of W6){let a=q.v[i],u=bdot(v,q.t,0,5);
  if(a>=0){if(u>0&&32*u*u>a*a*n)return false;}
  else if(u>=0||32*u*u<a*a*n)return false;
 }return true;}

function leech(){
 const bit=(x)=>{x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0xf0f0f0f)*0x1010101)>>>24};
 let Y=[],G=[0];for(let i=0;i<12;i++){let m=[0,2,4,5,6,10,11].reduce((x,j)=>x|1<<(i+j)%23,1<<23);G=G.concat(G.map(x=>x^m))}
 for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(let a of [-8,8])for(let b of [-8,8]){let x=new Int8Array(24);x[i]=a;x[j]=b;Y.push(x)}
 for(let q of G){
   if(bit(q)==8){let js=Array.from({length:24},(_,i)=>i).filter(i=>q>>i&1);
      for(let s=0;s<256;s++)if(!(bit(s)&1)){let x=new Int8Array(24);js.forEach((i,k)=>x[i]=s>>k&1?4:-4);Y.push(x)}}
   for(let k=0;k<24;k++){let x=Int8Array.from({length:24},(_,j)=>q>>j&1?2:-2);x[k]*=-3;Y.push(x)}
 }
 Y.sort((a,b)=>{for(let j=0;j<24;j++)if(a[j]!==b[j])return a[j]-b[j];return 0});if(Y.length!==196560)throw Error('shell');return Y
}
const m=Array.from(six[0].b,(a,i)=>(a+six[1].b[i])/2);
const V=[0,2,4].map(j=>Array.from(six[j].b,(a,i)=>a-m[i]));
const del=new Set(Array.from({length:32198},(_,i)=>di.getInt32(4*i,true)));
const EA=leech().filter((_,i)=>!del.has(i));
const C6=new Map(), C8=new Map();
const oldA=new Set(six.slice(0,4).map(x=>kt(x.b))),old6=new Set(six.map(x=>kt(x.b)));
const addw=(C,u,e)=>{let k=kt(u),o=C.get(k);if(o===undefined||e<o)C.set(k,e)};
function doW(r,eq=false){
 let a=64-bdot(m,r), p=V.map(v=>bdot(v,r)), q=eq?Array(5).fill(0):Array.from(r.slice(24),x=>4*x); const k=eq?'':kt(r.slice(0,24));
 if(!oldA.has(k))addw(C6,[p[0],p[1],q[0],q[1]+q[2],q[3],q[4]],a);
 if(!old6.has(k))addw(C8,[...p,...q],a);
}
for(let r of EA)doW(r,true);for(let r of ADA)doW(r);
function red(C){const gcd=(a,b)=>b?gcd(b,a%b):Math.abs(a);const A=new Map();let b=0;for(let [r,c] of C){let z=[r[0],r[1],r[2]-r[5],r[3],r[4]],g=z.reduce(gcd,c);if(!g)continue; let q=z.map(x=>x/g);c/=g;if(q.every(x=>!x))continue; if(2*(q.reduce((a,b)=>a+b*b,0)+q[4]*q[4])<=c*c)continue;let ky=q.join(',');let cur=A.get(ky);if(cur===undefined||c<cur)A.set(ky,c)}return [...A].map(([k,b])=>[k.split(',').map(Number),b]);}
const D=(A)=>{if(A.length==0)return 1;let b=0;for(let j=0;j<A.length;j++)if(A[0][j]!==0)b+=(j%2?-1:1)*A[0][j]*D(A.slice(1).map(v=>v.filter((_,i)=>i!=j)));return b;};
function choose(n,k){let V=[];function f(a,j){if(a.length==k){V.push(a);return}for(let i=j;i<=n-k+a.length;i++)f(a.concat([i]),i+1)}f([],0);return V}
function verts(P){
 const n=P[0].length-1,L=[],S=new Set();
 for(let zz of choose(P.length,n)){let v=zz.map(i=>P[i]), b=v.map(r=>r[n]), A=v.map(r=>r.slice(0,n)), q=D(A);if(!q)continue;
 let a=Array.from({length:n},(_,k)=>D(A.map((r,j)=>r.map((x,i)=>i==k?b[j]:x))));if(q<0){q=-q;a=a.map(t=>-t)}
 if(P.some(r=>a.reduce((s,x,i)=>s+x*r[i],0)>q*r[n]))continue;
 let ky=a.map(x=>x/q).join(",");if(!S.has(ky)){S.add(ky);L.push([a,q]);}
 }return L;
}
const B=[[1,1,0,0,0,1],[0,0,1,1,0,1],[1,0,1,0,0,1],[0,1,1,0,0,1],[0,0,0,0,1,1],...Array.from({length:5},(_,i)=>Array.from({length:6},(_,j)=>i==j?-1:0))];
const v29=[...Array(2)].map((_,i)=>{const P=B.concat([i?[1,1,0,1,1,2]:[0,0,1,2,1,2]]);let vs=verts(P), V=vs.map(([a,q])=>[2*a.reduce((t,v)=>t+v*v,0)-a[4]*a[4],q*q,a.map(x=>x/q)]);return {N:V.length,top:V.sort((a,b)=>b[0]/b[1]-a[0]/a[1]).slice(0,8),bad:V.filter(a=>a[0]>4*a[1]),eq:V.filter(a=>a[0]==4*a[1]).map(a=>a[2])}});

assert.equal(EA.length,164362);assert.equal(C6.size,576);assert.equal(C8.size,1426);
assert.equal(ADA.filter(r=>oldA.has(kt(r.slice(0,24)))).length,4);
assert.equal(ADA.filter(r=>old6.has(kt(r.slice(0,24)))).length,8);
assert.equal(bdot(m,m),64);
for(let i=0;i<3;i++){assert.equal(bdot(V[i],m),0);for(let j=0;j<3;j++)assert.equal(bdot(V[i],V[j]),i==j?32:0);}
for(let i=0;i<6;i++)assert.deepEqual(Array.from(six[i].b),m.map((a,j)=>a+(i%2?-1:1)*V[i>>1][j]));
const tt=(r,val)=>{let u=phases.get(kt(r));assert.ok(u);for(let x of u.v)assert.ok(x<=val)};
for(let s of [-1,1]){
 for(let t of [-1,1]){
  tt([4*s,0,0,0,4*t],0);
  tt([0,0,4*s,4*t,0],16);
 }
 for(let j of [1,2,3]){let v=[0,0,0,0,0];v[j]=8*s;tt(v,32);}
 tt([0,4*s,0,4*s,0],16);
}
const b6 = [
 [0,4,4,0,0],[0,-4,-4,0,0],[0,4,-4,0,0],[0,-4,4,0,0],[0,4,0,-4,0],[0,-4,0,4,0]
];for(let r of b6)for(let i=0;i<6;i++)assert.ok(allowed(r,32,i));
for(let s of [16,-16])assert.equal(C6.get([0,0,s,0,0,s].join(',')),0);
const raw=[...C6].map(([k,c])=>[k.split(',').map(Number),c]);
const min=red(raw), nn=new Set(min.map(([r,c])=>r.join(',')+':'+c));
const req=(r,c)=>assert.ok(nn.has(r.join(',')+':'+c),r.join(',')+':'+c);
for(let a of [-1,1])for(let b of [-1,1]){
 req([a,b,0,0,0],1);req([0,0,a,b,0],1);req([a,0,b,0,0],1);req([0,a,b,0,0],1);
 for(let c of [-1,1]){req([0,0,a,2,c],2);req([a,b,0,-1,c],2);}
}
for(let a of [-1,1])req([0,0,0,0,a],1);
for(let o of v29){assert.equal(o.N,19);assert.equal(o.bad.length,0);assert.deepEqual(o.eq.map(row=>row.map(x=>x||0)).sort(),[[0,1,0,1,0],[1,0,0,1,0]]);}
if(process.argv[3]){
 let txt=readFileSync(process.argv[3],'utf8').trim().split(/\r?\n/).slice(1);
 const norm=[64n,32n,32n,32n,16n,16n,16n,16n,16n];
 const Q=txt.map(l=>l.trim().split(/\s+/).map((s,i)=>BigInt(Math.round(Number(s)*1e10/Math.sqrt(Number(norm[i]))))));
 assert.equal(Q.length,9);assert.ok(Q.every(r=>r.length==9&&r[0]>=0n));
 const dot=(a,b)=>a.reduce((s,g,i)=>s+g*b[i],0n);
 const W=[...C8].map(([p,c])=>[64-c,...p.split(',').map(Number)].map(BigInt));
 const ns=Q.map(r=>dot(r,r.map((x,i)=>x*norm[i])));
 for(let i=0;i<Q.length;i++){assert.ok(ns[i]>0n);for(let r of W){let t=dot(r,Q[i]);assert.ok(t<=0n||t*t<=32n*ns[i],'external inequality');}
  for(let j=0;j<i;j++){let t=dot(Q[i],Q[j].map((x,k)=>x*norm[k]));assert.ok(t<=0n||4n*t*t<=ns[i]*ns[j],'new-new inequality');}}
 console.log(JSON.stringify({verified:true,externalFrozen:232222,count:232231,strict:true,integerRays:Q.map(r=>r.map(String))}));
}
console.log(JSON.stringify({fixedBodyContinuousMax:8,keepBTwiceAxesMax:4,keptShell:EA.length,c6Normals:C6.size,c8Normals:C8.size,positiveCandidateChecked:!!process.argv[3]}));


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

const DS=[...BODY.values()];
const othp=OTHER.map(t=>({b:t.slice(0,24),d:Int16Array.from(rootpal,r=>Math.abs(bdot(r,t.slice(24),0,5)))}));
const dd=(a,b)=>{let d=0;for(let j=0;j<24;j++)d+=a[j]*b[j];return d;}
const all=(1<<20)-1;
const ODom=DS.map(o=>{
 let m=all;
 for(let a of othp){
  if(!m)break;
  let d=dd(o.b,a.b);
  if(d<=32)continue;let s=64-d;
  for(let ss=m;ss;ss&=ss-1){let j=31-Math.clz32(ss&-ss);if(a.d[j]>s)m&=~(1<<j);}
 }
 return m;
});

const Dom=ODom; DS.forEach((x,i)=>x.i=i);
const dot24=(a,b)=>{let d=0;for(let j=0;j<24;j++)d+=a[j]*b[j];return d;};
const RB=rootpal.map(r=>rootpal.map(s=>Math.abs(bdot(r,s,0,5))));const DD=DB.map(x=>[x.b,x.k,x.i]);const SS=SB.map(x=>[x.b,x.k,x.i, x.ts[0]]);
let Q=[];
for(let [b,k,i,s0] of SS) {
 if(!Dom[i])continue;
 const N=[];
 for(let [a,k1,j] of DD){let d=dot24(a,b);if(d>32)N.push([j,k1,64-d]);}
 let m=Dom[i];for(let [a,k1,j,t] of SS)if(j!=i){let p=64-dot24(a,b);if(p<32)for(let ms=m;ms;ms&=ms-1){let z=31-Math.clz32(ms&-ms);if(Math.abs(bdot(t,rootpal[z],0,5))>p)m&=~(1<<z)}}
 for(let ms=m;ms;ms&=ms-1){let z=31-Math.clz32(ms&-ms);let blockers=N.filter(([j,k1,p])=>RB[z][k1]>p).map(([j])=>j);let locked=blockers.filter(j=>(Dom[j]&(Dom[j]-1))==0).length;Q.push({i,z,fr:blockers,locked,N:N.length});}
}

const direct=Q.filter(r=>r.fr.every(j=>(DS[j].k>>1)==(r.z>>1))), par=new Int32Array(DS.length).map((_,i)=>i), size=new Int32Array(DS.length).fill(1);
DS.forEach((o,i)=>o.i=i);
function find(x){while(par[x]!=x){par[x]=par[par[x]];x=par[x];}return x}
function union(a,b){let x=find(a),y=find(b);if(x===y)return;if(size[x]<size[y])[x,y]=[y,x];par[y]=x;size[x]+=size[y]}
const col=Array.from({length:20},()=>[]), BC= new Map;
for(let o of DB)col[o.k].push(o);
let edge=0;
for(let s=0;s<20;s+=2)for(let a of col[s])for(let b of col[s+1])if(dd(a.b,b.b)>32){edge++;union(a.i,b.i);}
let aux=0,pin=0;for(let o of DB){let key=find(o.i),f=BC.get(key);if(!f){f={a:false,p:new Set,hold:[],s:o.k>>1};BC.set(key,f)}f.hold.push(o.i);
 if(!(Dom[o.i]>>(o.k^1)&1)){f.a=true;aux++;}
 for(let x of SB)if((x.k^1)==o.k&&dd(o.b,x.b)>32){f.p.add(x.i);pin++;}
}
const count=Array(10).fill(0), big=Array(10).fill(0);for(let o of BC.values()){count[o.s]++;big[o.s]=Math.max(big[o.s],o.hold.length);}
let H={},WIN=[],wit=[];
for(let q of direct){
 let N=col[q.z].concat(col[q.z^1]).filter(b=>dd(b.b,DS[q.i].b)>32);
 let M=new Map(),bad=[];
 for(let o of N){let c=find(o.i),f=+(o.k==q.z);if(M.has(c)&&M.get(c)!==f)bad.push('phase');M.set(c,f);}
 for(let [c,f] of M)if(f){let r=BC.get(c);if(r.a)bad.push('aux');if([...r.p].some(i=>i!=q.i))bad.push('single');}
 const k=[...new Set(bad)].sort().join('|')||'win';H[k]=(H[k]||0)+1;
 if(k=='win')WIN.push({q,changes:[...M].filter(([,f])=>f).flatMap(([i])=>BC.get(i).hold)});
 if(q.i==31061&&q.z==9){let L=new Map();for(let o of N){const c=find(o.i);if(!L.has(c))L.set(c,[]);L.get(c).push([o.i,+(o.k==q.z)]);}wit=[...L].map(([i,l])=>[i,l,BC.get(i).hold.length,BC.get(i).a,[...BC.get(i).p]]);}
}

const CB=[...BC.keys()], cmap=new Map(CB.map((x,i)=>[x,i]));
const NS=SB.length;
const dot5=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]+a[4]*b[4];
const rb=rootpal.map(r=>rootpal.map(s=>dot5(r,s)));
const d29domain=SB.map((s,j)=>{
 let N=DB.map(o=>[o,dd(o.b,s.b)]).filter(x=>x[1]>32);let excl=OTHER.map(o=>[o,dd(o,s.b)]).filter(x=>x[1]>32);
 let modes=[];
 for(let k=0;k<20;k++){
  let p=new Map(),bad=false;
  for(let [o,b] of N){
   if((k>>1)!=(o.k>>1)){if(Math.abs(rb[k][o.k])+b>64){bad=true;break}continue;}
   if(b>64){bad=true;break;}
   let a=find(o.i),f=k==o.k?1:0;if(p.has(a)&&p.get(a)!=f||f&&BC.get(a).a){bad=true;break;}
   p.set(a,f);
  }
  if(bad)continue;
  const pm=[...p].map(([a,b])=>[cmap.get(a),b]);
  const t=rootpal[k];
  const allowSign=[-1,1].filter(z=>excl.every(([o,b])=>b+z*dot5(t,o.slice(24))<=64));
  for(let sg of allowSign){modes.push({k,sg,pm,old:k==s.k&&sg*dot5(t,s.ts[0])>0});}
  if(allowSign.length==2){modes.push({k,sg:0,pm,old:false});}
 }
 return modes;
});
const Iad=SB.map(o=>o.i);
const E31=SB.map((s,j)=>SB.map(u=>64-dd(s.b,u.b)));
const Flags=CB.map(c=>BC.get(c).a?0:-1);
const U={all:d29domain,E31,rb,target:d29domain.flatMap((S,i)=>S.flatMap((q,j)=>q.sg==0?[[i,j]]:[]))}, D=U.all, N=D.length;
function match(i,a,j,b){let x=D[i][a],y=D[j][b],d=U.rb[x.k][y.k];if((!x.sg||!y.sg?Math.abs(d):x.sg*y.sg*d)>U.E31[i][j])return false;
for(let [k,v] of x.pm)for(let [l,w] of y.pm)if(k==l&&v!=w)return false;
return true}
const Adj=Array.from({length:N},()=>[]), G=Array.from({length:N},()=>[]);
for(let i=0;i<N;i++)for(let j=0;j<i;j++){
 let M=D[i].map((x,a)=>D[j].map((y,b)=>match(i,a,j,b)));if(M.every(r=>r.every(Boolean)))continue;
 G[i][j]=M;G[j][i]=D[j].map((y,b)=>M.map(r=>r[b]));
 Adj[i].push(j);Adj[j].push(i);
}
const bas=D.map(t=>t.flatMap((q,i)=>q.sg?[i]:[]));
function ac(I,j){
const R=bas.slice();R[I]=[j];let need=Array.from({length:N},(_,i)=>i),flag=new Array(N).fill(true),at=0,removed=0;
while(at<need.length){let k=need[at++];flag[k]=false;
 for(let g of Adj[k]){let a=R[g],r=G[g][k];let b=a.filter(x=>R[k].some(y=>r[x][y]));if(!b.length)return {bad:true,removed};if(b.length!=a.length){removed+=a.length-b.length;R[g]=b;if(!flag[g]){need.push(g);flag[g]=true}}}
}
return {bad:false,removed,R};
}

function ruleAC(R){
 let A=R.slice(), todo=A.map((_,i)=>i),forQ=A.map(()=>true),i=0;
 while(i<todo.length){let j=todo[i++];forQ[j]=false;for(let k of Adj[j]){
  let z=A[k], L=G[k][j], u=z.filter(x=>A[j].some(y=>L[x][y]));if(!u.length)return null;
  if(u.length!==z.length){A[k]=u;if(!forQ[k]){forQ[k]=true;todo.push(k)}}
 }}
 return A;
}

const CERT=[[4,2,0],[5,2,0],[5,5,0],[5,8,0],[5,11,0],[6,2,[5,[[3,0],[4,0]]]],[6,5,[5,[[0,0],[1,0]]]],[6,8,0],[6,11,0],[7,2,0],[9,2,0],[11,2,0],[12,2,[5,[[3,0],[4,0]]]],[12,5,[5,[[0,0],[1,0]]]],[12,8,0],[12,11,0],[13,2,0],[13,5,0],[13,8,0],[13,11,0],[14,2,0],[14,5,0],[15,2,0],[15,5,0],[16,2,0],[16,5,0],[16,8,0],[16,11,0],[17,2,[5,[[3,0],[4,0]]]],[17,5,[5,[[0,0],[1,0]]]],[17,8,0],[17,11,0],[18,2,0],[18,5,0],[19,2,0],[19,5,0],[20,2,0],[21,2,0],[22,2,0],[23,2,[5,[[3,0],[4,0]]]],[23,5,[5,[[0,0],[1,0]]]],[23,8,0],[23,11,0],[24,2,0],[24,5,0],[24,8,0],[24,11,0],[25,2,0],[27,2,0],[34,2,0],[34,5,0],[35,2,0],[36,2,0],[37,2,0],[45,2,0],[46,2,0],[53,2,0],[53,5,0],[54,2,0],[54,5,0],[55,2,0],[56,2,0],[65,2,0],[66,2,0],[67,2,0],[68,2,0],[68,5,0],[69,2,0],[71,2,0],[71,5,0],[72,2,[84,[[3,0],[4,0]]]],[72,5,[84,[[0,0],[1,0]]]],[72,8,0],[72,11,0],[73,2,[84,[[3,0],[4,0]]]],[73,5,[84,[[0,0],[1,0]]]],[73,8,0],[73,11,0],[74,2,0],[74,5,0],[75,2,0],[76,2,0],[77,2,0],[77,5,0],[78,2,[79,[[3,0],[4,0]]]],[78,5,[79,[[0,0],[1,0]]]],[78,8,0],[78,11,0],[79,2,[78,[[3,0],[4,0]]]],[79,5,[78,[[0,0],[1,0]]]],[79,8,0],[79,11,0],[80,2,0],[80,5,0],[81,2,0],[82,2,0],[82,5,0],[83,2,0],[84,2,0],[84,5,0],[84,8,0],[84,11,0],[85,2,[84,[[3,0],[4,0]]]],[85,5,[84,[[0,0],[1,0]]]],[85,8,0],[85,11,0],[86,2,0],[86,5,0],[87,2,0],[88,2,0],[89,2,0],[89,5,0],[90,2,0],[90,5,0],[91,2,0],[91,5,0],[92,2,0],[92,5,0],[93,2,0],[93,5,0],[94,2,0],[94,5,0],[95,2,0],[95,5,0],[96,2,[84,[[3,0],[4,0]]]],[96,5,[84,[[0,0],[1,0]]]],[96,8,0],[96,11,0],[97,2,0],[97,5,0],[97,8,0],[97,11,0],[98,2,0],[98,5,0],[99,2,[78,[[3,0],[4,0]]]],[99,5,[78,[[0,0],[1,0]]]],[99,8,0],[99,11,0],[100,2,[78,[[3,0],[4,0]]]],[100,5,[78,[[0,0],[1,0]]]],[100,8,0],[100,11,0],[101,2,0],[101,5,0],[102,2,0],[102,5,0],[103,2,0],[103,5,0],[104,2,0],[105,2,0],[106,2,[78,[[3,0],[4,0]]]],[106,5,[78,[[0,0],[1,0]]]],[106,8,0],[106,11,0],[107,2,[78,[[3,0],[4,0]]]],[107,5,[78,[[0,0],[1,0]]]],[107,8,0],[107,11,0],[108,2,0],[108,5,0],[109,2,0],[109,5,0],[110,2,0],[111,2,[84,[[3,0],[4,0]]]],[111,5,[84,[[0,0],[1,0]]]],[111,8,0],[111,11,0],[112,2,[84,[[3,0],[4,0]]]],[112,5,[84,[[0,0],[1,0]]]],[112,8,0],[112,11,0],[113,2,0],[113,5,0],[114,2,0],[114,5,0],[115,2,0],[116,2,0],[116,5,0],[117,2,0],[118,2,0],[118,5,0],[119,2,[78,[[3,0],[4,0]]]],[119,5,[78,[[0,0],[1,0]]]],[119,8,0],[119,11,0],[120,2,[79,[[3,0],[4,0]]]],[120,5,[79,[[0,0],[1,0]]]],[120,8,0],[120,11,0],[121,2,0],[121,5,0],[122,2,0],[124,2,0],[125,2,0],[126,2,0],[127,2,0],[134,2,0],[134,5,0],[137,2,0],[138,2,0],[145,2,0],[146,2,0],[147,2,0],[148,2,0],[149,2,0],[150,2,0],[157,2,0],[157,5,0],[158,2,0],[158,5,0],[160,2,0],[161,2,0],[162,2,0],[169,2,0],[170,2,0],[170,5,0],[172,2,0],[173,2,0],[173,5,0],[174,2,0],[174,5,0],[175,2,0],[176,2,0],[177,2,0],[178,2,0],[178,5,0],[179,2,0],[179,5,0],[180,2,0],[180,5,0],[181,2,0],[181,5,0],[182,2,0],[183,2,0],[185,2,0],[186,2,0],[186,5,0],[187,2,0],[187,5,0],[188,2,0]];
const E=Object.entries(H).sort();
assert.deepEqual(E,[["aux|single",30],["phase",18],["phase|single",18],["single",28]]);
assert.equal(BODY.size,32549); assert.equal(DB.length,32356); assert.equal(SB.length,193); assert.equal(OTHER.length,2963); assert.equal(DS.length,DB.length+SB.length);
assert.equal(BC.size,122); assert.equal(edge,542785); assert.equal(Q.length,448); assert.equal(direct.length,94);
for(let x of DB){assert.ok(rootpal[x.k]);assert.ok(x.ts[0].every((v,i)=>v===-x.ts[1][i]));}
for(let r of OTHER){const nonzero=r.slice(0,24).some(Boolean);const M=Math.max(...rootpal.map(t=>Math.abs(bdot(r.slice(24),t,0,5))));assert.ok(M<=(nonzero?32:64));}
assert.equal(U.target.length,240);assert.equal(ALLCOUNT(),853);
function ALLCOUNT(){return d29domain.reduce((s,D)=>s+D.length,0)}
assert.deepEqual(CERT.map(x=>x.slice(0,2)),U.target);
let visits=0,roots=0;
function verify(R,T){visits++;const A=ruleAC(R);if(T===0){assert.ok(A===null);return;}assert.ok(A!==null);const [i,V]=T;assert.equal(A[i].length,2);assert.deepEqual(V.map(x=>x[0]),A[i]);for(let [v,sub] of V){const X=A.slice();X[i]=[v];verify(X,sub);}}
for(let [i,k,p] of CERT){const R=bas.slice();R[i]=[k];verify(R,p);if(p!==0)roots++;}
assert.equal(visits,312);assert.equal(roots,36);
console.log(JSON.stringify({body:DS.length,doubles:DB.length,singles:SB.length,aux:OTHER.length,edges:edge,components:BC.size,unary_locks:aux,old_singleton_pins:pin,one_change_candidates:Q.length,pure_phase_candidates:direct.length,independent_refutations:visits,split_roots:roots,states:ALLCOUNT(),targets:U.target.length,new_lower_bound:false}));


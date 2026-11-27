
import { readFileSync } from "node:fs";
import { createHash } from "node:crypto";
import { gunzipSync } from "node:zlib";
import { join } from "node:path";
const dir=process.argv[2]??".";
function readPinned(rel,h){
 const u=readFileSync(join(dir,rel));
 if(createHash("sha256").update(u).digest("hex")!==h)throw Error("bad SHA256 "+rel);
 return u;
}
const extraText=gunzipSync(Buffer.from(readPinned("verifications/improved/dim60-61-rm80-puncture/extra60.json.gz.b64","8420af6106dcc6a6a2bda865593813e14e36f76dc547bac2accfcce320d05ac2").toString(),"base64")).toString();
const transText=gunzipSync(Buffer.from(readPinned("verifications/improved/dim60-61-rm80-guard-transfer/transfer60.txt.gz.b64","7fb3b5ed5b309bc87f05df08985e45e0f9d14c240835f4404b6912b748dcf65e").toString(),"base64")).toString();
const cubicText=readPinned("verifications/improved/dim62-67-rm2-cubic-cosets/witness.txt","5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f").toString();

const allcub=cubicText.trim().split(/\s+/).map(Number);
function B(x){x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (Math.imul((x+(x>>>4))&0x0f0f0f0f,0x01010101)>>>24);}
const PC=(u,v)=>B(u)+B(v);
const mon=m=>{let a=0,b=0;for(let i=0;i<64;i++)if((i&m)===m){if(i<32)a|=1<<i;else b|=1<<(i-32);}return [a,b];};
const mon2=Array.from({length:64},(_,i)=>i).filter(i=>B(i)<=2).map(mon);
const mon3=Array.from({length:64},(_,i)=>i).filter(i=>B(i)===3).map(mon);
const cube=allcub.map(c=>mon3.reduce(([a,b],[x,y],i)=>(c>>i&1)?[a^x,b^y]:[a,b],[0,0]));
const coord=(d,p)=>(d[p>>5]>>>(p&31))&1;
function project(v,points){return points.reduce((a,p,i)=>a|(coord(v,p)<<i),0);}
function mask(points){return points.reduce(([a,b],p)=>p<32?[a|(1<<p),b]:[a,b|(1<<(p-32))],[0,0]);}
const xor=(d,e)=>[d[0]^e[0],d[1]^e[1]];
const uint=n=>n>>>0;
function zerorank(rows,sites){
 const piv=Array(sites.length).fill(null), ker=[];
 function dec(x,e){
  for(let i=0;i<sites.length;i++)if(e>>i&1){if(!piv[i])return [x,e,i];x=xor(x,piv[i][0]);e^=piv[i][1];}
  return [x,e,sites.length];
 }
 for(const r of rows){const [x,e,i]=dec(r,project(r,sites));if(e)piv[i]=[x,e];else ker.push(x);}
 return {ker,piv,project:(r)=>dec(r,project(r,sites))};
}
function census(F=[],P=[]){
 const {ker,piv,project:fix}=zerorank(mon2,F);
 if(piv.some(x=>!x))throw "F not full";
 const T=1<<ker.length,lo=Int32Array.from(ker,v=>v[0]),hi=Int32Array.from(ker,v=>v[1]);
 const p=mask(P),res=[],cnt=Array(P.length+1).fill(0);
 for(let q=0;q<cube.length;q++){
  const [r,e]=fix(cube[q]);if(e)throw "can't fix F";
  let a=r[0],b=r[1];
  for(let i=0;i<T;i++){
    if(i){const k=31-Math.clz32(i&-i);a^=lo[k];b^=hi[k];}
    if(PC(a,b)===16){const c=PC(a&p[0],b&p[1]);cnt[c]++;res.push([a,b]);}
  }
 }
 return {cnt,res,rank:piv.length};
}
const hull=u=>u.reduce((p,m)=>[...p,...p.map(x=>x^m)],[0]);
const site=hull([3,20,41]);
const F=[0,20,41,61],P=[3,23,42,62];
const assert=(x,msg)=>{if(!x)throw Error(msg);};
const w=(v)=>v.reduce((s,x)=>s+x*x,0), dot=(u,v)=>u.reduce((s,x,i)=>s+x*v[i],0);
const gate=(d,n,m)=>d<=0||4*d*d<=n*m;
const slots=[...F,...P],iOf=new Map(slots.map((x,i)=>[x,i]));
const pol=slots.map((_,i)=>Array.from({length:4},(_,k)=>k===(i&3)?(i<4?-1:1):0));
const used=new Set;
const D=[];
function entry(label,h,expectParity){
 const N=56+w(h);
 assert(N<=64&&B(label)%2===expectParity,"prefix norms/parity");
 assert(!used.has(label),"duplicate prefix");used.add(label);
 D.push({label,h,N});
}
for(let i=0;i<8;i++)for(let j=i+1;j<8;j++)if((i&3)!==(j&3))
 entry((1<<i)|(1<<j),pol[i].map((x,k)=>-2*(x+pol[j][k])),0);
assert(D.length===24,"D4 tags");
for(let i=0;i<8;i++)entry(1<<i,pol[i].map(x=>-x),1);
entry(0,pol[4],0);
for(const a of D)for(const b of D)
 assert(gate(24+2*B(a.label^b.label)+dot(a.h,b.h),a.N,b.N),"dense gate");
const codeGenerators=[...mon2,...cube];
for(const u of codeGenerators)for(const v of codeGenerators)
 assert(!(PC(u[0]&v[0],u[1]&v[1])%2),"parent self orthogonality");
const ranks=zerorank(mon2,slots);
assert(ranks.ker.length===15&&ranks.piv.filter(Boolean).length===7,"marked rank");
assert(mon2.every(x=>B(project(x,slots))%2===0),"marked even");
const even=cube.filter(x=>!(B(project(x,slots))%2)).length;
assert(even===64&&allcub.length===80&&new Set(allcub).size===80,"coset parity ledger");
const den=[0,0];
for(const d of D)den[B(d.label)%2]++;
const dense=(even*den[0]+(80-even)*den[1])*2**15;
assert(dense===56623104,"dense count");
function atransfer(v,pattern){
 return P.map((j,k)=>(((pattern>>k)&1)?dot(pol[iOf.get(j)],v.h):0));
}
for(const d of D){
 for(let r=0;r<16;r++){
  const c=B(r),a=atransfer(d,r);
  const z0=4-c+a.reduce((s,x)=>s+Math.abs(x),0);
  assert(gate(z0,d.N,4),"four-site to dense");
  let z1=c===3 ? (16-c+a.reduce((s,t)=>s+t,0)) : z0+12;
  if(!gate(z1,d.N,16)){
    for(let k=0;k<4;k++)if(r>>k&1){
      const sign=(d.label>>(k+4)&1)?-1:1;
      assert(c===3 ? sign===1 : sign*a[k]>=1,"parity mismatch coefficient");
    }
    z1-=2;
  }
  assert(gate(z1,d.N,16),"sparse to dense");
 }
 assert(gate(1,d.N,1)&&d.h.every(x=>gate(Math.abs(x),d.N,1)),"pole to dense");
}
assert(gate(8,16,16)&&gate(4,16,4)&&gate(2,4,4)&&gate(1,16,1)&&gate(1,4,1),"analytic sparse gates");
function hamming16(){
 const S=new Array(32).fill(0);
 let low=16;
 for(let j=0;j<1<<16;j++){
   let sig=0;
   for(let k=0;k<16;k++)if(j>>k&1)sig^=16+k;
   S[sig]++;if(j&&sig===0)low=Math.min(low,B(j));
 }
 assert(S[16]===2048&&low===4,"odd Hamming bank 2048");
}
hamming16();

function stock(V,extraText,transText){
 const nn=64n, key=(u,v)=>(BigInt(u>>>0))|(BigInt(v>>>0)<<32n),
   unpack=w=>[Number(w&0xffffffffn)|0,Number(w>>32n)|0];
 const E=JSON.parse(extraText.replace(/\d+/g,s=>JSON.stringify(s)));
 const T=transText.trim().split(/\s*(?:\r\n|\n)\s*/);
 const [pn,qn]=T[0].split(/\s+/).map(Number);
 assert(pn===336&&qn===880&&T.length===1+pn+qn,"transfer shape");
 const ws=[...E.extra,...T.slice(1,1+pn)].map(BigInt);
 assert(E.extra.length===672&&new Set(ws).size===1008&&E.guards.length===2232,"guarded shape");
 const SS=[...V.res,...ws.map(unpack)];const n=SS.length;
 const idx=new Map(SS.map(([a,b],i)=>[key(a,b),i]));
 assert(n===13013&&idx.size===n,"complete distinct stock");
 const mm=mask(P), ff=mask(F);
 const hist=Array(5).fill(0); const freq=(a,b)=>PC(a&mm[0],b&mm[1]);
 assert(SS.every(([a,b])=>PC(a,b)===16&&PC(a&ff[0],b&ff[1])===0),"transferred masks");
 assert(SS.every(([a,b])=>codeGenerators.every(([u,v])=>PC(a&u,b&v)%2===0)),"all sparse support parity");
 SS.forEach(([a,b])=>hist[freq(a,b)]++);
 const g0=E.guards,g1=T.slice(1+pn).map(row=>row.split(/\s+/));
 const guards=new Set,S=Array.from({length:n},(_,i)=>[[key(...SS[i]),1]]);
 for(const [da,db,dg,ct] of [...g0,...g1]){
   const a=BigInt(da),b=BigInt(db),g=BigInt(dg),t=Number(ct);
   const i=idx.get(a),j=idx.get(b);
   assert(i!==undefined&&j!==undefined&&i!==j&&(t===0||t===1),"literal guard identifiers");
   assert(g>0n&&((a&b)&g)===g&&PC(...unpack(g))===8,"literal joint check");
   const k=Math.min(i,j)*n+Math.max(i,j);assert(!guards.has(k),"guard uniqueness");guards.add(k);
   S[i].push([g,t^1]);S[j].push([g,t]);
 }
 assert(guards.size===3112,"guard population");
 let hits=0;
 const ss0=Int32Array.from(SS,x=>x[0]),ss1=Int32Array.from(SS,x=>x[1]);
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
   const c=PC(ss0[i]&ss0[j],ss1[i]&ss1[j]);
   if(c<=8)continue;
   assert(c<=10&&guards.has(i*n+j),"uncovered actual support overlap");
   hits++;
 }
 assert(hits===guards.size,"all guards are necessary links");
 for(let i=0;i<S.length;i++){
   const [a,b]=SS[i],full=key(a,b),H=new Map([[0n,0n]]),gg=[];
   for(const [r,t] of S[i]){
      if(H.has(r)){assert(H.get(r)===BigInt(t),"inconsistent signed checks");continue;}
      if(r!==full)gg.push([r,t]);
      const old=[...H];for(const [x,y] of old)H.set(x^r,y^BigInt(t));
   }
   assert(H.size===(2<<gg.length)&&gg.length<=4,"bounded guard rank");
   assert([...H.keys()].every(z=>z===0n||z===full||PC(...unpack(z))===8),"all homogeneous checks");
   const us=Array.from({length:1<<gg.length},()=>[]);
   for(let x=0;x<64;x++)if((full>>BigInt(x))&1n){
     let v=0;for(let j=0;j<gg.length;j++)v|=Number((gg[j][0]>>BigInt(x))&1n)<<j;
     us[v].push(x);
   }
   const b0=16/us.length;
   assert(us.every(a=>a.length===b0),"balanced evaluation columns");
   // Deterministically complete the affine Hamming dual: low bits the
   // old independent check values, high bits index within each chart.
   const fullCols=new Array(64).fill(-1);
   us.forEach((qs,j)=>qs.forEach((x,k)=>{fullCols[x]=j+(k<<gg.length);}));
   const pm=Array.from({length:64},(_,x)=>fullCols[x]).filter(x=>x>=0);
   assert(new Set(pm).size===16&&pm.every(x=>x>=0&&x<16),"explicit Hamming completion");
   // Three marked + signs specify independent linear conditions: no
   // dual word supported there. This is also checked on the actual
   // three Hamming columns, augmented by the overall parity row.
   if(freq(a,b)===3){
     const mt=P.filter(x=>fullCols[x]>=0).map(x=>16+fullCols[x]);
     assert(mt.length===3&&new Set(mt).size===3,"actual marked triple");
     let rank=0;
     for(let h=0;h<32;h++){
       const v=mt.reduce((r,z,j)=>r|((B(h&z)%2)<<j),0);
       if(!v)rank++;
     }
     assert(rank===4,"rank-three marked restriction of row-space pairing");
     // Equivalently, all three fixed real signs cost 3 independent
     // coordinate constraints on the primal sign code, since the dual
     // check space has minimum nonzero weight 8. The primal loses 3 bits.
   }
 }
 return {hist,guards:hits,orig:V.cnt};
}

const V=census(F,P);
assert(V.cnt.join()==="5261,1924,4456,244,120","parent weight 16 census");
const ext=stock(V,extraText,transText), hist=ext.hist;
assert(hist.join()==="5315,2735,4519,324,120","final sparse bins");
const planes=new Set;
for(let i=0;i<64;i++)for(let j=i+1;j<64;j++)for(let k=j+1;k<64;k++){
 const l=i^j^k;if(l>k&&[i,j,k,l].every(x=>!F.includes(x)))planes.add([i,j,k,l].join());
}
assert(planes.size===7995,"affine 4-site layer");
const sparse=hist.reduce((s,z,i)=>s+z*(i===3?256:2048),0),quad=planes.size*16,poles=120,total=dense+sparse+quad+poles;
assert(total===82821160,"total");
console.log({status:"PASS",dimension:60,evenCosets:even,prefixCount:den,dense,sparse,guards:ext.guards,sparseBins:hist,quad,poles,total,previous:69704624,gain:total-69704624});

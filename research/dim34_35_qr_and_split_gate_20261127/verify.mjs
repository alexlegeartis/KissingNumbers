import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {resolve} from 'node:path';
import {inflateRawSync} from 'node:zlib';
const K=(( )=>{
const K0=['4000152727','2000C0DE0B','1000AA229D','8009F5CD6','4004FAE6B','200ED9AAD','100BC80CE','805E4067','40E56DAB','20B8FB4D','1096303E','84B181F','4EFC197','2BDAD53','1949B31'].map(s=>BigInt("0x"+s));
const R=['0','ADCFB','F7795F','740E29','940DFF','3B66AB','6CD301','4709D4','98A926','D23688'].map(s=>BigInt("0x"+s));
const W=[0n];for(const a of K0)W.push(...W.map(w=>w^a));
if(new Set(W).size!==32768)throw Error('rank');
const p32=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
const wpop=x=>p32(Number(x&0xffffffffn))+p32(Number(x>>32n));
const bits=x=>{const b=[];for(let i=0;i<39;i++)if(x>>BigInt(i)&1n)b.push(i);return b;};
const Bad4=new Uint8Array(39**4);
function idx(a,b,c,d){return ((a*39+b)*39+c)*39+d;}
for(const w of W)if(wpop(w)===12){const B=bits(w);for(let i=0;i<12;i++)for(let j=i+1;j<12;j++)for(let k=j+1;k<12;k++)for(let l=k+1;l<12;l++)Bad4[idx(B[i],B[j],B[k],B[l])]=1;}

const safe=[];let five=0;
for(let a=0;a<39;a++)for(let b=a+1;b<39;b++)for(let c=b+1;c<39;c++)for(let d=c+1;d<39;d++){
 if(Bad4[idx(a,b,c,d)])continue;safe.push([a,b,c,d]);
 for(let e=d+1;e<39;e++)if(!Bad4[idx(a,b,c,e)]&&!Bad4[idx(a,b,d,e)]&&!Bad4[idx(a,c,d,e)]&&!Bad4[idx(b,c,d,e)])five++;
}
return {safe,five};
})();
assert.equal(K.safe.length,23);assert.equal(K.five,0);
const Q=((M)=>{
const G0=['4000152727','2000C0DE0B','1000AA229D','8009F5CD6','4004FAE6B','200ED9AAD','100BC80CE','805E4067','40E56DAB','20B8FB4D','1096303E','84B181F','4EFC197','2BDAD53','1949B31'].map(s=>BigInt("0x"+s));
const P=Array(39).fill(0n);
for(let x of G0){while(x){let j=x.toString(2).length-1;if(!P[j]){P[j]=x;break;}x^=P[j];}if(!x)throw Error('rank');}
const F=Array.from({length:39},(_,i)=>i).filter(i=>!P[i]);
function col(i){
 let x=1n<<BigInt(i);
 for(let j=38;j>=0;j--)if(P[j]&&(x>>BigInt(j)&1n))x^=P[j];
 return F.reduce((m,j,k)=>m|((x>>BigInt(j)&1n)?1<<k:0),0);
}
const cols=Array.from({length:39},(_,i)=>col(i));if(F.length!==24)throw Error('quotient rank');
const N=1<<24,D=new Uint8Array(N);D.fill(80);D[0]=0;
for(let x of cols){
  const b=x?2**(31-Math.clz32(x)):0;if(!b)continue;
  for(let a=0;a<N;a+=2*b)for(let j=a;j<a+b;j++){
   const k=j^x,u=D[j],v=D[k]; D[j]=Math.min(u,v+1);D[k]=Math.min(v,u+1);
  }
}
let hist=Array(40).fill(0),far7=[],far8=[];
for(let i=0;i<N;i++){
 if(D[i]>=hist.length)throw Error('not expanded');
 hist[D[i]]++;
 if(D[i]>=8)far7.push(i);if(D[i]>=9)far8.push(i);
}
const post=M.map(m=>{
 let sh=[0];for(let j of m)sh.push(...sh.map(s=>s^cols[j]));
 if(new Set(sh).size!==16)throw Error('mark invariant');
 let f8=0,f9=0,b=[];
 for(let i of far7)if(sh.every(s=>D[i^s]>=8))f8++;
 for(let i of far8)if(sh.every(s=>D[i^s]>=9)){f9++;if(b.length<10)b.push(i);}
 return {m,n8:f8/16,n9:f9/16,b};
});
return {hist,post,cols};
})(K.safe);
assert.deepEqual(Q.hist.filter(x=>x),[1,39,741,9139,82251,575757,2836720,6944617,5468341,859047,554,9]);
assert(Q.post.every(x=>x.n9===0&&x.n8>0));
console.log('QR: 23 safe four-punctures, all covering radius eight, and no safe five-puncture.');
const input=(p,h)=>{const b=readFileSync(resolve(process.argv[2]??'.','verifications/improved',p));assert.equal(createHash('sha256').update(b).digest('hex'),h);return b;};
const decode=(src)=>{const bytes=new Uint8Array(src);
const dv=new DataView(bytes.buffer), ds=(x,n)=>String.fromCharCode(...bytes.slice(x,x+n));
const files={};
for(let p=0;p+30<bytes.length&&dv.getUint32(p,true)===0x04034b50;){
 const t=dv.getUint16(p+8,true),fn=dv.getUint16(p+26,true),ex=dv.getUint16(p+28,true);
 const name=ds(p+30,fn).replace(/\.npy$/,'');
 let sz=dv.getUint32(p+22,true),cs=dv.getUint32(p+18,true);
 if(sz===0xffffffff){const q=p+30+fn;if(dv.getUint16(q,true)!==1)throw Error('zip64');sz=Number(dv.getBigUint64(q+4,true));cs=Number(dv.getBigUint64(q+12,true));}
 p+=30+fn+ex;const S=bytes.slice(p,p+cs);p+=cs;
 const C=t===0?S:t===8?new Uint8Array(inflateRawSync(S)):null;if(C===null||C.length!==sz)throw Error('zip method');
 const v=new DataView(C.buffer),len=v.getUint16(8,true),h=String.fromCharCode(...C.slice(10,10+len)),match=h.match(/'descr': '(.+?)'.*'shape': \(([^)]+)\)/s);
 if(!match||!h.includes("'fortran_order': False"))throw Error('NPY parser');
 const shape=match[2].split(',').filter(x=>x.trim()).map(Number),k=Number(match[1].slice(2)),off=10+len;
 const getter=k===8?j=>Number(v.getBigInt64(j,true)):k===4?j=>v.getInt32(j,true):k===2?j=>v.getInt16(j,true):k===1?j=>v.getInt8(j):()=>{throw Error('np dtype')};
 const A=Array.from({length:shape.reduce((a,b)=>a*b,1)},(_,i)=>getter(off+k*i));
 files[name]=shape.length===1?A:Array.from({length:shape[0]},(_,i)=>A.slice(i*shape[1],(i+1)*shape[1]));
}
return files;
};
const code=JSON.parse(input('dim35-shared-root-heads/code_data.json','0f44ca2ff564941d58f92cbb0c47ada139612196fd3d69a547d5318f78009d05'));
const wit=JSON.parse(input('dim35-octad-switches/witness.json','3e8a212b8a8ee4c6d81a45134dae67b99fbfbc4a4f990fc174c42af987e01e16'));
const positions=decode(input('dim35-shared-root-heads/shared_positive.npz','9d4b35ac5d18f40011282ec2e83753bdffee947351584a77a910baef875768cf')).dense_positions;
const data=decode(input('dim35-joint-star-root/replacement.npz','ffb73e3082a0ed5a1f94ce8e5c32febb0683e151bcd1aee3838eb52a6e76320c'));
assert.deepEqual([data.stars.length,data.roots.length,data.heads.length],[262,2158,12]);
const D=((code,wit,S,positions)=>{
const row=Object.values(S).flat();
const D35=35,LIM=Math.sqrt(256/3),B=1.05,HP=Math.sqrt(16/3-B*B);
const dp=positions.slice(),D=dp.map((x,i)=>[x,i]),sp=Array.from({length:35},(_,i)=>i).filter(x=>!dp.includes(x));
const id=Array(35).fill(-1);D.forEach(([x,i])=>id[x]=i);
const mask=L=>L.reduce((a,j)=>a+(2**j),0);
const B32=L=>L.reduce((a,j)=>a|(1<<j),0)>>>0;
const gs=code.dense_generator,shift=code.dense_coset_shift;
const wt=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
let A=code.full_parity_supports.map(([m])=>m),norm=row.map(g=>g.reduce((s,x)=>s+x*x,0));
for(const [i,old,nex] of wit.octads){if(mask(old)!==A[i])throw Error('source octad');A[i]=mask(nex);}
const O=A.map((m,index)=>{
const o=Array.from({length:35},(_,i)=>i).filter(j=>Math.floor(m/2**j)%2);
 if(o.length!==8)throw Error('octad');
 const t=sp.filter(j=>o.includes(j)),tail=B32(t.map(i=>sp.indexOf(i))),body=B32(o.filter(j=>id[j]>=0).map(j=>id[j]));
 let allow=3;
 for(let k=0;k<row.length;k++){
  const g=row[k],N=norm[k];let s=0,mi=6,h=0;
  for(const j of o){const v=g[j],a=Math.abs(v);s+=a;if(a<mi)mi=a;h^=+(v<0);}
  if(4*(2*s)**2>32*N)allow&=~(1<<h);
  if(4*(2*s-4*mi)**2>32*N)allow&=~(1<<(h^1));
  if(!allow)throw Error('prior phase');
 }
 const eq=gs.reduce((a,w,i)=>a|((wt(w&body)&1)<<i),0);
 return {index,m,o,t,tail,body,allow,eq,rhs:wt(shift&body)&1};
});
function eqbasis(E){
 const U=Array(17),v=new Uint8Array(17);let r=0;
 for(let [a,b] of E){
  while(a){const j=31-Math.clz32(a);
   if(U[j]){a^=U[j];b^=v[j];}
   else{U[j]=a;v[j]=b;r++;break;}
  }
  if(!a&&b)return null;
 }
 let x=0,free=[];
 for(let j=0;j<17;j++)if(U[j]){
  if((wt(U[j]&x)&1)!==v[j])x|=1<<j;
 }else free.push(j);
 return {U,v,x,r,free};
}
function proj(k,U){
 let x=0;for(let j=0;j<17;j++)if((k>>j)&1)x^=U[j];
 return x>>>0;
}
const fixedCount=O.filter(x=>x.allow===1||x.allow===2).length;
return {O,sp,B,HP,dp,row,norm,gs,shift,wt,eqbasis,proj,fixedCount,LIM};
})(code,wit,data,positions);

const {O,sp,dp,eqbasis,proj,gs}=D;
function geo(o,t0,t1,r,i,j){
 let b0=o.o.includes(t0),b1=o.o.includes(t1),r0=o.o.includes(r),a=+b0+ +b1;
 let c=+o.o.includes(i),d=+o.o.includes(j),q=c+d;
 if(a===2&&(!r0&&q===2||r0&&q===1))return -2;
 if(r0||a===0||a===1&&q!==1)return -1;
 if(!q||c)return 0;
 return a===1?1:0;
}

let tested=0,geoCount=0,inconsistent=0;
for(let a=0;a<3;a++)for(let c=a+1;c<3;c++){
 const p=sp[a],q=sp[c],r=sp[3-a-c];
 for(let i=0;i<32;i++)for(let j=i+1;j<32;j++){
  tested++;const Q=O.map(o=>geo(o,p,q,r,dp[i],dp[j]));
  if(Q.includes(-2))continue;
  geoCount++;const V=O.filter((_,k)=>Q[k]>=0), eq0=V.map(o=>[o.eq,0]);
  assert.equal(eqbasis(eq0).r,17);
  for(let m=0;m<4;m++){
   const phase=o=>((m&1)*o.o.includes(p))^((m>>1)*o.o.includes(q));
   const E=eqbasis(V.map(o=>[o.eq,phase(o)]));
   if(m===0||m===3){assert(E);assert.equal(proj(E.x,gs)>>>0,m?0xffffffff:0);}
   else assert.equal(E,null);
   const C=V.filter(o=>o.allow!==3).map(o=>[o.eq,o.rhs^Q[o.index]^phase(o)^(o.allow===1?1:0)]);
   if(eqbasis(C)===null)inconsistent++;
  }
 }
}
assert.deepEqual([tested,geoCount,inconsistent,D.fixedCount],[1488,19,76,201]);
console.log('35D split: 1488 geometries, 19 geometrically legal, all 76 old-parity systems inconsistent; every shape rank 17.');

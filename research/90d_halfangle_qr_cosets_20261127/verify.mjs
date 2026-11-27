// QR90 half-angle barriers: Node >=20, no package dependencies.
import assert from 'node:assert/strict';
const q=x=>(x%3+3)%3, sum=x=>x.reduce((a,b)=>a+b,0);
const pc=a=>{a=a-((a>>>1)&0x55555555);a=(a&0x33333333)+((a>>>2)&0x33333333);return Math.imul((a+(a>>>4))&0x0f0f0f0f,0x01010101)>>>24;};
const poly=a=>{const z=Array(24).fill(0); for(const [i,j] of a)z[i]=q(j);return z;};
const g=poly([[23,1],[20,1],[18,-1],[17,1],[15,-1],[12,1],[10,1],[9,-1],[7,-1],[6,-1],[5,-1],[2,-1],[1,-1],[0,-1]]);
const h=poly([[23,1],[22,1],[21,1],[18,1],[17,1],[16,1],[14,1],[13,-1],[11,-1],[8,1],[6,-1],[5,1],[3,-1],[0,-1]]);
const mul=(a,b)=>{const c=Array(a.length+b.length-1).fill(0);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]=q(c[i+j]+a[i]*b[j]);return c;};
const f=mul([2,1],g),z=mul(f,h);
assert.deepEqual(z,[2,...Array(46).fill(0),1]);
assert(q(sum(g))&&q(sum(h)));
let ord=1,p=3;while(p!==1){++ord;p=3*p%47;} assert.equal(ord,23);
// g generates one of the two quadratic-residue cyclic codes via (X-1)g and the all-one vector.
const E=Array.from({length:24},(_,i)=>i===23?Array(48).fill(1):Array.from({length:48},(_,j)=>f[j-i]??0));
for(let i=0;i<24;i++){
 let j=i;while(j<24&&!E[j][i])++j;assert(j<24);
 [E[i],E[j]]=[E[j],E[i]];
 if(E[i][i]===2)E[i]=E[i].map(a=>q(2*a));
 for(let t=0;t<24;t++)if(t!==i&&E[t][i]){const s=E[t][i];E[t]=E[t].map((a,k)=>q(a-s*E[i][k]));}
}
const R=E.map((row,i)=>{for(let j=0;j<24;j++)assert.equal(row[j],Number(j===i));return row.slice(24);});
const IR=R.map((r,i)=>R.map(row=>q(-row[i])));
for(let i=0;i<24;i++)for(let j=0;j<24;j++)assert.equal(q(sum(R[i].map((a,k)=>a*IR[k][j]))),Number(i===j));
const bits=(R,s)=>R.map(row=>sum(row.map((a,k)=>a===s?1<<k:0)));
const G=Array.from({length:27},()=>[]);
function go(xa,xb,ya,yb,i,w,A,B,rev){
 if(w+pc(ya|yb)===15){
  let a=rev?ya:xa,b=rev?yb:xb,c=rev?xa:ya,d=rev?xb:yb,j=0;
  for(let n=2;n>=0;n--)j=3*j+((a>>>n)&1)+2*((b>>>n)&1);
  G[j].push([a&~7,b&~7,c,d]);
 }
 if(w===7)return;const nn=~(ya|yb);
 for(let j=i;j<24;j++){
  const bit=1<<j,a=A[j],b=B[j],mm=~(a|b);
  go(xa|bit,xb,(ya&mm)|(nn&a)|(yb&b),(yb&mm)|(nn&b)|(ya&a),j+1,w+1,A,B,rev);
  go(xa,xb|bit,(ya&mm)|(nn&b)|(yb&a),(yb&mm)|(nn&a)|(ya&b),j+1,w+1,A,B,rev);
 }
}
go(0,0,0,0,0,0,bits(R,1),bits(R,2),false);
go(0,0,0,0,0,0,bits(IR,1),bits(IR,2),true);
assert.equal(sum(G.map(x=>x.length)),415104);
assert.equal(G[0].length,130944);
const maxLog={0:[0,1],1:[0,1],2:[0,1],3:[0,1],4:[0,1]};
function maybeMax(i,num,n){if(num*maxLog[i][1]>n*maxLog[i][0])maxLog[i]=[num,n];}
for(let m=1;m<27;m++){
 const markWeight=[m%3,Math.floor(m/3)%3,Math.floor(m/9)].filter(Boolean).length;
 const w=15-markWeight,L=G[m],n=L.length;
 assert.equal(n,[0,31680,6930,1365][markWeight]);
 // Positive bad offsets at zero relative rotation already affinely span C.
 // In systematic coordinates, deleted marks are fixed and the last 21 left trits identify all C differences.
 const left=([a,b])=>Array.from({length:21},(_,i)=>((a>>>(i+3))&1)+2*((b>>>(i+3))&1));
 const ref=left(L[0]),B=Array(21).fill(null);let rank=0;
 for(let t=1;t<L.length&&rank<21;t++){
  const y=left(L[t]).map((a,i)=>q(a-ref[i]));
  for(let i=0;i<21;i++)if(y[i]){
   if(B[i]){const c=y[i];for(let j=i;j<21;j++)y[j]=q(y[j]-c*B[i][j]);}
   else{const c=y[i];for(let j=i;j<21;j++)y[j]=q(c*y[j]);B[i]=y;++rank;break;}
  }
 }
 assert.equal(rank,21);
 const ct=new Int32Array(45),sp=new Int32Array(45);
 const cov=new Int32Array(45*45), cub=new Int32Array(45**3);
 const N1=new Int32Array(90), P1=new Int32Array(90*45), M1=new Int32Array(90*45);
 const N2=new Int32Array(90*90),P2=new Int32Array(90*90*45),M2=new Int32Array(90*90*45);
 for(const [a,b,c,d] of L){
  const I=[],S=[],K=[];
  for(let i=3;i<48;i++){
   const p=(i<24?a>>>i:c>>>(i-24))&1,r=(i<24?b>>>i:d>>>(i-24))&1;
   assert(!(p&&r));
   if(p||r){I.push(i-3); S.push(p-r); K.push(2*(i-3)+r); ++ct[i-3]; sp[i-3]+=p-r;}
  }
  assert.equal(I.length,w);
  for(let i=0;i<w;i++){
   ++N1[K[i]];
   for(let t=0;t<w;t++)if(t!==i)++(S[t]===1?P1:M1)[45*K[i]+I[t]];
   for(let j=i+1;j<w;j++){
    if(markWeight===1){
     const p=I[i]*45+I[j],x=S[i]*S[j];
     cov[p]+=x;for(let k=j+1;k<w;k++)cub[p*45+I[k]]+=x*S[k];
    }
    const u=90*K[i]+K[j]; ++N2[u];
    for(let t=0;t<w;t++)if(t!==i&&t!==j)++(S[t]===1?P2:M2)[45*u+I[t]];
   }
  }
 }
 if(markWeight===1){
  assert(sp.every(a=>a===0));
  const abs=xs=>xs.reduce((a,b)=>a+Math.abs(b),0);
  const v=w*n-2*abs(cov),t=6*abs(cub);
  assert(v>2*n&&Number.isSafeInteger(v*v)&&v*v-4*v*n-2*t*n>0);
 }else{
  const dbl=w*n+sp.reduce((a,b)=>a+Math.abs(b),0);
  // wrong-sign count <=6 and weighted loss <=19 (indeed <=18 for three marks).
  assert(dbl<14*n);
  maybeMax(0,2*w*n+dbl,2*n);
 }
 for(let u=0;u<90;u++){
  const n=N1[u];assert(n>0);
  let c=0;
  for(let i=0;i<45;i++){const p=P1[45*u+i],m=M1[45*u+i];c+=p+m+Math.max(p,m);}
  assert(c<20*n); maybeMax(1,c,n);
 }
 let pairs=0;
 for(let u=0;u<90;u++)for(let v=(u>>1)*2+2;v<90;v++){
  ++pairs;const a=90*u+v,n=N2[a];assert(n>0);
  let c=0,d1=0,d2=0;
  for(let i=0;i<45;i++){
   const p=P2[45*a+i],m=M2[45*a+i],z=n-p-m;
   const t=p+m+Math.max(p,m);c+=t;
   if(i===(u>>1)||i===(v>>1))continue;
   const d=Math.max(p+z+Math.max(p,z),m+z+Math.max(m,z))-t;
   if(d>=d1){d2=d1;d1=d;}else if(d>d2)d2=d;
  }
  assert(c<20*n && c+d1<20*n);
  maybeMax(2,c,n);maybeMax(3,c+d1,n);
  if(markWeight===2){assert(c+d1+d2<20*n);maybeMax(4,c+d1+d2,n);}
 }
 assert.equal(pairs,3960);
}
assert(3*26**2>45**2);
console.log('AUDIT PASS QR15=415104; 26 nonzero mark classes; all signs and rotations at 30 degrees; all defects through three; two-mark class through four; unrotated positive bad offsets full rank',maxLog);

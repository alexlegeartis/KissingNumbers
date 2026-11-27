// Save as verify.mjs; Node 20+: node verify.mjs certificate.json
import {readFileSync} from "node:fs";
import {createHash} from "node:crypto";
const original=readFileSync(process.argv[2]);
if(createHash("sha256").update(original).digest("hex")!=="c4598a0dfc49365ae2b082745454c38b53b76703e5eef5dafbd71e2d06c59180")throw Error("current skeleton source identity");
const C=JSON.parse(original);
function loadCT(C) {
 const G=C.gram_12, M=C.sublattice_matrix, n=12, range=n=>Array.from({length:n},(_,i)=>i);
 const need=(x,label)=>{if(!x)throw Error(label)};
 const mod=x=>(x%3+3)%3, pw=range(6).map(i=>3**i);
 const E=range(n).map(i=>range(n).map(j=>+(i===j)-M[i][j]));
 const product=(A,B)=>A.map(r=>range(n).map(j=>r.reduce((v,x,i)=>v+x*B[i][j],0)));
 const transpose=A=>range(n).map(i=>A.map(r=>r[i]));
 const eq=(A,B)=>A.every((r,i)=>r.every((x,j)=>x===B[i][j]));
 need(eq(product(product(M,G),transpose(M)),G.map(r=>r.map(x=>3*x))),"axial similarity");
 need(eq(product(product(E,G),transpose(E)),G),"Eisenstein rotation isometry");
 need(eq(product(E,E),E.map((r,i)=>r.map((x,j)=>-x-+(i===j)))),"Eisenstein order-three identity");

 const gc=(a,b)=>b?gc(b,a%b):Math.abs(a),lcm=(a,b)=>a/gc(a,b)*b;
 function F(a,b=1){if(b<0){a=-a;b=-b;}const g=gc(a,b);need(Number.isSafeInteger(a)&&Number.isSafeInteger(b)&&b>0,"rational exact");return [a/g,b/g];}
 const add=([a,b],[c,d])=>F(a*d+b*c,b*d),sub=([a,b],[c,d])=>F(a*d-b*c,b*d),mul=([a,b],[c,d])=>F(a*c,b*d),div=([a,b],[c,d])=>F(a*d,b*c);
 let L=range(n).map(i=>range(n).map(j=>F(+(i===j)))),D=[];
 for(let j=0;j<n;j++){
  let z=F(G[j][j]);for(let k=0;k<j;k++)z=sub(z,mul(D[k],mul(L[j][k],L[j][k])));
  need(z[0]>0,"PD");D.push(z);
  for(let i=j+1;i<n;i++){let t=F(G[i][j]);for(let k=0;k<j;k++)t=sub(t,mul(D[k],mul(L[j][k],L[i][k])));L[i][j]=div(t,z);}
 }
 const denL=range(n).map(j=>L.reduce((p,x)=>lcm(p,x[j][1]),1));
 const pp=range(n).map(j=>D[j][1]*denL[j]**2), J=pp.reduce(lcm,1);
 const Q=range(n).map(j=>J/pp[j]*D[j][0]), K=range(n).map(j=>range(n).map(i=>L[i][j][0]*(denL[j]/L[i][j][1])));
 need(Number.isSafeInteger(8*J),"global integer bound");
 for(let a=0;a<n;a++)for(let b=0;b<n;b++)need(K.reduce((s,k,i)=>s+Q[i]*k[a]*k[b],0)===J*G[a][b],"LDL equals G");
 const GP=v=>G.map(r=>r.reduce((s,g,j)=>s+g*v[j],0));
 const ip=(x,gy)=>x.reduce((s,v,j)=>s+v*gy[j],0);
 const label=v=>pw.reduce((s,p,i)=>s+p*mod(v[C.label_first[i]]+v[C.label_second[i]]),0);
 const R=v=>v.map((t,j)=>t-v.reduce((s,x,i)=>s+x*M[i][j],0));
 const neg=k=>pw.reduce((s,p)=>s+p*mod(-(Math.floor(k/p)%3)),0);
 const diff=(a,b)=>pw.reduce((s,p)=>s+p*mod(Math.floor(a/p)%3-Math.floor(b/p)%3),0);
 const v=range(n).fill(0), S=range(729).map(()=>[]),sh=range(9).map(()=>[]);
 const chosen=[];
 function isqrt(x){let y=Math.floor(Math.sqrt(x));while((y+1)**2<=x)y++;while(y*y>x)y--;return y;}
 function walk(j,rem){
  if(j<0){const u=v.slice(), N=ip(u,GP(u)), ell=label(u);need(J*N===8*J-rem,"valid leaf");sh[N].push(u);S[ell].push([N,u]);return;}
  const s=K[j].reduce((a,k,i)=>a+(i>j?k*v[i]:0),0), p=isqrt(Math.floor(rem/Q[j])), a=K[j][j];
  for(let t=Math.ceil((-p-s)/a);t<=Math.floor((p-s)/a);t++) {v[j]=t;walk(j-1,rem-Q[j]*(a*t+s)**2);}
 }
 walk(11,8*J);
 need(JSON.stringify(sh.map(a=>a.length))===JSON.stringify([1,0,0,0,756,0,4032,0,20412]),"K12 short shells");
 const mn=S.map(x=>Math.min(...x.map(u=>u[0]))), banks=S.map((x,i)=>x.filter(u=>u[0]===mn[i]).map(u=>u[1]));
 const reps=banks.map((b,i)=>{if(mn[i]!==4)return null;let a=b[0];return [a,R(a),R(R(a))];});
 for(let i=0;i<729;i++)if(mn[i]===4)need(new Set(reps[i].map(v=>v.join(","))).size===3&&reps[i].every(v=>label(v)===i&&ip(v,GP(v))===4),"well-formed phase numbering");
 return {G,M,ip,GP,label,R,neg,diff,sh,bank:banks,mn,reps,J,D,L,Q,K};
}
function verifyRadial(C, X) {
 const {R,GP,ip,sh,bank,mn,neg,label}=X, range=n=>Array.from({length:n},(_,i)=>i), need=(x,m)=>{if(!x)throw Error(m);};
 const labs=range(729).filter(c=>mn[c]===6), B=labs.length, V=labs.flatMap(c=>bank[c]);
 need(B===224&&V.length===4032,"norm-six labels/rows");
 const by=Int16Array.from(range(729),()=>-1);labs.forEach((c,i)=>by[c]=i);
 const Tr=labs.map(c=>{
  const seen=new Set(),out=[],inBank=new Set(bank[c].map(t=>t.join(",")));
  for(const u of bank[c])if(!seen.has(u.join(","))){
   const r=[u,R(u),R(R(u))];need(r.every(x=>inBank.has(x.join(","))),"phase closure");
   r.forEach(x=>seen.add(x.join(",")));
   for(let i=0;i<3;i++)for(let j=0;j<3;j++)need(ip(r[i],GP(r[j]))===(i===j?6:-3),"trine");
   need(!out.some(t=>r.some(x=>t.some(y=>ip(x,GP(y))!==0))),"orthogonal trines");out.push(r);
  }need(out.length===6,"six trines");return out;
 });
 const maps=range(3).map(()=>range(3).map(()=>new Int16Array(B).fill(-1)));
 for(const c of C.outer_cosets){
  if(c.filter(j=>mn[j]===6).length!==2||!c.includes(0))continue;
  const [u,v]=range(3).filter(i=>c[i]!==0);
  need(maps[u][v][by[c[u]]]===-1&&maps[v][u][by[c[v]]]===-1,"one full Cartesian bank per label");
  maps[u][v][by[c[u]]]=by[c[v]];maps[v][u][by[c[v]]]=by[c[u]];
 }
 const N=labs.map(c=>by[neg(c)]), ordered=[];
 for(let h=0;h<3;h++)for(let e=0;e<3;e++)if(h!==e){
  let f=maps[h][e];need(new Set(f).size===B&&f.every((j,i)=>j>=0&&f[N[i]]===N[j]),"complete matching preserving negation");
  ordered.push(f);
 }
 // All pairwise projector square traces are exact rationals with denominator 54^2.
 const pr=Tr.flat().map(x=>x[0]),prR=Tr.flat().map(x=>x[1]),prg=pr.map(GP);
 const Q=range(B).map(()=>new Int32Array(B));
 for(let i=0;i<pr.length;i++)for(let j=0;j<=i;j++){
  const a=ip(pr[i],prg[j]),b=ip(prR[i],prg[j]),t=(a*a+b*b+(a+b)**2-9)**2;
  Q[i/6|0][j/6|0]+=t;if(i!==j)Q[j/6|0][i/6|0]+=t;
 }
 const Rr=[], seen=new Set();
 for(const r of sh[4])if(!seen.has(r.join(","))){
  const tri=[r,R(r),R(R(r))];for(const v of tri){seen.add(v.join(","));seen.add(v.map(n=>-n).join(","));}
  Rr.push(r);
 }
 need(Rr.length===126,"126 complex root lines");
 const U=range(B).map(()=>new Uint32Array(4)), Bad=range(B).map(()=>[new Uint32Array(4),new Uint32Array(4)]);
 for(const [i,r] of Rr.entries()){
  const s=R(r),bit=(1<<(i%32))>>>0,k=i>>5;
  for(let j=0;j<B;j++){
    if(range(6).every(t=>{const g=prg[j*6+t],a=ip(r,g),b=ip(s,g);return a*a+b*b+(a+b)**2===6;}))U[j][k]=(U[j][k]|bit)>>>0;
    // Pairings with both signs. R-closure means the same failure for all phases of one trine.
    for(const w of bank[labs[j]]){
     const a=ip(w,GP(r));
     if(a>1)Bad[j][0][k]=(Bad[j][0][k]|bit)>>>0;
     if(a< -1)Bad[j][1][k]=(Bad[j][1][k]|bit)>>>0;
    }
  }
 }
 const gh=V.map(GP);const poss=mn.flatMap((m,i)=>m===8?bank[i]:[]);
 need(poss.length===20412,"norm-eight rows"); const handled=new Set();let orbit=0, incidence=0;const allSignPossible={};
 for(const v of poss){
  const A=new Int8Array(B).fill(-8),vs=v;
  for(let i=0;i<gh.length;i++){
   const w=gh[i],p=vs[0]*w[0]+vs[1]*w[1]+vs[2]*w[2]+vs[3]*w[3]+vs[4]*w[4]+vs[5]*w[5]+vs[6]*w[6]+vs[7]*w[7]+vs[8]*w[8]+vs[9]*w[9]+vs[10]*w[10]+vs[11]*w[11];
   const j=i/18|0;if(p>A[j])A[j]=p;
  }
  const ab=range(B).filter(i=>A[i]===5);
  need(ab.length===32&&ab.every(i=>A[N[i]]===4),"32 contacts, paired 4 and 5");
  const key=ab.join(",");if(handled.has(key))continue;handled.add(key);orbit++;
  for(const f of ordered){
   const S=ab.map(j=>f[j]), ls=new Uint32Array(4).fill(0xffffffff),
      fail=[new Uint32Array(4),new Uint32Array(4)];let tr2=0;
   for(const j of S){
    for(let k=0;k<4;k++){ls[k]&=U[j][k];for(let h=0;h<2;h++)fail[h][k]|=Bad[j][h][k];}
    for(const i of S)tr2+=Q[i][j];
   }
   need(tr2===2488320,"projector frame trace square");
   const axes=range(126).filter(i=>ls[i>>5]>>>(i%32)&1);
   need(axes.length===6,"six simultaneous unbiased root lines");
   for(let i=0;i<6;i++)for(let j=0;j<i;j++){
     const x=Rr[axes[i]],y=GP(Rr[axes[j]]);
     need(ip(x,y)===0&&ip(R(x),y)===0,"root complex lines pairwise orthogonal");
   }
   let possible=0;for(const i of axes)for(let s=0;s<2;s++)if(!((fail[s][i>>5]>>>(i%32))&1))possible++;
   need(possible===1,"unscaled body unique oriented root trine");
   incidence++;
  }
 }
 need(orbit===6804&&incidence===40824,"full contact pattern replay");
 // Every retained norm-six fibre has the negative partner used in the real plane proof.
 const out={sourceNorm6:V.length,sourceNorm8:poss.length,allOrderings:ordered.length,contactPatterns:orbit,
  auditedShadowPairs:incidence,projectorTrace:160,projectorTraceSquare:"2488320/2916",
  unbiasedRootLines:6,unscaledAllowedOrientedTrines:1};
 return out;
}
const res=verifyRadial(C,loadCT(C));console.log(JSON.stringify(res));

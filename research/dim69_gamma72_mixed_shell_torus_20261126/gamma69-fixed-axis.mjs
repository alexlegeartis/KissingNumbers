
/* Independent consequence within one actual 2D4 normal space. Run node FILE.
   This regenerates exact arithmetic rather than trusting saved LP values.
   Source premises and geometry are in lower56_69d_discriminant_torus_20261126.md. */
"use strict";
const theta=["1","0","0","0","6218175600","15281788354560","9026867482214400","1989179450818560000","213006159759990870000","13144087517631410995200","525100718690287495741440","14756609779472604266496000","310160311536865273422120000","5108106498702230684344320000","68347693198345293752360140800"];
const makeM=function(thetaAll){
const n=4,S=[8,10,12];const theta=thetaAll;
const q=h=>h.reduce((s,x)=>s+x*x,0);
const H=[[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]];
const hfun=h=>H.map(r=>r.reduce((t,x,i)=>t+x*h[i],0)/2);
const sort=h=>h.map(Math.abs).sort((a,b)=>b-a);
const lex=(a,b)=>{for(let i=0;i<n;i++)if(a[i]!==b[i])return a[i]-b[i];return 0};
function minrep(h){
 let v=h.map(x=>Math.floor((x+4)/8));let r=h.map((x,i)=>x-8*v[i]);let k=v.reduce((a,b)=>a+b,0)&1;
 if(k){let max=-1;for(let i=0;i<n;i++)if(max<0||Math.abs(r[i])>Math.abs(r[max]))max=i;
 r[max]+=r[max]>0?-8:8;}
 return r;
}
function minQ(h){return q(minrep(h));}
function key(h){let v=h.map(x=>Math.floor(x/8));let r=h.map((x,i)=>x-8*v[i]);let k=v.reduce((a,b)=>a+b,0)&1;return r.join()+":"+k;}

function orbit(h){
 const v=h.map(x=>Math.floor((x+4)/8));const r=h.map((x,i)=>Math.abs(x-8*v[i]));
 return [r[0],...r.slice(1).sort((a,b)=>b-a),r.includes(4)?0:v.reduce((a,b)=>a+b,0)&1];
}
let low=new Map(),P=[],M=new Map(),U=[];
function elk(s,h,r,m){
 if(r+m<128)return;
 let rep=minrep(h),ca=orbit(h);
 let id=r+":"+ca.join();if(!M.has(id)){M.set(id,U.length);U.push({r,can:ca,s0:(r+m)/16});}
 let i=M.get(id);P.push({s,h,r,i,m,origin:q(h)===m});
 if(q(h)===m){let k=r+":"+key(h),prev=low.get(k);if(!prev||lex(h,prev.h)<0)low.set(k,{s,h,r,i,m,origin:true});}
}
let R0=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let a of [-8,8])for(let b of [-8,8]){let r=Array(n).fill(0);r[i]=a;r[j]=b;R0.push([8,r]);}
const p0=Array(n).fill(0);
for(let s of S){
function rec(k,a,b,parity){
 if(k===n){
  if(a>16*s)return;
  if(a===16*s){if(p0.every(x=>x%8===0)&&p0.reduce((x,y)=>x+y,0)%16===0)R0.push([s,p0.slice()]);return;}
  let h=p0.slice(),r=16*s-a,m=minQ(h);elk(s,h,r,m);return;
 }
 for(let t=-s;t<=s;t++)if((k===0||((t&1)===parity))&&Math.abs(t)+b<=s){
  p0[k]=t;rec(k+1,a+t*t,Math.max(b,Math.abs(t)),t&1);
 }
}
rec(0,0,0,0);
}
let C=[...low.values()].filter(z=>z.r<=128);
let ps=[];
for(let t=0;t<=5;t++)for(let i=0;i<=t;i++)for(let j=0;j<=t-i;j++)for(let k=0;k<=j;k++){let ell=t-i-j-k;if(ell>=0&&ell<=k)ps.push([i,j,k,ell]);}
const roots=R0.map(x=>x[1]);const A=[],B=[],raw=[];
const f=(h,p)=>h.reduce((x,v,j)=>x*v**(2*p[j]),1);
for(let s of S)for(let p of ps){
 let row=U.map(()=>0);for(let z of P)if(z.s===s)row[z.i]+=f(z.h,p);
 let l=p.reduce((a,b)=>a+b,0);let num=BigInt(theta[s/2])*(16n*BigInt(s))**BigInt(l),den=1n;
 for(let a of p)for(let j=1;j<=a;j++)num*=BigInt(2*j-1);
 for(let j=0;j<l;j++)den*=BigInt(72+2*j);
 if(s===8)for(let h of roots)num-=den*BigInt(f(h,p));
 A.push(row);B.push(Number(num)/Number(den));raw.push({s,p,num:String(num),den:String(den)});
}
return {U,P,C,A,B,raw,roots,theta,low:R0, D4:true};
};
const harm=function(model){
const factor=[1n];for(let i=1;i<=90;i++)factor.push(factor.at(-1)*BigInt(i));
const gcd=(a,b)=>{if(a<0n)a=-a;if(b<0n)b=-b;while(b){let c=a%b;a=b;b=c;}return a;};
const lcm=(a,b)=>a/gcd(a,b)*b;
const dim=72;
const a0=model.A.map(r=>r.map(BigInt));
const b0=model.raw.map(z=>{let n=BigInt(z.num),d=BigInt(z.den);if(n%d)throw Error("bad exact");return n/d;});
const N=model.C.length;
const Ushell=[8,10,12];
const baseTables=[];for(let ss of Ushell){
 let I=new Map();for(let z of model.P.filter(p=>p.s===ss)){let k=z.h.map(Math.abs).join();let t=k+":"+z.i;if(!I.has(t))I.set(t,[z.h.map(Math.abs),z.i,0n]);I.get(t)[2]++;}
 baseTables.push([...I.values()]);
}
const cache=new Map();
function mom(sh,exps){
 let p=[exps[0],...exps.slice(1).sort((a,b)=>b-a)],key=sh+":"+p.join();if(cache.has(key))return cache.get(key);
 let d=model.U.map(()=>0n),e=0n;
 for(let [h,i,m] of baseTables[sh])d[i]+=m*h.reduce((s,x,j)=>s*BigInt(x)**BigInt(2*p[j]),1n);
 if(sh===0)for(let h of model.roots)e+=h.reduce((s,x,j)=>s*BigInt(x)**BigInt(2*p[j]),1n);
 let r=[d,e];cache.set(key,r);return r;
}
const M3=7;const one=Array.from({length:M3},(_,i)=>i?0n:1n);
const mul=(a,b)=>a.map((_,i)=>a.slice(0,i+1).reduce((z,v,j)=>z+v*b[i-j],0n));
const pow=(a,e)=>{let r=one;for(let i=0;i<e;i++)r=mul(a,r);return r};
const sig=(d,a)=>Array.from({length:d},(_,i)=>i+1).filter(i=>d%i===0).reduce((s,i)=>s+BigInt(i)**BigInt(a),0n);
const ex=a=>one.map((_,i)=>i?(a===3?240n:-504n)*sig(i,a):1n), E4=ex(3),E6=ex(5);
let Del=one.map((_,i)=>i===1?1n:0n);for(let j=1;j<7;j++){let a=one.map((x,i)=>i===j?-1n:x);Del=mul(Del,pow(a,24));}
function nullModular(m){
 const rows=[];let present=[];
 for(let i=0;i<=2;i++){
  let w=2*m-12-12*i,b=(w%4+4)%4===2?1:0,a=(w-6*b)/4;
  if(a>=0){let y=mul(pow(Del,4+i),mul(pow(E4,a),pow(E6,b))).slice(4,7);present.push([i,y]);}
 }
 for(let i=0;i<=2;i++)if(!present.some(x=>x[0]===i)){
  let v=[0n,0n,0n];v[i]=1n;
  for(let [j,y] of present.slice().reverse())if(j<i)v[j]=-(v[0]*y[0]+v[1]*y[1]+v[2]*y[2]);
  rows.push(v);
 }
 return rows;
}
const modularRows=[];
let meta=[];
for(let m=6;m<=24;m++){
 let lambdas=nullModular(m);if(lambdas.length===0)continue;
 let numer=[1n],denom=[1n];
 for(let i=1;i<=m;i++){
  let nu=-numer[i-1]*BigInt((2*m-2*i+2)*(2*m-2*i+1));
  let de=denom[i-1]*BigInt(2*i*(4*m-2*i+dim-2));
  let d=gcd(nu,de);numer.push(nu/d);denom.push(de/d);
 }
 let cden=denom.reduce(lcm,1n),cc=numer.map((z,i)=>z*cden/denom[i]);
 for(let a=m;a>=0;a--)for(let b=m-a;b>=0;b--)for(let c=Math.min(b,m-a-b);c>=0;c--){let d=m-a-b-c;if(d<0||d>c)continue;
  let rr=[],qq=[];
  for(let sh=0;sh<3;sh++){
   let arr=model.U.map(()=>0n), rhs=0n,s=Ushell[sh];
   let pre=cc.map((co,i)=>co*BigInt(16*s)**BigInt(i)*factor[i]*factor[2*m-2*i]);
   for(let ia=0;ia<=a;ia++)for(let ib=0;ib<=b;ib++)for(let ic=0;ic<=c;ic++)for(let id=0;id<=d;id++){
    let i=ia+ib+ic+id,pp=[a-ia,b-ib,c-ic,d-id];
    let w=pre[i]/(factor[ia]*factor[ib]*factor[ic]*factor[id]*pp.reduce((t,x)=>t*factor[2*x],1n));
    let [ro,re]=mom(sh,pp);
    for(let k=0;k<arr.length;k++)arr[k]+=w*ro[k];
    rhs+=w*re;
   }
   rr.push(arr);qq.push(rhs);
  }
  for(let z of lambdas){
   let row=rr[0].map((x,k)=>x*z[0]+rr[1][k]*z[1]+rr[2][k]*z[2]);let rhs=-(z[0]*qq[0]+z[1]*qq[1]+z[2]*qq[2]);let d=row.reduce(gcd,gcd(rhs,0n));if(d===0n)continue;
   modularRows.push([row.map(x=>x/d),rhs/d]);meta.push([m,a,b,c,m-a-b-c,z.map(String)]);
  }
 }
}
const A=[...a0,...modularRows.map(x=>x[0])],B=[...b0,...modularRows.map(x=>x[1])];
return {A:A.map(r=>r.map(String)),B:B.map(String),meta,modularCount:modularRows.length};
};
const mod=function(harm){
let A=harm.A.map(r=>r.map(BigInt)),B=harm.B.map(BigInt);const m=A.length,n=A[0].length;
function lap(p){p=BigInt(p);const res=x=>{const y=x%p;return y<0n?y+p:y}; const inv=x=>{let y=p-2n,s=1n;for(;y;y>>=1n,x=x*x%p)if(y&1n)s=s*x%p;return s;};
let W=A.map((r,i)=>[...r.map(res),res(48n*B[i])]);let labels=W.map((_,i)=>i);let ids=[],J=[];
for(let j=0;j<n;j++){let i=W.findIndex((row,i)=>i>=ids.length&&row[j]);if(i<0)continue;let rk=ids.length;[W[i],W[rk]]=[W[rk],W[i]];[labels[i],labels[rk]]=[labels[rk],labels[i]];J.push(labels[rk]);let z=inv(W[rk][j]);for(let k=j;k<=n;k++)W[rk][k]=W[rk][k]*z%p;
for(let v=0;v<W.length;v++){if(v===rk||!W[v][j])continue;let x=W[v][j];for(let k=j;k<=n;k++)W[v][k]=res(W[v][k]-x*W[rk][k]);}ids.push(j);}
let consistent=W.slice(ids.length).every(x=>!x[n]);
const fixed=ids.filter((j,i)=>W[i].slice(0,n).every((x,k)=>!x||ids.includes(k)));
return {p,ids,J,W,consistent,fixed,fv:fixed.map(j=>W[ids.indexOf(j)][n])};}
const out=[1000003].map(lap);
return {out:out.map(x=>({...x,p:String(x.p),W:x.W.map(row=>row.map(String)),fv:x.fv.map(String)}))};
};
const affine=function(harm,rr){
const F=(n,d=1n)=>{if(n===0n)return [0n,1n];if(d<0n){d=-d;n=-n;}let z=G(n,d);return [n/z,d/z]};
const G=(n,d)=>{if(n<0n)n=-n;if(d<0n)d=-d;while(d){let a=n%d;n=d;d=a;}return n};
const plus=([a,b],[c,d])=>F(a*d+b*c,b*d);
const sub=([a,b],[c,d])=>F(a*d-b*c,b*d);
const mul=([a,b],[c,d])=>F(a*c,b*d);
const div=([a,b],[c,d])=>F(a*d,b*c);
let A=harm.A.map(r=>r.map(BigInt)),B=harm.B.map(BigInt);const n=A[0].length;
const rank=rr.ids.length, ix=rr.ids.concat(Array.from({length:n},(_,i)=>i).filter(i=>!rr.ids.includes(i)));
const G0=rr.J.map(j=>ix.map(i=>A[j][i]).concat(B[j]));
let prev=1n,maxl=0;
for(let i=0;i<rank;i++){
 let k=G0.findIndex((z,j)=>j>=i&&z[i]);if(k<0)throw Error("sing "+i);
 [G0[k],G0[i]]=[G0[i],G0[k]];let pivot=G0[i][i];
 for(let j=i+1;j<rank;j++){
  let x=G0[j][i];
  for(let k=i+1;k<=n;k++){
   let t=G0[j][k]*pivot-x*G0[i][k];
   if(t%prev!==0n)throw Error("bareiss");
   G0[j][k]=t/prev;
  }
  G0[j][i]=0n;
 }
 prev=pivot;
}
let gsol=Array.from({length:n},(_,i)=>Array.from({length:n-rank+1},(_,j)=>F(i>=rank&&j===i-rank?1n:0n)));
for(let j=rank-1;j>=0;j--){
 for(let k=0;k<=n-rank;k++){
  let x=k===n-rank?[G0[j][n],1n]:[0n,1n];
  for(let i=j+1;i<n;i++)if(G0[j][i])x=sub(x,mul([G0[j][i],1n],gsol[i][k]));
  gsol[j][k]=div(x,[G0[j][j],1n]);
 }
}
const s=Array.from({length:n},(_,i)=>gsol[ix.indexOf(i)]), zero=[0n,1n];
for(let i=0;i<A.length;i++){
  for(let k=0;k<=n-rank;k++){
   let c=s.reduce((z,x,j)=>A[i][j]?plus(z,mul(x[k],[A[i][j],1n])):z,zero);
   if(c[1]!==1n||c[0]!== (k===n-rank?B[i]:0n))throw Error("bad "+JSON.stringify([i,k,c.map(String)]));
 }}
return {s:s.map(r=>r.map(x=>x.map(String))),free:ix.slice(rank)};
};
const exactB=function(e){
const p=(a,b)=>a+b;
const G=(n,d)=>{if(n<0n)n=-n;if(d<0n)d=-d;while(d){let a=n%d;n=d;d=a;}return n};
const F=(n,d=1n)=>{if(n===0n)return [0n,1n];if(d===0n)throw Error("0-den");if(d<0n){d=-d;n=-n;}let z=G(n,d);return [n/z,d/z]};
const plus=([a,b],[c,d])=>F(a*d+b*c,b*d);
const sub=([a,b],[c,d])=>F(a*d-b*c,b*d);
const mul=([a,b],[c,d])=>F(a*c,b*d);
const div=([a,b],[c,d])=>F(a*d,b*c);
const S=e.s.map(r=>r.map(x=>x.map(BigInt))),n=e.free.length;
function exactValue(j,sgn,act){
 const a=S[j].map(([x,d])=>[BigInt(sgn)*x,d]); // I want proof n_j for lower, -n_j for upper.
 let M=Array.from({length:n},(_,i)=>act.map(r=>S[r][i]).concat([a[i]]));
 let k=0,J=[];for(let i=0;i<act.length;i++){
  let j=M.findIndex((r,p)=>p>=k&&r[i][0]);if(j<0)continue;
  [M[k],M[j]]=[M[j],M[k]];const x=M[k][i];M[k]=M[k].map(z=>div(z,x));
  for(let j=0;j<M.length;j++)if(j!==k&&M[j][i][0]){
   let z=M[j][i];M[j]=M[j].map((x,l)=>sub(x,mul(z,M[k][l])));
  }
  J.push(i);k++;
 }
 if(!M.slice(k).every(row=>!row.at(-1)[0]))return null;
 let q=act.map((_,i)=>J.includes(i)?M[J.indexOf(i)].at(-1):[0n,1n]);
 if(q.some(x=>x[0]<0n))return null;
 let dif=a.map((v,i)=>sub(v,act.reduce((s,j,k)=>plus(s,mul(q[k],S[j][i])),[0n,1n])));
 if(dif.slice(0,n).some(v=>v[0]))return null;
 let bound=dif[n];if(sgn===-1)bound=[-bound[0],bound[1]];
 return {i:j,sgn,bound:bound.map(String),expr:act.map((k,i)=>[k,q[i].map(String)]).filter(x=>x[1][0]!=="0")};
}
return exactValue;
};
function need(x,w){if(!x)throw Error(w);}
const m=makeM(theta), H=harm(m);
const col=m.U.map((_,i)=>i).sort((i,j)=>m.U[j].s0-m.U[i].s0||m.U[j].r-m.U[i].r);
const rr=mod({A:H.A.map(r=>col.map(i=>r[i])),B:H.B}).out[0];
need(m.U.length===162&&rr.consistent&&rr.ids.length===149,"fixed-axis rank");
const e=affine(H,{...rr,ids:rr.ids.map(i=>col[i])});
need(e.free.join()==="0,1,2,3,6,7,10,11,12,24,34,41,44","free basis");
need(m.U[0].s0===8&&m.U[54].s0===10&&m.U[108].s0===12&&[0,54,108].every(i=>m.U[i].can.join()==="0,0,0,0,1"),"actual axial hole cells");
need(m.U[5].s0===8&&m.U[5].r===92&&m.U[5].can.join()==="2,0,0,0,1","actual norm8 wing-difference cell");
need(e.s[54].every(x=>x[1]==="1")&&e.s[54].map(x=>x[0]).join()==="-64,-64,16,0,-64,-48,0,0,0,0,-16,0,0,49664","deep empty implication");
const G=(a,b)=>{if(a<0n)a=-a;if(b<0n)b=-b;while(b){let c=a%b;a=b;b=c;}return a};
const F=(n,d=1n)=>{let a=G(n,d);return [n/a,d/a]};
const plus=([a,b],[c,d])=>F(a*d+b*c,b*d);
const mul=([a,b],[c,d])=>F(a*c,b*d);
const cols=e.free.filter(i=>i!==0&&i!==2);
const ev=e.s.map(r=>r.map(x=>x.map(BigInt)));
const shift={1:4,6:4,7:3,34:1};
const S={free:cols,s:ev.map(r=>[...cols.map(i=>plus(r[e.free.indexOf(i)],mul(r[2],[BigInt(shift[i]??0),1n]))),plus(r.at(-1),mul(r[2],[-3104n,1n]))].map(x=>x.map(String)))};
const solve=exactB(S);
const cuts=[{"i":5,"sgn":1,"bound":["3264","1"],"expr":[[2,["16","1"]],[42,["6","1"]],[3,["12","1"]],[41,["8","1"]]]},{"i":33,"sgn":1,"bound":["629070656","5"],"expr":[[2,["7664","5"]],[4,["288","5"]],[13,["2448","5"]],[20,["256","5"]],[25,["1728","5"]],[42,["498","5"]],[64,["288","5"]],[3,["11124","5"]],[11,["24768","5"]],[24,["10296","5"]],[41,["3928","5"]]]},{"i":53,"sgn":1,"bound":["366016472","1"],"expr":[[2,["640","1"]],[1,["2048","1"]],[3,["96","1"]],[6,["2048","1"]],[10,["8784","1"]],[11,["1728","1"]],[12,["216","1"]],[24,["576","1"]],[41,["32","1"]],[44,["72","1"]]]},{"i":108,"sgn":1,"bound":["1825685184","5"],"expr":[[2,["3456","5"]],[4,["768","5"]],[20,["512","5"]],[42,["768","5"]],[3,["96","1"]],[10,["8784","1"]],[11,["1728","1"]],[12,["216","1"]],[24,["576","1"]],[41,["32","1"]],[44,["72","1"]]]}];
for(const old of cuts){
 const z=solve(old.i,old.sgn,old.expr.map(x=>x[0]));
 need(JSON.stringify(z)===JSON.stringify(old),"per-axis rational cut");
}
const A=solve(5,1,[2,42,3,41]);
need(A.bound[0]==="3264"&&A.bound[1]==="1","critical nonempty fiber");
// Independent proof inspection: these are equality coefficients, not just a dual optimization outcome.
need(A.expr.map(([i,[a,b]])=>i+":"+a+"/"+b).join()==="2:16/1,42:6/1,3:12/1,41:8/1","critical identity");
const ret={classes:162,rank:149,scope:"actual deep-hole and fixed-axis populations; remaining auxiliary counts may be orientation averages",
 simultaneous_absence_at_s8_s10_implies:{minimum_norm8_axis_difference:3264,certificate:A.expr},
 conditional_actual_exact_anchors:[366016472,365137037],
 conditional_actual_wing:125814132,
 use:"The roughly 857M universal construction using two radial-8 anchors and a fixed axial wing fails its required third shell absence; the theorem does not establish that conflict endpoints are actually both present."};
if(typeof report==="function")report(ret);else console.log(JSON.stringify(ret,null,2));


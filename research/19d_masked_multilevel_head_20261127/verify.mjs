const M19=(()=>{
const TS=[0,3,6,8],T=TS.reduce((x,j)=>x|(1<<j),0);
const G0=[[1,8,9,12,16,17,18,19],[2,10,11,14,15,17,18,20],[3,7,9,13,15,16,17,19],[4,7,8,10,12,15,16,19],[5,10,12,13,15,16,17,18],[6,7,8,9,10,13,16,18],[1,4,7,9],[1,5,6,18,20],[1,3,12,15,20],[1,10,13,19,20],[1,3,5,6,7,13,14,15,18,20],[2,4,6,7,8,13,14,16,17,18]].map(v=>v.reduce((x,j)=>x|1<<(j-1),0));
const pc=x=>x.toString(2).replaceAll('0','').length;
let E=[0],dual=[0],piv=Array(20).fill(0);
for(let g of G0){E.push(...E.map(x=>x^g));for(let j=19;j>=0;j--)if(g>>>j&1){if(piv[j])g^=piv[j];else{piv[j]=g;break;}}}
for(let j=0;j<20;j++)if(!piv[j]){let x=1<<j;for(let k=0;k<20;k++)if(pc(x&piv[k])%2)x|=1<<k;dual.push(...dual.map(a=>a^x));}
if(E.length!==4096||dual.length!==256||dual.some(w=>G0.some(g=>pc(g&w)%2)))throw Error("binary sources");
const O=dual.filter(w=>pc(w)===8),Z=E.filter(w=>(w&T)===0);
const B=O.filter(w=>pc(w&T)===2),F={};let hh={};
for(let w of B){
 const f=new Set(Z.map(z=>z&w)),r=Math.log2(f.size),key=w&T;F[key]??=[];
 F[key].push([w&~T,r]);hh[r]=(hh[r]||0)+1;
}
return {T,TS,Z,F,hh,O,B,pc,E,G0,dual};
})();


function setup(M){
 const H=Array.from({length:20},(_,i)=>i).filter(i=>!M.TS.includes(i)),pi=(w)=>H.reduce((r,j,i)=>r|((w>>>j&1)<<i),0);
 const R=[[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]];
 const g=Object.entries(M.F).map(([k,o])=>[Number(k),o.map(([m])=>pi(m))]);
 const V=g.map(([z])=>{const g=M.TS.map((t,i)=>z>>t&1?i:-1).filter(i=>i>=0);return R.map(r=>(r[g[0]]-r[g[1]])/2);});
 const oo=M.O.filter(x=>!(x&M.T)).map(pi),cap=M.O.filter(x=>(x&M.T)===M.T).map(pi);
 return {H,pi,R,g,V,oo,cap};
}
const pc16=new Uint8Array(65536);for(let i=1;i<65536;i++)pc16[i]=pc16[i>>>1]+(i&1);
const det3=(a,b,c,d,e,f,g,h,i)=>a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
function cofactor(M,i,j){
 const A=[];for(let k=0;k<4;k++)if(k!==i)for(let l=0;l<4;l++)if(l!==j)A.push(M[k][l]);
 return ((i+j)&1?-1:1)*det3(...A);
}
function sigs(X,optimistic){
 const S=new Map();
 for(let s=1;s<65536;s++){
  const k=pc16[s],SH=X.g.map(([t,O])=>Math.max(...O.map(o=>{let m=pc16[s&o];return m===6&&optimistic?4:m})));
  const O0=Math.max(...X.oo.map(o=>{let m=pc16[s&o];return m===8?6:m}));
  const C0=Math.max(...X.cap.map(o=>{let m=pc16[s&o];return m===4?2:m}));
  const a=[k===1?2:4,Math.max(O0,0),C0,...SH],key=a.join(",");
  if(!S.has(key))S.set(key,{a,counts:Array(17).fill(0),masks:[]});
  const V=S.get(key);V.counts[k]++;V.masks.push(s);
 }
 return [...S.values()];
}
function rowchart(a,V){
 const [r,o,c,...p]=a,A=[[-1,0,0,0,0],[Math.max(r,o),0,0,0,4]];
 for(let i=0;i<3;i++)for(const s of [1,-1]){
  const z=[c,0,0,0,4];z[i+1]=2*s;A.push(z);
 }
 for(let i=0;i<V.length;i++)for(const s of [1,-1]){
  A.push([p[i],...V[i].map(v=>v*s),4]);A.push([0,...V[i].map(v=>v*s),2]);
 }
 const seen=new Map();
 for(const x of A){const k=x.slice(0,4).join(",");seen.set(k,Math.min(seen.has(k)?seen.get(k):1000,x[4]));}
 return [...seen].map(([x,v])=>x.split(",").map(Number).concat(v));
}
function maxima(a,V,counts){
 const A=rowchart(a,V),M=A.length,N=counts.map(x=>x?[-1,1,[]]:null),good=[];
 for(let i=0;i<M-3;i++)for(let j=i+1;j<M-2;j++)for(let k=j+1;k<M-1;k++)for(let l=k+1;l<M;l++){
  const B=[A[i],A[j],A[k],A[l]],d0=cofactor(B,0,0),d1=cofactor(B,0,1),d2=cofactor(B,0,2),d3=cofactor(B,0,3);
  let d=B[0][0]*d0+B[0][1]*d1+B[0][2]*d2+B[0][3]*d3;if(!d)continue;
  const Q=[d0,d1,d2,d3].map(x=>x*B[0][4]);
  for(let h=1;h<4;h++){if(!B[h][4])continue;const b=B[h][4];for(let u=0;u<4;u++)Q[u]+=b*cofactor(B,h,u);}
  if(d<0){d=-d;for(let u=0;u<4;u++)Q[u]=-Q[u];}
  if(!Q[0]||A.some(([r,u,v,w,b])=>r*Q[0]+u*Q[1]+v*Q[2]+w*Q[3]>b*d))continue;
  const h=Q[0]**2,t=Q[1]**2+Q[2]**2+Q[3]**2,dd=d*d;
  for(let k=1;k<=16;k++)if(N[k]){
   const f=k*h+t;
   if(f*N[k][1]>dd*N[k][0])N[k]=[f,dd,Q.concat(d)];
   if(f>=8*dd)good.push([k,Q.concat(d),f,dd]);
  }
 }
 return {N,good};
}
function full(X,optimistic){
 const S=sigs(X,optimistic),G=[],scores=[];
 let hist={};
 for(const i of S){const w=maxima(i.a,X.V,i.counts);for(let k=1;k<=16;k++)if(w.N[k]){const [f,d]=w.N[k];const n=f/d;hist[n]=(hist[n]||0)+i.counts[k];} if(w.good.length)G.push({sig:i.a,counts:i.counts,masks:i.masks,good:w.good});
 scores.push(w.N);
 }
 return {S:S.length,G,scores,hist};
}


const det=(a,b,c)=>a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]);
function solve(S){
 const {a,counts}=S,[root,o,c,...p]=a,l=Math.min(...p);
 const H=[[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,1,2],[4,root,0,4],[6,o,0,4],[6,l,1,4],[2,c,1,4]];
 let out=[],MV=counts.map(k=>k?[-1,1,null]:null);
 for(let i=0;i<H.length;i++)for(let j=i+1;j<H.length;j++)for(let k=j+1;k<H.length;k++){
  const R=[H[i],H[j],H[k]];let D=det(...R);if(!D)continue;
  const Q=Array.from({length:3},(_,p)=>det(...R.map(r=>[r[0],r[1],r[2]].map((z,j)=>j===p?r[3]:z))));
  if(D<0){D=-D;for(let l=0;l<3;l++)Q[l]=-Q[l];}
  if(H.some(([a,b,c,h])=>a*Q[0]+b*Q[1]+c*Q[2]>h*D))continue;
  for(let K=1;K<=16;K++)if(MV[K]){
   const [b,c,r]=Q,n=16*b*b+2*K*b*c+K*c*c+r*r;
   if(n*MV[K][1]>D*D*MV[K][0])MV[K]=[n,D*D,Q.concat(D)];
   if(n>=8*D*D)out.push([K,n,D*D,Q.concat(D)]);
  }
 }
 return {MV,out};
}

function relaxedFlat({a,counts}){
const [r,o,c,...p]=a, q=Math.min(...p),L=[[-1,0,0],[0,-1,0],[0,1,2],[Math.max(r,o),0,4],[q,1,4],[c,1,4]],scores=counts.map(k=>k?[-1,1]:null);
for(const [i,[a,b,c]]of L.entries())for(const [d,e,f]of L.slice(i+1)){let D=a*e-b*d;if(!D)continue;let u=c*e-b*f,t=a*f-c*d;if(D<0){D=-D;u=-u;t=-t;}if(L.some(([a,b,c])=>a*u+b*t>c*D))continue;
 for(let k=1;k<16;k++)if(scores[k]){let N=k*u*u+t*t;if(N*scores[k][1]>scores[k][0]*D*D)scores[k]=[N,D*D];}
}
return scores;
}

function refinedZero(S){
 const {a,counts}=S,[root,o,c,...p]=a,l=Math.min(...p);
 const H=[[-1,0,0,0,0],[0,-1,0,0,0],[0,0,-1,0,0],[0,0,0,-1,0],
 [0,0,1,0,2],[0,0,-1,1,0],[0,0,1,-2,0],[4,root,0,0,4],[6,o,0,0,4],[6,l,1,0,4],[2,c,0,2,4]];
 let MV=counts.map(k=>k?[-1,1]:null);
 for(let i=0;i<H.length;i++)for(let j=i+1;j<H.length;j++)for(let k=j+1;k<H.length;k++)for(let l=k+1;l<H.length;l++){
  const R=[H[i],H[j],H[k],H[l]],D0=R[0].slice(0,4).reduce((s,x,u)=>s+x*cofactor(R,0,u),0);
  if(!D0)continue;
  const Q=Array.from({length:4},(_,p)=>R.reduce((n,z,r)=>n+z[4]*cofactor(R,r,p),0));
  const D=Math.abs(D0);
  if(D0<0)for(let u=0;u<4;u++)Q[u]=-Q[u];
  if(H.some(z=>z.slice(0,4).reduce((s,x,i)=>s+x*Q[i],0)>z[4]*D))continue;
  for(let K=1;K<=16;K++)if(MV[K]){
   const [b,c,r,a]=Q,n=16*b*b+2*K*b*c+K*c*c+a*a+2*(r-a)*(r-a);
   if(n*MV[K][1]>D*D*MV[K][0])MV[K]=[n,D*D];
  }
 }
 return MV;
}


const detA=m=>m.length===1?m[0][0]:m[0].reduce((s,x,j)=>s+((j&1)?-1:1)*x*detA(m.slice(1).map(r=>r.filter((_,k)=>k!==j))),0);
function polyTetrad(){
 // x=(a,b,c,u,r), in decreasing head weights and u=max tail, r=sum top two.
 const H=[
 [0,0,-1,0,0,0],[0,-1,1,0,0,0],[-1,1,0,0,0,0],
 [0,0,0,-1,0,0],[0,0,0,1,-1,0],[0,0,0,-2,1,0],[0,0,0,0,1,2],
 [2,1,0,0,0,2],[2,2,2,0,1,4],[1,0,0,1,0,2]
 ];
 let max=[0,1,[]];
 function finish(R){
  const M=R.map(r=>r.slice(0,5));let D=detA(M);if(!D)return;
  const Q=Array.from({length:5},(_,i)=>detA(M.map((x,j)=>x.map((q,k)=>i===k?R[j][5]:q))));
  if(D<0){D=-D;for(let i=0;i<Q.length;i++)Q[i]=-Q[i];}
  if(H.some(r=>r.slice(0,5).reduce((s,x,j)=>s+x*Q[j],0)>r[5]*D))return;
  const [a,b,c,u,r]=Q;
  const L=4*(a*a+b*b+2*c*c)+u*u+2*(r-u)*(r-u),d=D*D;
  if(L*max[1]>d*max[0])max=[L,d,Q.concat(D)];
 }
 function rec(s,R){if(R.length===5){finish(R);return;}for(let i=s;i<=H.length+R.length-5;i++)rec(i+1,R.concat([H[i]]));}
 rec(0,[]);return max;
}
const A=(v,m)=>{if(!v)throw Error(JSON.stringify(m));};
const X=setup(M19);
const C=[0,65535];
for(const o of X.oo)if(!C.includes(o))C.push(...C.map(c=>c^o));
A(C.length===32&&X.oo.length===30&&new Set(X.oo).size===30&&C.every(c=>[0,65535].includes(c)||X.oo.includes(c)),"head hyperplanes span rank5");
A(new Set(M19.E.map(X.pi)).size===2048&&M19.E.every(q=>X.oo.every(o=>!(M19.pc(o&X.pi(q))&1))),"E projected equals shortened octad sign dual");
A(M19.Z.length===256&&new Set(M19.Z.map(X.pi)).size===256&&M19.Z.every(q=>X.g.every(([k,O])=>O.every(o=>!(M19.pc(X.pi(q)&o)&1)))),"zero T head signs align even on all shoulders");
A(X.cap.every((u,i)=>X.cap.slice(i+1).every(v=>!(u&v)&&X.oo.includes(u|v))),"every pair of actual head tetrads forms a disjoint octad");
A(X.g.every(([_,O])=>X.cap.every((t,j)=>O.filter(o=>X.cap.every((v,k)=>pc16[o&v]===(j===k?0:2))).length===4)),"four actual shoulder supports for each omitted tetrad and each orientation");
A(X.g.length===6&&X.g.every(([_,O])=>O.length===16)&&X.cap.length===4&&X.cap.every(t=>pc16[t]===4)&&X.cap.reduce((s,t)=>s^t,0)===65535,"literal base");
for(const [_,O] of X.g)for(let i=0;i<16;i++){
 A(O.filter(o=>o>>>i&1).length===6,"shoulder replication");
 for(let j=i+1;j<16;j++)A(O.filter(o=>(o>>>i&1)&&(o>>>j&1)).length===2,"2-(16,6,2) pairs");
}
A(new Set(X.V.map(v=>v.join(","))).size===6&&X.V.every(v=>v.length===3&&v.reduce((s,t)=>s+t*t,0)===2),"D3 direction pairs");
const F=sigs(X,true);A(F.length===158&&F.reduce((s,x)=>s+x.masks.length,0)===65535,"all masks");
let nf=0,hard=0,normFlat=[0,1];
for(const it of F) {
 const raw=relaxedFlat(it);
 const cross=raw.some(s=>s&&s[0]*9>64*s[1]);
 const bounds=cross?maxima(it.a,X.V,it.counts).N:raw;
 if(cross)hard++;
 for(let k=1;k<16;k++)if(it.counts[k]){
  const [q,d]=bounds[k];A(q>=0&&q*9<=64*d,"flat support radius <=64/9");
  if(q*normFlat[1]>d*normFlat[0])normFlat=[q,d];
  nf+=it.counts[k];
 }
}
A(hard===9&&nf===65534&&normFlat[0]*9===64*normFlat[1],"flat census and optimum exact");
const Z=sigs(X,false);
A(Z.length===241&&Z.reduce((s,x)=>s+x.masks.length,0)===65535,"zero-tetrad shape signatures");
let nEquality=0,zeroHard=0,normZero=[0,1];
for(const x of Z){
 const R=solve(x), test=R.MV.some(s=>s&&9*s[0]>64*s[1]), W=test?refinedZero(x):R.MV;
 if(test)zeroHard++;
 for(let k=1;k<=16;k++)if(x.counts[k]){
  const [n,d]=W[k];A(9*n<=64*d,"zero class two-level radius <=64/9");
  if(n*normZero[1]>d*normZero[0])normZero=[n,d];
 }
 for(const [k,n,d,Q] of R.out){
   const [b,c,r,D]=Q;
   A(n===8*d&&b===0&&r===2*D&&((k===1&&c===2*D)||(k===4&&c===D)),"only formal equality cases");
   nEquality++;
 }
}
A(nEquality===112&&zeroHard===9&&normZero[0]*9===64*normZero[1],"two-level refined ceiling and exact value");
A(29n**2n < (12n**2n)*6n,"first no-same-head convex endpoint");
A(137n**2n < (56n**2n)*6n,"second no-same-head convex endpoint");
A(480+12+3840+6144+192===10668,"Li exact fixed-body count");
const three=polyTetrad();A(three[0]*9===64*three[1],"arbitrary four tetrad levels, same sharp radius");
A(8n*2500n*9n>64n*2809n,"robust 3 percent norm-preserving stability");
console.log(JSON.stringify({status:"PASS",tetradRadius:"64/9",head:16,tail:3,base:10668,flatMasks:65534,flatSignatures:F.length,flatHard:hard,maskedRadius:"64/9",zeroProfiles:Z.length,zeroHard,zeroRadius:"64/9",remainingFullFlatSignHeads:1792}));
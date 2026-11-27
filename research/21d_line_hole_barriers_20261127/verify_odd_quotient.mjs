// Exact 21D Li order-four odd quotient local barrier, Node or V8.

const lineLists=[[4,5,6,8,10],[2,3,6,7,15],[1,4,9,14,15],[3,9,10,11,16],[1,6,12,13,16],[7,8,9,12,17],[4,7,11,13,18],[2,10,12,14,18],[1,3,5,17,18],[2,5,9,13,19],[6,11,14,17,19],[8,15,16,18,19],[1,2,8,11,20],[5,7,14,16,20],[10,13,15,17,20],[3,4,12,19,20],[3,8,13,14,21],[5,11,12,15,21],[2,4,16,17,21],[1,7,10,19,21],[6,9,18,20,21]];
const tags=[11,13,8,15,3,1,12,5,10,4,7,1,2,3,2,1,9,6,14,3,2];
const lineMasks=lineLists.map(x=>x.reduce((s,j)=>s|(1<<(j-1)),0));
const H=[0], eta=[0], M=new Map([[0,0]]);
for(let i=0;i<21;i++){let x=lineMasks[i];if(!M.has(x)){let n=H.length;for(let j=0;j<n;j++){M.set(H[j]^x,j+n);H.push(H[j]^x);eta.push(eta[j]^tags[i]);}}}
if(H.length!==1024||lineMasks.some((x,i)=>eta[M.get(x)]!==tags[i]))throw Error("geometry rank/tag");
const L=lineMasks.map(x=>M.get(x)), holes=H.flatMap((x,i)=>eta[i]===0?[i]:[]);
const LC=L.filter((x,i)=>tags.filter(y=>y===tags[i]).length===1);
if(LC.length!==12||holes.length!==64)throw Error("phase");
const state={};
const pc=x=>{x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0xf0f0f0f)*0x1010101)>>>24};
const subspaces=[];
for(let i=1;i<64;i++)for(let j=i+1;j<64;j++){let a=holes[i],b=holes[j],c=a^b;if(c<b)continue;subspaces.push([0,a,b,c]);}
if(subspaces.length!==651)throw Error("rank2 subspaces");
function qdata(W){
 const Q=new Uint16Array(1024);Q.fill(65535);let R=[];
 for(let i=0;i<1024;i++)if(Q[i]===65535){for(let d of W)Q[i^d]=R.length;R.push(i);}
 const old=R.flatMap((x,j)=>eta[x]===0?[j]:[]), newv=R.flatMap((x,j)=>eta[x]>0?[j]:[]);
 const nb=R.map(x=>[...new Set([Q[x],...L.map(g=>Q[x^g])])]);
 const P=old.map(i=>[i,...LC.map(g=>Q[R[i]^g])]);
 const A=new Uint16Array(newv.length*old.length), D=new Uint16Array(newv.length*old.length);
 for(let i=0;i<newv.length;i++){const a=nb[newv[i]],u=new Set(a);for(let j=0;j<old.length;j++){
  let q=0;for(let k=0;k<13;k++)if(u.has(P[j][k]))q|=1<<k;
  A[i*old.length+j]=q;
  D[i*old.length+j]=pc(q);
 }}
 return {W,Q,R,old,newv,nb,P,A,D};
}
function qc(W,check3=1){
 const z=qdata(W),{old,newv,nb,P,A,D}=z,K=old.length,N=newv.length,all=8191;
 let thresholdK=0,Kloss=0,j2=0,confirmed=0,good=null;
 let maxOwn=0,maxJoint=0;
 for(let i=0;i<N;i++){for(let j=0;j<K;j++)if(D[i*K+j]>maxOwn)maxOwn=D[i*K+j];}
 for(let i=0;i<N;i++)for(let j=i+1;j<N;j++){
  let support=0,p=K*i,q=K*j,M=0;
  for(let k=0;k<K;k++)if((A[p+k]|A[q+k])===all){support++;M|=1<<k;}
  if(support>maxJoint)maxJoint=support;
  if(support<3)continue;j2++;
  const ii=[];for(let k=0;k<K;k++)if(M>>k&1)ii.push(k);
  for(let a=0;a<ii.length&&!good;a++)for(let b=a+1;b<ii.length&&!good;b++)for(let c=b+1;c<ii.length&&!good;c++){
   let V=old.filter((_,k)=>k!==ii[a]&&k!==ii[b]&&k!==ii[c]).concat([newv[i],newv[j]]);
   let covered=new Set(V.flatMap(k=>nb[k]));
   if(covered.size===z.R.length){confirmed++;good={W,V,new:[newv[i],newv[j]],removed:[old[ii[a]],old[ii[b]],old[ii[c]]]};}
  }
 }
 return {maxOwn,maxJoint,j2,confirmed,good};
}

function q3(W){
 const z=qdata(W),{old,newv,A,D,nb,R}=z,K=old.length,N=newv.length,all=8191, zero=old.indexOf(z.Q[0]);
 if(zero!==0)throw Error('normalize origin');
 const mask=Array(8192),over=[];
 for(let i=0;i<N;i++)if(A[i*K]){
  const a=A[i*K];for(let b=a;b;b=(b-1)&a)(mask[b]??=[]).push(i);
  over.push(i);
 }
 let exact=0,trips=0,max=0,win=null;
 for(let i of over)for(let j of over)if(j>i){
  let left=all^(A[i*K]|A[j*K]);if(!left)throw Error("pair covers");
  for(let k of mask[left]??[])if(k>j){trips++;let f=[];
   for(let x=0;x<K;x++)if((A[i*K+x]|A[j*K+x]|A[k*K+x])===all)f.push(x);
   if(f.length>max)max=f.length;
   if(f.length<4)continue;
   for(let a=1;a<f.length;a++)for(let b=a+1;b<f.length;b++)for(let c=b+1;c<f.length;c++){
    let V=old.filter((_,v)=>v!==0&&v!==f[a]&&v!==f[b]&&v!==f[c]).concat([newv[i],newv[j],newv[k]]);
    if(new Set(V.flatMap(v=>nb[v])).size===R.length){exact++;win={W,V,trips,new:[newv[i],newv[j],newv[k]],rem:[old[0],old[f[a]],old[f[b]],old[f[c]]]};}
   }
  }
 }
 return {win,exact,trips,max};
}

let wCt=0,triples=0,hist={};
for(let W of subspaces)if(W.some(x=>pc(H[x])&1)){
 let r=qc(W),t=q3(W);wCt++;
 if(r.maxJoint!==0||r.good||t.win||t.max>1)throw Error('quotient obstruction');
 triples+=t.trips;hist[t.max]=(hist[t.max]||0)+1;
}
if(wCt!==496||triples!==222||hist[0]!==468||hist[1]!==28)throw Error('quotient scan counts');
console.log(JSON.stringify({oddOrderFourSubgroups:wCt,triplesCoverOrigin:triples,maxMandatoryOwnerHist:hist,claim:'no 2->1, 3->2, or 4->3 swaps'}));

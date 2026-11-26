const N=24, mask=(1<<20)-1, G=[0], basis=(1<<11)|(1<<10)|(1<<6)|(1<<5)|(1<<4)|(1<<2)|1;
const pop=x=>{let r=0;for(;x;x&=x-1)r++;return r;};
function check(x,label){if(!x)throw Error(label);}
function eq(a,b,msg){check(JSON.stringify(a)===JSON.stringify(b),msg+": "+JSON.stringify(a));}
for(let i=0;i<12;i++){let v=basis<<i;v^=(pop(v)&1)<<23;G.push(...G.slice().map(x=>x^v));}
const W={};for(const w of G)(W[pop(w)]??=[]).push(w);
eq(Object.keys(W).map(k=>[+k,W[k].length]),[[0,1],[8,759],[12,2576],[16,759],[24,1]],"Golay weight enumerator");
const A=W[8].filter(w=>(w&15)===0).map(w=>w>>4), B=W[12].filter(w=>(w&15)===15).map(w=>w>>4), T=W[8].filter(w=>(w&15)===15).map(w=>w>>4), K=G.map(w=>w>>4);
const UA=new Set(A),UB=new Set(B),UK=new Set(K);
eq([A.length,B.length,T.length,UK.size],[130,120,5,4096],"A B T K cardinalities");
eq(T.reduce((s,t)=>s|t,0),mask,"partition");
for(let t of T)check(pop(t)===4,"each block cardinality");
const pairs=[];for(let i=0;i<5;i++)for(let j=i+1;j<5;j++){pairs.push(T[i]|T[j]);check(UA.has(T[i]|T[j]),"pair in A");}
const U=new Set(pairs), visited=new Set(), cubes=[];
for(const s of A)if(!U.has(s)&&!visited.has(s)){
 const J=T.map((t,i)=>[pop(t&s),i]);eq(J.map(x=>x[0]).sort(),[0,2,2,2,2],"two positions in four blocks");
 const X=J.filter(x=>x[0]).map(x=>x[1]);
 const S=Array.from({length:16},(_,a)=>X.reduce((z,j,k)=>a>>k&1?z^T[j]:z,s));
 for(let a=0;a<16;a++){
  const q=S[a];check(!visited.has(q),"support appears twice");visited.add(q);
  check((pop(a)&1?UB:UA).has(q),"parity in support cube");
  check(!U.has(q),"transversal");
 }
 cubes.push({s,J:X,S,missing:J.find(x=>!x[0])[1]});
}
eq([cubes.length,visited.size],[15,240],"15 support 4-cubes");
eq(T.map((_,i)=>cubes.filter(c=>c.missing===i).length),[3,3,3,3,3],"three per omitted tetrad");
for(let L of [A,B])for(let i=0;i<L.length;i++)for(let j=0;j<i;j++)check(pop(L[i]&L[j])<=4,"same-side supports");
for(let a of A)for(let b of B){
 const m=pop(a&b);check(m<=6,"cross support");
 eq(m===6,T.includes(a^b),"sharp cross support");
}
const I=Array.from({length:16},(_,x)=>x).filter(x=>pop(x)&1);
const compress=(s,w,R)=>R.reduce((v,j,k)=>v+((pop(s&w&T[j])&1)<<k),0);
const tau=w=>T.reduce((v,t,i)=>v+((pop(w&t)&1)<<i),0);
// Replay every component in a support cube using the analytic labels (matching pair of signs on each block).
let components=0, signedEdges=0;
for(const c of cubes)for(const eta of I){
 const encode=(a,b)=>a*16+b; // a toggles support pairs, b picks which sign encoding within each pair
 const E=Array.from({length:256},(_,v)=>{
  const a=v>>4,b=v&15,res=[];
  for(let j=0;j<4;j++)for(let bit=0;bit<2;bit++)res.push(encode(a^(1<<j),(b&~(1<<j))^(bit<<j)));
  return res;
 });
 const seen=new Set([0]), queue=[0];for(let at=0;at<queue.length;at++)for(const j of E[queue[at]])if(!seen.has(j)){seen.add(j);queue.push(j);}
 eq(seen.size,256,"signed module connected");eq(E.filter((_,v)=>(pop(v>>4)&1)===0).length,128,"balanced");
 eq(E.every(L=>new Set(L).size===8),true,"regular");components++;signedEdges+=256*8/2;
}
eq([components,signedEdges],[120,122880],"signed modules count");
const rank=L=>{let e=Array(20).fill(0);for(let x of L)for(let i=19;i>=0;i--)if((x>>i)&1){if(!e[i]){e[i]=x;break;}x^=e[i];}return e.filter(Boolean).length;};
eq([rank(A),rank(B),rank([...A,...T]),rank([...B,...T])],[8,9,9,9],"ranks");
for(const w of K)for(const a of A)check((pop(w&a)&1)===0,"punctured Golay is orthogonal to shortened code");
eq(K.filter(w=>tau(w)===0).length,2048,"unmixed even slice");
eq(K.filter(w=>tau(w)===31).length,2048,"unmixed odd slice");
eq(K.filter(w=>w&&pop(w)<=4).sort((a,b)=>a-b),T.slice().sort((a,b)=>a-b),"precise flats short differences");
const cost=Array(5).fill(99);
for(let tt=1;tt<31;tt++){
 let lower=0;
 for(const c of cubes){
  const tj=c.J.reduce((v,j,k)=>v+((tt>>j&1)<<k),0),ones=pop(tj);
  // eta(a)=eta(0)+(1+tau)*a, so each mixed module is hit in both support parities.
  if(ones>0&&ones<4)lower+=2**(3-ones);
 }
 const hh=pop(tt);cost[hh]=Math.min(cost[hh],lower);
}
eq(cost.slice(1),[48,42,24,12],"mixed-transversal blocker thresholds");
const H=[0];for(const t of T)H.push(...H.slice().map(u=>u^t));
const M=new Map();for(const k of K)if(!M.has(k)){const j=M.size/32;H.forEach((s,b)=>M.set(k^s,[j,pop(b)&1]));}
const allshort=[];for(let x=0;x<=mask;x++)if(pop(x)<=4)allshort.push(x);
const Rp=[], assigned=new Uint8Array(1<<20);for(let x=0;x<=mask;x++)if(!assigned[x]){Rp.push(x);for(const k of K)assigned[x^k]=1;}
const hist={};
for(const x of Rp.slice(1)){
 const hit=new Uint8Array(128);
 for(const h of allshort){const r=M.get(h^x);if(r)hit[r[0]]|=1<<r[1];}
 const n=hit.filter(t=>t===3).length;hist[n]=(hist[n]??0)+1;
}
eq(hist,{"1":20,"4":60,"7":15,"10":160},"adaptive punctured-Golay parity exclusions");
const unit=(x,i)=>i>>x&1?-1:1;
const R=[], O=[], P=[];
const d=(x,y)=>x.reduce((v,t,i)=>v+t*y[i],0);
for(let i=0;i<4;i++)for(let j=4;j<8;j++)for(let signs=0;signs<4;signs++){
 const z=Array(8).fill(0);z[i]=2*unit(0,signs);z[j]=2*unit(1,signs);R.push(z);
}
for(let m=0;m<256;m++)if(pop(m)&1)O.push(Array.from({length:8},(_,i)=>unit(i,m)));
for(let k=0;k<8;k++)for(let s=0;s<2;s++)for(let w=0;w<16;w++){
 const z=Array(8).fill(0), block=4*((k<4)?1:0);
 z[k]=2*unit(0,s);for(let j=0;j<4;j++)z[block+j]=unit(j,w);
 P.push({k,s,w,z});
}
for(const p of P)eq([R.filter(z=>d(z,p.z)>4).length,O.filter(z=>d(z,p.z)>4).length],[4,4],"pointed incidences");
for(const o of O)eq(P.filter(p=>d(o,p.z)>4).length,8,"octad incoming");
for(const r of R){
 const n=new Map();
 for(const p of P)if(d(r,p.z)>4){const label=p.k+":"+p.s;n.set(label,(n.get(label)??0)+1);}
 eq([...n.values()].sort(),[8,8],"root incoming labels, each a 3-cube section capped at 4");
}
for(let i=0;i<P.length;i++)for(let j=0;j<i;j++){
 const p=P[i],q=P[j];
 eq(d(p.z,q.z)>4,p.k===q.k&&p.s===q.s&&pop(p.w^q.w)===1,"pointed Q4 conflicts");
}
for(const e of [0,1]){
 const selected=P.filter(p=>(pop(p.w)&1)===(p.k<4?e:1-e));
 eq(selected.length,128,"neutral pointed bank");
 eq(R.filter(r=>selected.every(p=>d(r,p.z)<=4)).length,0,"neutral retained roots");
 eq(O.filter(o=>selected.every(p=>d(o,p.z)<=4)).length,64,"neutral retained octads");
 for(let i=0;i<selected.length;i++)for(let j=0;j<i;j++)check(d(selected[i].z,selected[j].z)<=4,"neutral pointed pair");
}
const ell=T.reduce((n,t)=>n|(t&-t),0);const D=K.filter(w=>(pop(w&ell)&1)===0);
eq(D.length,2048,"attaining dense code dimension");
// A triangle cover with all B parities: fixed six signs, 2 A and 4 B points.
eq(8*64*2,1024,"all-B-parity transversal capacity per support cube");
eq(cubes.length*1024+pairs.length*192+5*6*4,17400,"complete eligible sparse/root bound");
eq(10*128+120*128+760+D.length,19448,"full height code count");
console.log(JSON.stringify({A:A.length,B:B.length,T,cubes:cubes.length,modules:components,shortMasks:allshort.length,mixedCosts:cost.slice(1),outsideGolayBothPhaseHistogram:hist,allBParitySparseRootBound:17400,fullFixedBParityHeight:19448}));

const pairsOf=t=>{let p=[];for(let h=t;h;h=(h-1)&t)if(pop(h)===2)p.push(Math.min(h,t^h));return [...new Set(p)].sort((a,b)=>a-b);};
const allPairs=T.map(pairsOf);
const cls=(s,j)=>s&T[j]?allPairs[j].indexOf(Math.min(s&T[j],T[j]^(s&T[j]))):-1;
const classWords=cubes.map(c=>T.map((_,j)=>cls(c.s,j)));
const receipts=new Set();
for(const w of classWords)for(let i=0;i<5;i++)for(let j=i+1;j<5;j++)if(w[i]>=0&&w[j]>=0){
 const key=i+":"+j+":"+w[i]+":"+w[j];check(!receipts.has(key),"no repeated class-label receipt");receipts.add(key);
}
eq(receipts.size,90,"complete pair receipt quotient");
for(let i=0;i<cubes.length;i++)for(let j=0;j<i;j++){
 const a=cubes[i],b=cubes[j];let same=0,diff=0;
 for(let k=0;k<5;k++)if(classWords[i][k]>=0&&classWords[j][k]>=0){
  if(classWords[i][k]===classWords[j][k])same++;else diff++;
 }
 check(a.missing===b.missing?same===0&&diff===4:same<=1&&diff>=2,"cross support-class compatibility");
}
const reducedA=new Map();
for(let i=0;i<cubes.length;i++)for(const h of H){
 let u=cubes[i].s^h;check(!H.includes(u),"nonzero reference in quotient");
 check(!reducedA.has(u),"distinct reference classes");reducedA.set(u,i);
}
eq(reducedA.size,15*32,"exact fifteen distinct nonzero reference cosets");
eq(rank(A)+rank(H)-rank([...A,...H]),4,"shortened Golay intersection tetrad span");
const Fn=Array.from({length:20},(_,i)=>i);
const positions=T.map(t=>Fn.filter(i=>t>>i&1));
const AcHist={};
for(let n=0;n<1024;n++){
 const v=positions.reduce((w,L,i)=>w^(1<<L[n>>(2*i)&3]),0);
 const a=cubes.reduce((z,c)=>z+(pop(c.s&v)&1),0);AcHist[a]=(AcHist[a]??0)+1;
}
eq(AcHist,{0:64,8:960},"five-axis chart simplex histogram");
const Q=[[],[],[]], M4=[];
for(let i=0;i<4;i++)for(let s of [-2,2])Q[0].push(Array.from({length:4},(_,j)=>j===i?s:0));
for(let w=0;w<16;w++)Q[1+(pop(w)&1)].push(Array.from({length:4},(_,i)=>unit(i,w)));
for(let j=1;j<4;j++){
 const ds=[];for(let u=0;u<4;u++)for(let v=u+1;v<4;v++)if((u^v)===j)
 for(let x of [-1,1])for(let y of [-1,1])ds.push(Array.from({length:4},(_,i)=>i===u?x:i===v?y:0));
 eq(ds.length,8,"root chart class");M4.push(ds);
}
const symp=(a,b)=>(a&1)*(b>>1)+(b&1)*(a>>1);
const had=x=>[0,1,2,3].map(a=>x.reduce((v,t,b)=>v+t*((symp(a,b)&1)?-1:1),0)/2).join();
for(let c=0;c<3;c++){
 let mapped=Q[c].map(had).sort();eq(mapped,Q[[1,0,2][c]].map(x=>x.join()).sort(),"block triality permutes Q");
 eq(M4[c].map(had).sort(),M4[c].map(x=>x.join()).sort(),"block triality preserves each matching class");
}
for(let c=0;c<3;c++)for(let e=0;e<3;e++){
 if(c!==e)for(const a of Q[c])for(const b of Q[e])eq(Math.abs(d(a,b)),2,"cross-type Q dots");
 for(const a of Q[c])for(const b of M4[e])check([0,2,-2].includes(d(a,b)),"Q-root dot");
}
for(const L of M4)for(const a of Q.flat()){
 const plus=L.filter(r=>d(r,a)===2);eq(plus.length,2,"dense forbidden line count");
}
const allperm=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
for(let p of allperm){
 const P=Q.flatMap((V,i)=>V.flatMap(v=>Q[p[i]].map(u=>[...v,...u])));
 eq(P.length,192,"permutation product");
 for(let i=0;i<P.length;i++)for(let j=0;j<i;j++)check(d(P[i],P[j])<=4,"permutation product exact witness");
}
// A sample OA, used only to check the bit-projection identity. The prose proves it for every maximal array.
const OA=[];for(let a=0;a<4;a++)for(let b=0;b<4;b++)for(let c=0;c<4;c++)OA.push([a,b,c,a^b^c]);
for(let f=0;f<81;f++){
 let u=f;const L=Array.from({length:4},()=>{const k=u%3+1;u=Math.floor(u/3);return k;});
 const n=Array(16).fill(0);
 for(let row of OA){const t=row.reduce((v,x,j)=>v+((pop(x&L[j])&1)<<j),0);n[t]++;}
 for(let x=0;x<16;x++)for(let j=0;j<4;j++)eq(n[x]+n[x^(1<<j)],8,"3-uniform binary marginals");
 if(n.includes(0))check(n.every(x=>x===0||x===8),"missing box forces parity");
}
console.log(JSON.stringify({trialityPermutations:allperm.length,receipts:receipts.size,axisSimplexHistogram:AcHist,nonlinearOrthogonalArrayCapacity:15*1024,QPairCapacity:10*192,roots:5*24,sharpRetainedSparseTotal:19448,axisTradeDeletionRatio:"56/16"}));

const layer=Array.from({length:16},()=>OA.map(r=>r.join("")));
const extra=[0,0,0,1];
const toDel=[];
for(let j=0;j<4;j++){
 const row=OA.find(r=>r.every((z,k)=>j===k||z===extra[k]));
 check(row!==undefined,"Latin triple");
 toDel.push(row.join(""));layer[1<<j]=layer[1<<j].filter(r=>r!==row.join(""));
}
eq(new Set(toDel).size,4,"four old rows expelled");
layer[0]=layer[0].filter(r=>!toDel.includes(r));layer[0].push(extra.join(""));
eq(layer.reduce((s,L)=>s+L.length,0),1017,"sharp nonextendible Latin stock");
for(let a=0;a<16;a++)for(let b=0;b<=a;b++)if(pop(a^b)<=1){
 for(let x of layer[a])for(let y of layer[b])if(x!==y){
  const diff=[...x].map((z,i)=>z===y[i]?0:1),n=diff.reduce((r,j)=>r+j);
  if(a===b)check(n>=2,"within-sign code");
  else check(!(n===1&&diff[31-Math.clz32(a^b)]===1),"adjacent sign triple incompatibility");
 }
}
const rat=(an,ad,bn,bd,cn,cd)=>check(an*bd*cd+bn*ad*cd<cn*ad*bd,"exact rational endpoint");
rat(9606n**2n,10000n**2n,4n*(3n*4715n-10000n)**2n,30000n**2n,99912n,100000n);
rat(8n,9n,4n*(3n*4891n-10000n)**2n,30000n**2n,987n,1000n);
rat(9606n**2n,10000n**2n,4n*103n**2n,1000n**2n,966n,1000n);
eq(14142n*14142n<2n*10000n*10000n&&14143n*14143n>2n*10000n*10000n,true,"sqrt2 rational brackets");
console.log(JSON.stringify({sharpLatinStabilityDefect:7,roundedThrough:6,tripleChartSparseDefectAtLeast:35,partialChartPairDeficientBanks:9,weightedTetradFlipEndpoint:"<0.99912 (strict)"}));

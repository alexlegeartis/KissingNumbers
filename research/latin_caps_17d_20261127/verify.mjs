// Standalone Node/V8 integer replay for Lower84 17D (no external dependencies).

const all=[0,1,2,3],Id=all; const gen=a=>a.length?a.flatMap((x,j)=>gen(a.filter((_,i)=>i!=j)).map(y=>[x,...y])):[[]];
const P=gen(all);
const pos=(a,m,n)=>+(a[m]==n);
const e=(A,B,u,v)=>{let [r,s]=all.filter(t=>t!=u&&t!=v);for(let m of all)for(let n of all)if(+(A[r][m]==B[r][n])+ +(A[s][m]==B[s][n])!=+(A[u][m]==B[u][n])+ +(A[v][m]==B[v][n]))return false;return true;};
const R=[[[Id,Id,Id,Id]]];for(let u=1;u<4;u++){let rr=[], [r,s]=all.filter(t=>t!=0&&t!=u);for(let z of P)for(let a of P){let b=[];for(let m of all){if(a[m]==m)b.push(z[m]);else if(a[m]==z[m])b.push(m);else break;}if(b.length<4|| new Set(b).size<4)continue;const A=[Id,null,null,null];A[u]=z;A[r]=a;A[s]=b;if(!e(R[0][0],A,0,u))throw Error('row');rr.push(A);}R.push(rr);}
const pc=new Uint8Array(65536);for(let i=1;i<65536;i++)pc[i]=pc[i>>1]+(i&1);
const rows=(rr)=>rr.flatMap((f,u)=>all.map(m=>f.reduce((mask,col,v)=>mask|(((v==u?15:0)^(1<<col[m]))<<(4*v)),0)));
const rank=rr=>{const B=new Uint16Array(16);let k=0; for(let x of rr){while(x){let i=31-Math.clz32(x); if(!B[i]){B[i]=x;k++;break;}x^=B[i];}}return k;};
const blocks=(()=>{let out=[];for(let a=0;a<16;a++)for(let b=a+1;b<16;b++)for(let c=b+1;c<16;c++)for(let d=c+1;d<16;d++)out.push(1<<a|1<<b|1<<c|1<<d);return out})();
const count={},repr={},dualStats={},cliqueExamples={};
for(let a of R[1])for(let b of R[2])if(e(a,b,1,2))for(let c of R[3])if(e(a,c,1,3)&&e(b,c,2,3)){let masks=rows([R[0][0],a,b,c]);let r=rank(masks);
 if(!count[r]){count[r]=0;repr[r]=masks;}
 count[r]++;
 if(r>=7){let bad=blocks.filter(x=>masks.every(h=>pc[x&h]%2==0));let pair=null;for(let i=0;i<bad.length&&!pair;i++)for(let j=i+1;j<bad.length&&!pair;j++)if(pc[bad[i]^bad[j]]==4)pair=[bad[i],bad[j]];let key=r+','+bad.length+','+!!pair;dualStats[key]=(dualStats[key]||0)+1;if(!cliqueExamples[key])cliqueExamples[key]={pair,masks};}
}
const result={sizes:R.map(x=>x.length),count,repr,dualStats,cliqueExamples};
const H=x=>{x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0xf0f0f0f)*0x1010101)>>>24};
const par=x=>H(x)&1, span=gs=>gs.reduce((a,g)=>a.concat(a.map(x=>x^g)),[0]);
const ebl=Array.from({length:16},(_,i)=>i);
const code=f=>ebl.reduce((a,i)=>a|(f(i)<<i),0);
const lin=[0,1,2,3].map(i=>code(x=>x>>i&1)), mon=[];for(let i=0;i<4;i++)for(let j=0;j<i;j++)mon.push(lin[i]&lin[j]);
const aff=span([65535,...lin]), reps=span(mon), rm2=reps.flatMap(p=>aff.map(a=>p^a)), oct=aff.filter(a=>H(a)==8);
const bmap=reps.map((v,i)=>aff.some(a=>H(a^v)==6)), bis=(a,b)=>par(reps[a]&reps[b]);
const s6=rm2.filter(w=>H(w)==6);

const W=[0,13,33,46,50,60,23,27],NR=W.flatMap(z=>aff.map(s=>s^reps[z]));
if(NR.length!==256||new Set(NR).size!==256)throw Error('NR count');
for(let i=0;i<256;i++)for(let j=0;j<i;j++)if(H(NR[i]^NR[j])<6)throw Error('NR min');
const freq={};const ty={};
for(let m=0;m<65536;m++)if(H(m)===6){const ct=new Uint16Array(64);for(let w of NR){let z=0,k=0;for(let i=0;i<16;i++)if(m>>i&1)z|=(w>>i&1)<<k++;
ct[z]++;}for(let n of ct)freq[n]=(freq[n]||0)+1;let k=ct.reduce((a,b)=>Math.min(a,b),256);ty[k]=(ty[k]||0)+1;}
const outNR={freq,ty};
let counts={},mx={},non7=[];
function dfs(cur,remain) {let k=cur.length;counts[k]=(counts[k]||0)+1;
  let extend=reps.flatMap((_,i)=>cur.every(j=>bmap[i^j])?[i]:[]);
  if(!extend.length){mx[k]=(mx[k]||0)+1;if(k>=7)non7.push(cur)}
  for(let j=0;j<remain.length;j++)dfs([...cur,remain[j]],remain.slice(j+1).filter(x=>bmap[x^remain[j]]))
}
dfs([0],reps.flatMap((_,i)=>bmap[i]?[i]:[]));
const bout={counts,mx,first:non7[0]};

const check=(s,m)=>{if(!s)throw Error(m);};
const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
check(eq(result.sizes,[1,53,53,53]),'complete reduced row menu');
check(eq(result.count,{6:16,7:72,8:60}),'biplane count');
check(eq(result.dualStats,{'7,28,true':72,'8,12,true':60}),'dual code and weight-four subgroup');
check(eq(outNR.ty,{2:448,4:7560})&&eq(outNR.freq,{2:14336,4:483840,6:14336}),'six-coordinate multiplicities');
check(eq(bout.counts,{1:1,2:28,3:168,4:336,5:280,6:168,7:56,8:8})&&eq(bout.mx,{4:56,8:8}),'bent quotient maximal cliques');
const sampleCross=[]; for(let u=0;u<4;u++)for(let k=0;k<4;k++){let r=0;for(let v=0;v<4;v++)r|=((u===v?15:0)^(1<<k))<<(4*v);sampleCross.push(r)}
check(rank(sampleCross)===6,'Cohn-Li cross span 6');
for(const [dim,masks] of Object.entries(result.repr)){check(rank(masks)===+dim,'literal representative rank');
 check(masks.length===16&&masks.every(h=>H(h)===6)&&masks.every((h,i)=>masks.slice(0,i).every(k=>H(k&h)===2)),'literal symmetric design');
 check(masks.reduce((k,m)=>k^m,0)===0,'even complete column degrees');
 check(rank([...masks,65535])===+dim,'all-one in incidence span'); }
const stdTypes={};
const hyperplanes=aff.filter(a=>H(a)===8);
for(let i=0;i<65536;i++)if(H(i)===6){
 let good=hyperplanes.every(h=>H(h&i)<=4); let q=rm2.includes(i);
 check(good===q,'RM2 is exact sparse support gate');
 if(good){let t=[0,1,2,3].map(v=>H(i>>(4*v)&15)).sort().join('');stdTypes[t]=(stdTypes[t]||0)+1;}
}
check(eq(stdTypes,{'1113':256,'0222':192}),'standard sparse supports');
// Integer polynomial below checked against P_j(1) and every attainable dot.
function gcd(a,b){return b?gcd(b,a%b):a<0n?-a:a;}
const F=(n,d=1n)=>{let q=gcd(n,d);return d<0n?[-n/q,-d/q]:[n/q,d/q];};
const plus=(a,b)=>F(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const mul=(a,n,d=1n)=>F(a[0]*n,a[1]*d);
const coeff=[55n,704n,3201n,8520n,13464n,11016n];
function polynomial(x){let p=[F(1n),F(x,6n)];for(let k=1n;k<5n;k++)p.push(mul(plus(mul(p[Number(k)],(2n*k+14n)*x,6n),mul(p[Number(k)-1],-k)),1n,k+14n));return p.reduce((a,v,i)=>plus(a,mul(v,coeff[i])),F(0n));}
for(let x=-6n;x<=2n;x++)check(polynomial(x)[0]<=0n,'sparse Delsarte sign');
check(eq(polynomial(6n).map(String),['36960','1']),'sparse Delsarte endpoint');
console.log(JSON.stringify({reducedDesigns:result.count,dualWitnesses:result.dualStats,NRsix:outNR.ty,quadraticMaximalCliques:bout.mx,standardSparse:stdTypes,universalSparseDelsarte:672}));

const T=[38192,22787,12437,857];
const WW=[0,1,2,3].map(v=>[0,1,2,3].filter(u=>u!==v).reduce((z,u)=>z|((1|(1<<(u^v)))<<(4*u)),0));
check(eq(T,WW)&&T.every(m=>H(m)===6&&aff.includes(m^T[0])),'identical bent checks');
const allowed=reps.flatMap((r,i)=>par(r&T[0])===0?[i]:[]);
let df=null;
function find6(c,rest){if(c.length===6){df=c;return true;}for(let j=0;j<rest.length;j++)if(find6([...c,rest[j]],rest.slice(j+1).filter(i=>bmap[i^rest[j]])))return true;return false;}
check(allowed.length===32&&find6([0],allowed.filter(i=>bmap[i])),'six explicit dense affine cosets');
const shift=allowed.find(i=>df.every(j=>!df.includes(i^j)));
check(shift!==undefined,'disjoint translated dense bank');
const points=[];
function add(mask,neg,amp,rail){
 const p=new Int8Array(17);for(let j=0;j<16;j++)if(mask>>j&1)p[j]=(neg>>j&1?-1:1)*amp;p[16]=rail;points.push(p);
}
for(let i=0;i<16;i++)for(let j=i+1;j<16;j++){let m=(1<<i)|(1<<j);for(let neg of [0,1<<i,1<<j,m])add(m,neg,6,0);}
for(let m of hyperplanes)for(let neg=m;neg;neg=(neg-1)&m)if(par(neg)===1)add(m,neg,3,0);
check(points.length===4320,'equator');
const pairbits=(u,r,b)=> {const cols=[0,1,2,3].filter(j=>!b?(j===0||j===r):(j!==0&&j!==r));return [u*4+cols[0],u*4+cols[1]];};
for(let omit=0;omit<4;omit++){
 const U=[0,1,2,3].filter(u=>u!==omit);
 for(let a=0;a<4;a++)for(let b=0;b<4;b++){const axis=[a,b,2^a^b];let m=0,core=0,op=[];
  for(let i=0;i<3;i++){const [p,q]=pairbits(U[i],U[i]^omit,axis[i]&1);op.push((1<<p)|(1<<q));m|=op[i];if(axis[i]>>1)core|=1<<q;}
  check(H(m)===6&&rm2.includes(m),'constructive sparse support');
  for(let g=0;g<8;g++){const neg=op.reduce((x,p,i)=>x^(g>>i&1?p:0),core);add(m,neg,3,3);add(m,neg,3,-3);}
 }
}
check(points.length===5344,'exact sparse additions');
for(let z of df)for(let a of aff){add(65535,a^reps[z],2,2);add(65535,a^reps[z^shift],2,-2);}
add(0,0,0,6);add(0,0,0,-6);
check(points.length===5730&&new Set(points.map(p=>p.join(','))).size===5730,'count/uniqueness');
let mxDot=-100;
for(let i=0;i<points.length;i++){const x=points[i]; let norm=2*x[16]*x[16];for(let k=0;k<16;k++)norm+=x[k]*x[k];check(norm===72,'norm');
for(let j=0;j<i;j++){ const y=points[j];let d=2*x[16]*y[16];for(let k=0;k<16;k++)d+=x[k]*y[k];if(d>36)throw Error('contact '+i+' '+j+' '+d);if(d>mxDot)mxDot=d;}}
console.log(JSON.stringify({alternateExact17D:points.length,mxDot,df,shift,pairCount:points.length*(points.length-1)/2}));

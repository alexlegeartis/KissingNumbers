'use strict';
// Zinoviev–Ericson 1154 (1999) literal 51-support equator and 144 binary cap.
// Binary and support masks transcribed from published Henry Cohn archive extract.
const assert=require('node:assert/strict');
const C=[1535,3711,3007,1055,303,4047,207,2575,567,343,2199,2151,1191,1607,903,3335,123,667,2331,2219,1579,843,3147,1419,435,3123,1235,2643,1811,739,1379,2851,2499,3715,2109,605,413,1133,685,2381,3213,1805,245,1333,3157,1685,2837,869,2469,3621,1477,2757,4089,1209,825,2265,1369,3609,489,2665,3369,1737,2953,9,2417,1649,2737,977,3473,3297,1953,3905,190,2142,798,622,3118,1358,2446,1678,4086,1142,2358,726,1430,3606,486,2726,1830,3270,2886,6,1338,2618,474,1626,3226,1258,2410,938,2762,3850,2290,882,1714,3410,2962,3682,3490,1986,380,1596,1244,2716,3356,2284,1452,2860,972,3660,2676,948,3252,2516,1876,1764,3428,3972,760,3192,2488,2904,1944,1896,3752,3528,1520,3888,48,3792,3040,1088,384,2560], Q=[531,99,387,3075,269,1045,165,2117,57,1097,2185,209,2321,1313,2593,833,1665,46,2070,198,774,90,650,2314,306,1170,2210,1570,1346,2626,156,588,3084,564,340,1124,2340,1412,2692,1304,2584,2152,1192,808,456,3120,1616,912,736,3264,3840];
const pc=x=>{let s=0;for(;x;x&=x-1)s++;return s};
assert.equal(new Set(C).size,144);assert.equal(new Set(Q).size,51);
assert(C.every((a,i)=>C.includes(a^4095)&&C.slice(i+1).every(b=>pc(a^b)>=4)));
assert(Q.every((a,i)=>pc(a)===4&&Q.slice(i+1).every(b=>pc(a&b)<=2)));
const twos=Array.from({length:4096},(_,i)=>i).filter(i=>pc(i)===2);
const threes=Array.from({length:4096},(_,i)=>i).filter(i=>pc(i)===3);
for(const p of twos){const s=Q.filter(x=>(x&p)===p).length;assert(s===4||s===5);}
for(const p of threes)assert.equal(new Set(C.map(x=>x&p)).size,8);
for(let i=0;i<12;i++){
 const m=1<<i, N=Q.filter(q=>q&m);assert.equal(N.length,17);
 const H=twos.filter(q=>!(q&m)&&!N.some(z=>(z&q)===q));
 assert.equal(H.length,4);
 const deg=Array.from({length:12},(_,j)=>H.filter(q=>q>>>j&1).length).sort((a,b)=>a-b);
 assert.deepEqual(deg,[0,0,0,0,0,0,0,0,2,2,2,2]);
}
const minimum={7:100,8:100,9:100};
for(let m=0;m<4096;m++){const t=pc(m);if(!(t in minimum))continue;
 minimum[t]=Math.min(minimum[t],Q.filter(q=>(m&q)===q).length);
}
assert.deepEqual(minimum,{7:2,8:5,9:12});
function ideals93(k,need){
 const comb=(n,r)=>{let a=[],o=[];function f(i){if(a.length===r){o.push([...a]);return;}for(;i<=n-r+a.length;i++){a.push(i);f(i+1);a.pop();}}f(0);return o};
 const Q=comb(k,4), gt=Q.map((v,i)=>Q.flatMap((u,j)=>j!==i&&u.every((a,h)=>a>=v[h])?[j]:[]));
 let I=[[]], N=[];
 for(let u=0;u<need;u++){const V=new Map();for(const S of I)for(let i=0;i<Q.length;i++)if(!S.includes(i)&&gt[i].every(j=>S.includes(j))){const Y=[...S,i].sort((a,b)=>a-b);V.set(Y.join(','),Y);}I=[...V.values()];N.push(I.length);}
 const ideals=I.map(T=>T.filter(i=>!T.some(j=>gt[j].includes(i))).map(i=>Q[i]));
 return {N,ideals};
}
function simplexNorm(n,extra) {
 const Z=extra.map(([r,b])=>{let z=0;return [r.map(x=>z+=x),b]});
 const E=Array.from({length:n},(_,i)=>i),D=Array.from({length:Z.length},(_,i)=>i);
 const choose=(A,k,f)=>{const B=[];function g(i){if(B.length===k){f(B);return;}for(;i<=A.length-k+B.length;i++){B.push(A[i]);g(i+1);B.pop();}}g(0);};
 const solve=(W,m)=>{
  let den=1n;
  for(let j=0;j<m-1;j++){
   let p=j;while(p<m&&W[p][j]===0n)p++;if(p===m)return null;if(p!==j)[W[p],W[j]]=[W[j],W[p]];
   const c=W[j][j];
   for(let i=j+1;i<m;i++){for(let k=j+1;k<=m;k++){const v=W[i][k]*c-W[i][j]*W[j][k];if(v%den)throw Error('Bareiss');W[i][k]=v/den;}W[i][j]=0n;}
   den=c;
  }
  const d=W[m-1][m-1];if(!d)return null;
  const N=Array(m).fill(0n);for(let i=m-1;i>=0;i--){let z=d*W[i][m];for(let j=i+1;j<m;j++)z-=W[i][j]*N[j];
   if(z%W[i][i])throw Error('Cramer');N[i]=z/W[i][i];}
  return d>0n?[d,N]:[-d,N.map(x=>-x)];
 };
 let best=[0n,1n],arg=[],eq=new Set(),checks=0;
 for(let k=1;k<=Math.min(n,D.length);k++)choose(E,k,J=>choose(D,k,I=>{
  ++checks;
  const W=I.map(i=>[...J.map(j=>BigInt(Z[i][0][j])),BigInt(Z[i][1])]),z=solve(W,k);if(!z)return;
  const [d,N]=z;if(N.some(x=>x<=0n))return;
  if(!Z.every(([r,b])=>J.reduce((sum,j,i)=>sum+BigInt(r[j])*N[i],0n)<=BigInt(b)*d))return;
  const x=E.map(i=>J.reduce((s,j,u)=>j>=i?s+N[u]:s,0n)),v=x.reduce((s,a)=>s+a*a,0n),b=d*d;
  if(v*best[1]>b*best[0]){best=[v,b];arg=x.map(t=>Number(t)/Number(d));}
  if(4n*v===3n*b)eq.add(x.map(t=>Number(t)/Number(d)).join(','));
 }));
 return {best:best.map(String),arg,eq:[...eq],checks};
}
function diffuseGen(n,Qextra) {
 const row=(list,b=1)=>[Array.from({length:n},(_,i)=>list.includes(i)?1:0),b];
 const cases=[[row([0,1,2])],[row([0,1,3])],[row([0,2,3]),row([0,1,4])]];
 const B=[[Array.from({length:n},(_,i)=>i===0?2:0),1],row([2,4,5,6]),row([0,1,n-2,n-1])];
 let out=[];
 for(const Q of Qextra)for(let c=0;c<3;c++)out.push({Q,c,...simplexNorm(n,[...B,...cases[c],...Q.map(t=>row(t))])});
 return out;
}
const eq=C=>Array.from({length:C},(_,i)=>i<3?0.5:0).join(',');
const src8=({n,Qextra})=>{
 const row=(list,b=1)=>[Array.from({length:n},(_,i)=>list.includes(i)?1:0),b];
 const cases=[[row([0,1,2])],[row([0,1,3])],[row([0,2,3]),row([0,1,4])]];
 const B=[[Array.from({length:n},(_,i)=>i===0?2:0),1],row([2,4,5,6])];
 return Qextra.flatMap(Q=>cases.map(cs=>simplexNorm(n,[...B,...cs,...Q.map(t=>row(t))])));
};
const r8=src8({n:8,Qextra:ideals93(8,5).ideals}),r9=diffuseGen(9,ideals93(9,12).ideals);
assert.equal(r8.length,15);assert.equal(r9.length,144);
for(const [n,R] of [[8,r8],[9,r9]])for(const o of R){
 const [a,b]=o.best.map(BigInt);
 assert(4n*a<=3n*b);assert(o.eq.every(s=>s===eq(n)));
}
console.log(JSON.stringify({minQuartets:minimum,systems:[r8.length,r9.length],equalities:[r8.filter(o=>o.eq.length).length,r9.filter(o=>o.eq.length).length]}));

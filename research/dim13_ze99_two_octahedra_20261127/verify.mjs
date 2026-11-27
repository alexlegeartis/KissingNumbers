#!/usr/bin/env node
import assert from "node:assert/strict";
// Equatorial support masks and cap signs transcribed from Cohn's 1154-point
// archive entry, read through Rafael Amichis Luengo's MIT data snapshot.
const S=[46,57,90,99,156,165,198,209,269,306,340,387,456,531,564,588,650,736,774,808,833,912,1045,1097,1124,1170,1192,1304,1313,1346,1412,1570,1616,1665,2070,2117,2152,2185,2210,2314,2321,2340,2584,2593,2626,2692,3075,3084,3120,3264,3840];
const A=[6,9,48,123,190,207,245,303,343,380,384,413,435,474,486,489,567,605,622,667,685,726,739,760,798,825,843,869,882,903,938,948,972,977,1055,1088,1133,1142,1191,1209,1235,1244,1258,1333,1338,1358,1369,1379,1419,1430,1452,1477,1520,1535,1579,1596,1607,1626,1649,1678,1685,1714,1737,1764,1805,1811,1830,1876,1896,1944,1953,1986,2109,2142,2151,2199,2219,2265,2284,2290,2331,2358,2381,2410,2417,2446,2469,2488,2499,2516,2560,2575,2618,2643,2665,2676,2716,2726,2737,2757,2762,2837,2851,2860,2886,2904,2953,2962,3007,3040,3118,3123,3147,3157,3192,3213,3226,3252,3270,3297,3335,3356,3369,3410,3428,3473,3490,3528,3606,3609,3621,3660,3682,3711,3715,3752,3792,3850,3888,3905,3972,4047,4086,4089];
const p=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
assert.equal(S.length,51);
assert.equal(new Set(S).size,S.length);
assert(S.every((s,i)=>p(s)===4&&S.slice(i+1).every(t=>p(s&t)<=2)));
assert.equal(A.length,144);
assert.equal(new Set(A).size,A.length);
assert(A.every((s,i)=>A.includes(4095^s)&&A.slice(i+1).every(t=>p(s^t)>=4)));
const degree=Array(12).fill(0);
for(let i=0;i<12;i++)for(let j=i+1;j<12;j++){
 const m=(1<<i)|(1<<j), r=S.filter(s=>(s&m)===m).length;
 assert(r===4||r===5); if(r===4){degree[i]++;degree[j]++;}
}
assert(degree.every(x=>x===4));
assert(degree.every((_,i)=>S.filter(s=>(s>>i)&1).length===17));
for(let Q of S)assert(S.some(V=>!(Q&V)&&S.includes(4095^Q^V)));
for(let a=1;a<=4;a++)for(let b=a;b<=4;b++)for(let c=b;c<=4;c++){
 assert(Math.min(8,a*b)<=a+b+c);
 assert(Math.min(16,2*a*b)-(a+b+c)<=7);
}
let n9=0,stats={},H=[];
for(let s=0;s<4096;s++){
 const r=p(s),Q=S.filter(t=>(t&s)===t);
 if(r===9){assert(Q.some((t,i)=>Q.slice(i+1).some(u=>!(u&t))));n9++;}
 if(r===3&&!S.some(t=>(t&s)===s))H.push(s);
 if(r===7)assert(Q.length);
}
assert.equal(n9,220);
const G=[[[0,10],[1,9],[2,3]],[[4,8],[5,11],[6,7]]];
assert.equal(new Set(G.flat(2)).size,12);
const H1=G.flatMap(([a,b,c])=>a.flatMap(i=>b.flatMap(j=>c.map(k=>(1<<i)|(1<<j)|(1<<k))))).sort((a,b)=>a-b);
assert.deepEqual(H,H1);
let populations={},hist={},match={};
for(let s=1;s<4096;s++){
 const m=Math.max(...S.map(t=>p(t&s)));
 if(3*m*m<=4*p(s))populations[p(s)]=(populations[p(s)]||0)+1;
}
for(let s=0;s<4096;s++){
 let d=12,E=[];
 for(let u of A){let r=p(u^s);if(r<d){d=r;E=[];}if(r===d)E.push(u^s);}
 hist[d]=(hist[d]||0)+1;
 if(d!==3)continue;
 assert.equal(E.length,15);
 const deg=Array(12).fill(0);
 for(let e of E)for(let i=0;i<12;i++)if(e>>i&1)deg[i]++;
 const type=E.map(e=>deg.flatMap((a,i)=>e>>i&1?[a]:[]).sort((a,b)=>a-b).join());
 const is2=type.includes("5,5,5");
 const table=is2?{"2,2,5":12,"4,4,4":12,"4,4,5":4,"5,5,5":0}:{"2,2,5":12,"4,4,4":9,"4,4,5":5,"4,5,5":1};
 assert(type.every(t=>Object.hasOwn(table,t)));
 for(let i=0;i<12;i++)assert.equal(E.reduce((v,e,j)=>v+((e>>i&1)?table[type[j]]:0),0),24);
 match[is2?"type2":"type1"]=(match[is2?"type2":"type1"]||0)+1;
}
assert.deepEqual(hist,{0:144,1:1728,2:1904,3:320});
assert.deepEqual(match,{type1:256,type2:64});
assert.deepEqual(populations,{1:12,3:16,12:1});
console.log(JSON.stringify({status:"PASS",quartets:S.length,nineSets:n9,holes:H,match,populations,hamming:hist}));

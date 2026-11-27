// Node 20+: node verify.mjs. Cohn/MIT 1154-entry, matching archived Lower74 literal masks.
const A=[6,9,48,123,190,207,245,303,343,380,384,413,435,474,486,489,567,605,622,667,685,726,739,760,798,825,843,869,882,903,938,948,972,977,1055,1088,1133,1142,1191,1209,1235,1244,1258,1333,1338,1358,1369,1379,1419,1430,1452,1477,1520,1535,1579,1596,1607,1626,1649,1678,1685,1714,1737,1764,1805,1811,1830,1876,1896,1944,1953,1986,2109,2142,2151,2199,2219,2265,2284,2290,2331,2358,2381,2410,2417,2446,2469,2488,2499,2516,2560,2575,2618,2643,2665,2676,2716,2726,2737,2757,2762,2837,2851,2860,2886,2904,2953,2962,3007,3040,3118,3123,3147,3157,3192,3213,3226,3252,3270,3297,3335,3356,3369,3410,3428,3473,3490,3528,3606,3609,3621,3660,3682,3711,3715,3752,3792,3850,3888,3905,3972,4047,4086,4089];
const P=[[0,3],[1,2],[4,5],[6,10],[7,8],[9,11]];
const need=(b,l)=>{if(!b)throw Error(l)};
const pc=x=>x?1+pc(x&(x-1)):0;
const sign=(z,i)=>1-2*((z>>i)&1);
need(A.length===144&&new Set(A).size===144&&A.every(w=>A.includes(w^4095)),"actual antipodal archive");
const U=A.filter(w=>w<(w^4095)), order=[];
for(let i=0;i<12;i++)for(let j=i+1;j<12;j++)order.push([i,j]);
const B=U.map(w=>{const h=Math.min(pc(w),12-pc(w));if(h===2){need(P.some(([i,j])=>(pc(w)===2?w:(w^4095))===(1<<i|1<<j)),"exceptional pair");return 15;}
 need(h===6,"balanced or exceptional");const r=P.reduce((a,[i,j])=>a+((w>>i^w>>j)&1),0);
 need([0,2,4,6].includes(r),"even splits");return [1,4,7,10][r>>1];});
need(B.reduce((a,b)=>a+b,0)===480,"positive tight frame total");
const Z=U.map(w=>order.map(([i,j])=>sign(w,i)*sign(w,j)));
need(order.every((_,j)=>Z.reduce((a,x,i)=>a+x[j]*B[i],0)===0),"tight second moment");
for(const [i,j] of P){
 const y=order.map(([a,b])=>+(a===i)+ +(b===i)- +(a===j)- +(b===j));
 need(Z.every(z=>z.reduce((s,x,k)=>s+x*y[k],0)===0),"six exact kernel vectors");
}
const mod=x=>((x%101)+101)%101,M=Z.map(z=>z.map(mod));let rk=0;
for(let j=0;j<66;j++){
 const k=M.findIndex((z,i)=>i>=rk&&z[j]);if(k<0)continue;[M[rk],M[k]]=[M[k],M[rk]];
 let c=1;while(c*M[rk][j]%101!==1)c++;for(let k=j;k<66;k++)M[rk][k]=mod(M[rk][k]*c);
 for(let i=rk+1;i<M.length;i++){let d=M[i][j];if(d)for(let k=j;k<66;k++)M[i][k]=mod(M[i][k]-d*M[rk][k]);}
 rk++;
}
need(rk===60,"rank at least 60 over Q");
let n=0;
function choose(x,i){if(x.length===5){
 const s=new Set(A.map(w=>x.reduce((a,j,k)=>a+((w>>j)&1)*(1<<k),0)));
 need(s.size===32,"full 5-coordinate projection");n++;return;
} for(;i<12;i++)choose([...x,i],i+1);}
choose([],0);need(n===792,"all support masks");
console.log(JSON.stringify({binaryCap:144,weightedTightSum:480,pairCharacterRank:rk,allFiveProjections:n,heightCeiling:"3/5",boundary:"24 signed axes",strictChargeSupportUpper:6}));

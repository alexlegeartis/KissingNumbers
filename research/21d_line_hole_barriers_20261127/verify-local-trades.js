"use strict";
const lines=[[4,5,6,8,10],[2,3,6,7,15],[1,4,9,14,15],[3,9,10,11,16],[1,6,12,13,16],[7,8,9,12,17],[4,7,11,13,18],[2,10,12,14,18],[1,3,5,17,18],[2,5,9,13,19],[6,11,14,17,19],[8,15,16,18,19],[1,2,8,11,20],[5,7,14,16,20],[10,13,15,17,20],[3,4,12,19,20],[3,8,13,14,21],[5,11,12,15,21],[2,4,16,17,21],[1,7,10,19,21],[6,9,18,20,21]];
const mask=L=>L.reduce((s,i)=>s|1<<(i-1),0);
function pop(x){let n=0;while(x){x&=x-1;n++;}return n;}
const G=lines.map(mask),xs=new Set(G),basis=[[3,4],[2,11],[2,5,6,10],[1,3,8,14]].map(mask);
const tag=x=>basis.reduce((s,b,j)=>s|((pop(x&b)&1)<<(3-j)),0);
let C=new Set([0]);for(const v of G)for(const w of [...C])C.add(v^w);
const H=new Set([...C].filter(x=>!tag(x)));
const mult=Array.from({length:16},(_,j)=>G.filter(x=>tag(x)===j).length);
if(C.size!==1024||H.size!==64||xs.size!==21||mult.join()!=='0,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1')throw Error('source');
const near=x=>G.map(g=>g^x);
const bad=(a,b)=>xs.has(a^b);
const out=[];
for(const t of [1,2,3]){
 const center=G.find(x=>tag(x)===t), R0=near(center).filter(x=>H.has(x));
 const psets=R0.map(h=>near(h).filter(x=>mult[tag(x)]===1));
 if(R0.length!==3||psets.some(v=>v.length!==12))throw Error('private sets');
 const flat=psets.flat();
 if(new Set(flat).size!==36)throw Error('private overlap');
 const opts=[];
 for(const y of C)if(mult[tag(y)]===1&&!bad(y,center)){
   let c=0n,counts=[0,0,0];
   for(let i=0;i<flat.length;i++)if(y===flat[i]||bad(y,flat[i])){c|=1n<<BigInt(i);counts[Math.floor(i/12)]++;}
   if(counts.every(q=>q===2))opts.push([y,c]);
 }

 if(opts.length!==0)throw Error('nonempty necessary candidate set');
 out.push({tag:t,privateNeighbors:flat.length,possibleSingletons:opts.length});
}
console.log({exact:out});

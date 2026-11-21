"use strict";
const fs = require("node:fs");
const path = require("node:path");
const DIM=24;
function pop(x){let n=0;for(;x;x&=x-1)n++;return n;}
function verify(list,mask){
 const g=[7,1,4,4,2,7,3,3,4,7,2,1];
 const rows=Array.from({length:12},(_,i)=>g.reduce((s,x,j)=>s|((x&1)<<(i+j)),1<<23));
 const words=new Uint32Array(4096);
 for(let i=1;i<4096;i++)words[i]=words[i&(i-1)]^rows[31-Math.clz32(i&-i)];
 if(new Set(words).size!==4096)throw Error("binary parent rank");
 const ws=Array.from(words,pop);
 const distro={};for(const w of ws)distro[w]=(distro[w]??0)+1;
 if(JSON.stringify(distro)!=='{"0":1,"8":759,"12":2576,"16":759,"24":1}')throw Error("Golay weight distribution");
 const H=list.map(id=>[words[id&4095],words[id>>>12]]);
 if(new Set(list).size!==list.length||H.some(([a,b])=>pop(a|b)!==12||((a|b)&mask)))throw Error("wrong outer support or duplicates");
 let min=Infinity,counts={};for(let i=0;i<H.length;i++)for(let j=0;j<i;j++){
  const [a,b]=H[i],[c,d]=H[j],d0=pop((a^c)|(b^d));
  if(d0<12)throw Error("outer distance");
  min=Math.min(min,d0);counts[d0]=(counts[d0]??0)+1;
 }
 return {patterns:H.length,min,counts,points:H.length*2**20};
}
function checkLocal(){
 const g=Array.from("210100100101",Number);
 const coc=Array.from({length:12},(_,i)=>g.reduce((x,a,j)=>x|(a<<(2*((i+j)%12))),0));
 const dwords=new Uint32Array(4096);
 for(let i=1;i<4096;i++)dwords[i]=dwords[i&(i-1)]^coc[31-Math.clz32(i&-i)];
 if(new Set(dwords).size!==4096)throw Error("dodeca rank");
 const distro={};
 for(const x of dwords){
  let n=0;for(let j=0;j<12;j++)if(x>>(2*j)&3)n++;
  distro[n]=(distro[n]??0)+1;
 }
 if(JSON.stringify(distro)!=='{"0":1,"6":396,"8":1485,"10":1980,"12":234}')throw Error("dodeca distance");
 const signs=[];
 for(let a=0;a<4096;a++){
  let syn=0;for(let j=0;j<12;j++)if(a>>j&1)syn^=j+1;
  if(!syn)signs.push(a);
 }
 if(signs.length!==256||signs.some(x=>x&&pop(x)<3))throw Error("sign code");
 const matches=[[[0,1],[2,3]],[[0,2],[1,3]],[[0,3],[1,2]]];
 const named=matches.map(q=>q.flatMap(([a,b])=>[-1,1].flatMap(s=>[-1,1].map(t=>Array.from({length:4},(_,i)=>i===a?2*s:i===b?2*t:0)))));
 const norm=(x,y)=>x.reduce((s,v,i)=>s+v*y[i],0);
 for(let a=0;a<3;a++)for(let b=0;b<3;b++)for(const x of named[a])for(const y of named[b]){
  if(norm(x,x)!==8||Math.abs(norm(x,y))>(a===b?8:4))throw Error("frame geometry");
 }
 const incumb=Array.from({length:16},(_,j)=>Array.from({length:4},(_,i)=>j>>i&1?1:-1));
 for(let i=0;i<4;i++)for(const v of [-2,2])incumb.push(Array.from({length:4},(_,j)=>i===j?v:0));
 if(named.flat().some(x=>incumb.some(y=>norm(x,y)>4)))throw Error("dense coupling");
 return {dodeca:dwords.length,signs:signs.length,inner:dwords.length*signs.length};
}
const lists = JSON.parse(fs.readFileSync(path.join(__dirname, "patterns.json"), "utf8"));
const expected = [[15,309,324009984],[3,938,983564288],[0,2583,2708471808]];
const local = checkLocal();
if (local.inner !== 2**20 || lists.length !== expected.length) throw Error("inner count or list count");
console.log(local);
for (let i = 0; i < expected.length; i++) {
 const [mask, n, points] = expected[i];
 const {zero, encoding} = lists[i];
 const tokens = encoding.split(",");
 if (zero !== mask || tokens.some(x => !/^[0-9a-z]+$/.test(x))) throw Error("pattern-list encoding");
 const ids = tokens.map(x => parseInt(x, 36));
 if (ids.some(x => !Number.isSafeInteger(x) || x < 0 || x >= 2**24)) throw Error("pattern ID range");
 const result = verify(ids, mask);
 if (result.patterns !== n || result.points !== points || result.min !== 12) throw Error("wrong certified counts");
 console.log(result);
}
const incumbent = {"92":29797911144,"93":42951379608,"94":73016281804,"95":109523504196,"96":206160530240};
const added = {"92":324009984,"93":983564288,"94":983564288,"95":983564288,"96":2708471808};
const totals = {"92":30121921128,"93":43934943896,"94":73999846092,"95":110507068484,"96":208869002048};
for (const d of Object.keys(incumbent)) {
 if (incumbent[d] + added[d] !== totals[d]) throw Error("incorrect combined count");
 console.log({dimension: Number(d), incumbent: incumbent[d], added: added[d], lowerBound: totals[d]});
}

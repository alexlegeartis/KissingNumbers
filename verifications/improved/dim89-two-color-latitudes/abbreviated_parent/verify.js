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
const data = JSON.parse(fs.readFileSync(path.join(__dirname, "patterns.json"), "utf8"));
const counts = [629,328,989], masks = [3,15,3], pivots = [3,-1,-1];
if (data.length !== 3 || checkLocal().inner !== 1048576) throw Error("small codes");
const g = [7,1,4,4,2,7,3,3,4,7,2,1];
const rows = Array.from({length:12},(_,i)=>g.reduce((s,x,j)=>s|((x&1)<<(i+j)),1<<23));
const linear = rows => rows.reduce((a,r)=>a.concat(a.map(w=>w^r)),[0]);
const golay = linear(rows);
const basis = Array.from({length:12},(_,i)=>[...("210100100101")].reduce((v,x,j)=>v|(Number(x)<<(2*((j+i)%12))),0));
const shortAxes = linear(basis).filter(a=>a>>>22===0);
if (shortAxes.length!==1024 || new Set(shortAxes).size!==1024) throw Error("short dodecacode size");
for (const a of shortAxes) if (a && pop((a|(a>>>1))&0x555555)<6) throw Error("short dodecacode distance");
const seeds = ["001122","002211","010212","020121","012021","011332","012313","013123","013231","021133","311111","200000","022222"];
const gray = [0,1,3,2];
const S=new Set();
for (const s of seeds) for (const m of [1,3]) for (let k=0;k<6;k++) {
 let v=0;
 for (let i=0;i<6;i++) v|=gray[(m*Number(s[(k+i)%6]))%4]<<(2*i);
 S.add(v>>>1);
}
const signs=[...S];
if (signs.length!==144 || signs.some(x=>x<0||x>=2048)) throw Error("Julin size");
for (let i=0;i<signs.length;i++) for (let j=0;j<i;j++) if (pop(signs[i]^signs[j])<3) throw Error("Julin distance");
const out=[];
for (let i=0;i<data.length;i++) {
 const {mask,pivot,ids}=data[i];
 if (mask!==masks[i]||pivot!==pivots[i]||ids.length!==counts[i]||ids.some(x=>!Number.isInteger(x)||x<0||x>=2**24)) throw Error("bad input");
 verify(ids,mask);
 if (pivot<0) {out.push(ids.length*2**20);continue;}
 const p=1<<pivot;
 const words=ids.map(x=>[golay[x%4096],golay[Math.floor(x/4096)]]);
 if (words.some(([a,b])=>!((a|b)&p))) throw Error("inactive pivot");
 for (let k=0;k<words.length;k++) for (let l=0;l<k;l++) {
  const [a,b]=words[k],[c,d]=words[l];
  if (pop(((a^c)|(b^d))&~p)<11) throw Error("punctured distance");
 }
 out.push(ids.length*shortAxes.length*signs.length);
}
const ten=JSON.parse(fs.readFileSync(path.join(__dirname,"weight10.json"),"utf8"));
const original11=data[0].ids.map(x=>[golay[x%4096]&~8,golay[Math.floor(x/4096)]&~8]);
const qdist=([a,b],[c,d])=>pop((a^c)|(b^d));
if(ten.length!==11||ten.some(r=>r.length!==2||r.some(x=>!Number.isInteger(x)||x<0||x>=2**24)||pop(r[0]|r[1])!==10||((r[0]|r[1])&3)))throw Error("weight10 pattern");
for(let i=0;i<ten.length;i++){
 for(let j=0;j<i;j++)if(qdist(ten[i],ten[j])<10)throw Error("distance within the 10-row layer");
 for(const old of original11)if(qdist(ten[i],old)<11)throw Error("distance between 10- and 11-row layers");
}
const axes10=shortAxes.filter(a=>a>>>20===0);
const signs10=signs.filter(a=>a>>>10===0);
if(axes10.length!==256||signs10.length!==72)throw Error("10-row code sizes");
if(4*40**2>80*88)throw Error("cross-radius cap");
out[0]+=ten.length*axes10.length*signs10.length;
const prev92={"92":29797911144,"93":42951379608,"94":73016281804,"95":109523504196};
const claims=JSON.parse(fs.readFileSync(path.join(__dirname,"claims.json"),"utf8"));
for (const {dimension:n,incumbent,added_over_incumbent:gain,lower_bound:value} of claims) {
 const original=n===88?12*8**10:n===89?13*8**10+2:prev92[n];
 const layer=out[n<92?0:n===92?1:2];
 if (original+layer!==value || incumbent+gain!==value) throw Error("counts");
 console.log(n,value,"+",gain,"over predecessor");
}

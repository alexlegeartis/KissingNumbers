"use strict";
const fs=require("node:fs"),path=require("node:path");
const p=[7,1,4,4,2,7,3,3,4,7,2,1],g=Array.from({length:12},(_,i)=>p.reduce((v,a,j)=>v|((a&1)<<(i+j)),1<<23));
const wt=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
const G=[0];for(const x of g)for(let i=0,k=G.length;i<k;i++)G.push(G[i]^x);
const tally={};for(const x of G)tally[wt(x)]=(tally[wt(x)]??0)+1;
if(JSON.stringify(tally)!=='{"0":1,"8":759,"12":2576,"16":759,"24":1}')throw Error("binary parent");
const octBy4=new Map();
for(const b of G)if(wt(b)===8)for(let t=b;t;t=(t-1)&b)if(wt(t)===4){
 const a=octBy4.get(t);if(a)a.push(b);else octBy4.set(t,[b]);
}
if(octBy4.size!==10626||[...octBy4.values()].some(x=>x.length!==5))throw Error("bad octad 4-design");
const tagsOf=(a,b)=>[a&~b,b&~a,a&b];
const code=Array.from("210100100101",Number);
const Dc=Array.from({length:12},(_,i)=>code.reduce((x,a,j)=>x|(a<<(2*((i+j)%12))),0));
const D=[0];for(const x of Dc)for(let i=0,k=D.length;i<k;i++)D.push(D[i]^x);
const qweight=(x,len)=>{let w=0;for(let i=0;i<len;i++)if(x>>(2*i)&3)w++;return w;};
if(new Set(D).size!==4096||D.some(x=>x&&qweight(x,12)<6))throw Error("bad dodecacode");
const H=[];
for(let b=0;b<4096;b++){
 let v=0;for(let i=0;i<12;i++)if(b>>i&1)v^=i+1;
 if(!v)H.push(b);
}
if(H.length!==256||H.some(x=>x&&wt(x)<3))throw Error("sign code");
const pairedIndices=[[4,5,6,7,8,9,10,11],[0,1,2,3,8,9,10,11],[0,1,2,3,4,5,6,7]];
const phases=[[0,5,10,15],[0,5,9,12],[0,5,10,15]];
const Project=[];
for(let u=0;u<3;u++){
 const a=pairedIndices[u],proj=x=>a.reduce((v,i,j)=>v|((x>>(2*i)&3)<<(2*j)),0);
 const gen=Dc.map(proj),P=new Uint32Array(16);
 for(let x of gen)for(let i=15;i>=0;i--)if(x>>i&1){if(P[i])x^=P[i];else{P[i]=x;break;}}
 const holes=Array.from({length:16},(_,i)=>i).filter(i=>!P[i]);
 if(holes.length!==4)throw Error("bad projected rank");
 const red=x=>{for(let i=15;i>=0;i--)if(P[i]&&(x>>i&1))x^=P[i];
   return holes.reduce((s,i,j)=>s|((x>>i&1)<<j),0);
 };
 const forbidden=new Set([0]);
 for(let i=0;i<8;i++)for(let a=1;a<=3;a++)forbidden.add(red(a<<(2*i)));
 for(let i=0;i<4;i++)for(let j=0;j<i;j++)if(forbidden.has(phases[u][i]^phases[u][j]))throw Error("phase distance");
 Project.push({proj,red,forbidden});
}
const full=x=>Project.reduce((s,{proj,red},i)=>s|(red(proj(x))<<(4*i)),0);
const labels=new Int32Array(4096);labels.fill(-1);labels[0]=0;
const messages=[0];
for(let bit=0;bit<24;bit++){
 const c=1<<bit,u=full(c);
 if(labels[u]>=0)continue;
 for(let i=0,k=messages.length;i<k;i++){
  const s=messages[i]^u;messages.push(s);labels[s]=labels[messages[i]]^c;
 }
}
if(labels.some(x=>x<0)||labels.some((x,i)=>full(x)!==i)||D.some(x=>full(x)!==0))throw Error("joint syndrome rank");
function layer(id,mask){
 const a=G[id&4095],b=G[id>>>12];
 if(wt(a|b)!==12||((a|b)&mask))throw Error("pattern weight or pinned zeros");
 const T=tagsOf(a,b),active=Array.from({length:24},(_,i)=>i).filter(i=>(a|b)>>i&1);
 if(wt(a)!==8||wt(b)!==8||wt(a^b)!==8)return {a,b,T,active,shift:0,mode:false};
 if(T.some(x=>wt(x)!==4))throw Error("bad triple partition");
 const ind=[],s=[];
 for(let u=0;u<3;u++){
  const x=T[(u+1)%3],y=T[(u+2)%3];
  const F=octBy4.get(x).map(o=>o^x).filter(w=>!(w&(x|y))).sort((a,b)=>a-b);
  if(F.length!==4||!F.includes(T[u])||F.some(w=>!G.includes(w|y)))throw Error("octad completions");
  const rank=F.indexOf(T[u]);ind.push(rank);s.push(phases[u][rank]);
 }
 const shift=labels[s[0]+(s[1]<<4)+(s[2]<<8)];
 const order=T.flatMap(t=>active.filter(i=>t>>i&1));
 return {a,b,T,active:order,shift,mode:true,ind};
}
function verify(ids,mask){
 if(new Set(ids).size!==ids.length)throw Error("duplicate pattern");
 const q=ids.map(id=>layer(id,mask));const counts={},rescuedGroups=[0,0,0];
 for(let i=0;i<q.length;i++)for(let j=0;j<i;j++){
  const x=q[i],y=q[j],dis=wt((x.a^y.a)|(x.b^y.b));counts[dis]=(counts[dis]??0)+1;
  if(dis>=12)continue;
  if(dis!==8||!x.mode||!y.mode)throw Error("illegal relaxed pair distance");
  const diff=x.T.map((t,u)=>t===y.T[u]?null:u).filter(u=>u!==null);
  if(diff.length!==1)throw Error("wrong shared tags");
  const u=diff[0],s=(x.a|x.b)&(y.a|y.b);
  if(wt(s)!==8||s!==(x.T[(u+1)%3]|x.T[(u+2)%3]))throw Error("wrong occupied overlap");
  if(pairedIndices[u].some(k=>x.active[k]!==y.active[k]))throw Error("coordinate alignment");
  const {proj,red,forbidden}=Project[u];
  if(forbidden.has(red(proj(x.shift^y.shift)))||x.ind[u]===y.ind[u])throw Error("shift distance");
  rescuedGroups[u]++;
 }
 const f=[[[0,1],[2,3]],[[0,2],[1,3]],[[0,3],[1,2]]];
 const ref=Array.from({length:16},(_,j)=>Array.from({length:4},(_,i)=>j>>i&1?1:-1));
 for(let i=0;i<4;i++)for(const a of [-2,2])ref.push(Array.from({length:4},(_,j)=>i===j?a:0));
 function dot(x,y){let s=0;for(let i=0;i<4;i++)s+=x[i]*y[i];return s;}
 const frames=f.map(x=>x.flatMap(([i,j])=>[-2,2].flatMap(a=>[-2,2].map(b=>Array.from({length:4},(_,k)=>k===i?a:k===j?b:0)))));
 for(let u=0;u<3;u++)for(const x of frames[u]){
  if(dot(x,x)!==8||ref.some(y=>Math.abs(dot(x,y))>4))throw Error("dense bound");
  for(let v=0;v<3;v++)if(v!==u&&frames[v].some(y=>Math.abs(dot(x,y))!==4))throw Error("cross matching");
 }
 return {mask,patterns:ids.length,points:ids.length*4096*256,counts,rescuedGroups};
}

const read=name=>JSON.parse(fs.readFileSync(path.join(__dirname,name),"utf8"));
const expected=[[15, 342], [3, 1064], [0, 2970]], raw=read("shifted_patterns.json"), results=new Map();
if(raw.length!==expected.length)throw Error("outer list count");
for(let i=0;i<raw.length;i++){
 const {mask,patterns}=raw[i], [m,n]=expected[i], tokens=patterns.split(",");
 if(mask!==m || tokens.length!==n || tokens.some(x=>!/^[0-9a-z]+$/.test(x)))throw Error("list structure");
 const ids=tokens.map(x=>parseInt(x,36));
 if(ids.some(x=>!Number.isInteger(x)||x<0||x>=2**24))throw Error("pattern ID range");
 const result=verify(ids,mask);
 if(result.points!==n*2**20)throw Error("wrong result");
 results.set(mask,result);
 console.log(result);
}
for(const c of read("claims.json")){
 const k=results.get(c.mask);
 if(!k || k.patterns!==c.patterns || c.inherited+k.points!==c.lower_bound || c.previous+c.increase!==c.lower_bound)throw Error("claimed count");
 console.log(c.dimension,c.lower_bound,c.increase);
}

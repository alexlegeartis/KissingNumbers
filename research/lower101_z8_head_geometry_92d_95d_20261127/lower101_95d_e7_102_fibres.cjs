'use strict';
const assert=require('node:assert/strict');
function e7fibres(){
 const A=[0,1,2,3],pairs=[[0,1],[0,2],[0,3],[2,3],[1,3],[1,2]],id=[0,1,2,4,5,6];
 const H=[],put=(c,p,h)=>H.push({c,p,h});
 for(let i=0;i<6;i++)for(let j=0;j<6;j++)put(0,[id[i],id[j]],[...A.map(a=>pairs[i].includes(a)?1:-1),...A.map(a=>pairs[j].includes(a)?1:-1)]);
 put(0,[3,3],[1,1,1,1,-1,-1,-1,-1]);put(0,[7,7],[-1,-1,-1,-1,1,1,1,1]);
 for(let i=0;i<4;i++)for(let j=0;j<4;j++){
  const s=[...A.map(a=>a===i?1:-1),...A.map(a=>a===j?-1:1)],t=[...A.map(a=>a===i?2:0),...A.map(a=>a===j?-2:0)];
  put(1,[i,j+4],s);put(1,[i+4,j],s.map(v=>-v));put(2,[i,j+4],t);put(2,[i+4,j],t.map(v=>-v));
 }
 const dot=(x,y)=>x.reduce((s,a,i)=>s+a*y[i],0),pd=(u,v)=>u===v?4:(u^4)===v?-4:0;
 const ok=(a,b)=>{const c=a.c===b.c,d=pd(a.p[0],b.p[0])+pd(a.p[1],b.p[1]),g=dot(a.h,b.h);
  if(c&&a.p.join(',')===b.p.join(','))return false;
  return g<=(c?d:4);
 };
 assert.equal(H.length,102);assert.equal(new Set(H.map(x=>x.h.join(','))).size,102);
 assert.deepEqual([0,1,2].map(c=>H.filter(a=>a.c===c).length),[38,32,32]);
 assert(H.every((a,i)=>a.h.reduce((s,v)=>s+v,0)===0&&dot(a.h,a.h)===8&&H.slice(i+1).every(b=>ok(a,b))));
 const omitted=[];
 for(let r=0;r<2;r++)for(let i=0;i<4;i++)for(let j=0;j<4;j++)if(i!==j){
  const t=Array(8).fill(0);t[r*4+i]=2;t[r*4+j]=-2;omitted.push(t);
 }
 assert.equal(omitted.length,24);
 let min=999;
 for(const h of omitted){
  let d=999;
  for(let c=0;c<3;c++)for(let x=0;x<8;x++)for(let y=0;y<8;y++){
   const bad=H.filter(a=>!ok({c,p:[x,y],h},a));
   if(bad.length<d)d=bad.length;
  }
  min=Math.min(min,d);
 }
 assert.equal(min,6);
 return {fibres:H.length,dense:H.length*8**10,withGolay:H.length*8**10+145004544,leastOneOmittedRemoval:min};
}
console.log(JSON.stringify(e7fibres()));
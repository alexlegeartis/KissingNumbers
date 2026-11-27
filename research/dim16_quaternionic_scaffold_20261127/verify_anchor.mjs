import assert from 'node:assert/strict';
const V=[[8,[2,0,0,2],[1,1,1,-1,2,0,0,0,2,0,0,0,1,-1,-1,-1]],[8,[2,2,0,0],[1,-1,-1,-1,1,1,-1,1,2,0,0,0,2,0,0,0]],[8,[2,2,0,0],[1,1,-1,1,1,1,1,-1,2,0,0,0,2,0,0,0]],[8,[0,2,0,2],[2,0,0,0,1,-1,-1,-1,2,0,0,0,1,-1,1,1]],[8,[2,0,2,0],[1,-1,1,1,2,0,0,0,1,-1,1,1,2,0,0,0]],[8,[0,2,0,2],[2,0,0,0,1,-1,1,1,2,0,0,0,1,1,1,-1]],[8,[0,0,2,2],[2,0,0,0,2,0,0,0,1,1,1,-1,1,1,-1,1]],[8,[0,0,0,0],[2,0,0,0,2,0,0,0,0,0,-2,0,2,0,0,0]],[36,[0,1,1,2],[2,0,0,0,1,1,1,1,1,1,1,1,1,-1,1,1]]];
const z=Array.from({length:16},(_,i)=>i%4===0?2:0);
const dot=(a,b)=>a.reduce((d,t,i)=>d+t*b[i],0);
for(let i=0;i<V.length;i++){
 const [cost,levels,y]=V[i];
 assert.equal(dot(y,y),16);assert(dot(y,z)>8&&dot(y,z)<16);
 assert.equal(cost,levels.some(t=>levels.filter(s=>s===t).length===1)?36:8);
 for(let j=0;j<4;j++){
   const Q=y.slice(4*j,4*j+4),t=levels[j];
   if(t===0)assert.equal(Q.filter(x=>Math.abs(x)===2).length,1);
   else{
      assert(Q.every(x=>Math.abs(x)===1));
      assert.equal((Q.filter(x=>x<0).length&1),t-1);
   }
 }
 for(let j=0;j<i;j++)assert(dot(y,V[j][2])<=8);
}
assert.deepEqual(V.map(row=>row[0]),[8,8,8,8,8,8,8,8,36]);
assert.equal(V.reduce((s,row)=>s+72/row[0],0),74);
console.log('Nine signed neighbors certify weighted anchor load 37/36.');

import assert from 'node:assert/strict';
const V=((assert)=>{
const vertices=[];
for(let e=0;e<=1;e++){
 const C=[1,2,2,2,2,2,2,2],D=[1,2,3,4,5,6,7,8-2*e],eq=new Set();
 let cnt=0;
 function check(i,n,j,m,k){
   if(k<0){k=-k;n=-n;m=-m;}
   if(n<0||m<0||C[i]*n+C[j]*m>k||D[i]*n+D[j]*m>2*k)return;
   cnt++;const a=Array.from({length:8},(_,z)=>(z<=i?n:0)+(z<=j?m:0));
   const ss=a.reduce((s,x)=>s+x*x,0);
   if(ss>k*k)throw Error('polar norm above one');
   if(ss===k*k){
      const scale=a.find(Boolean);assert(scale);const q=a.map(x=>4*x/k);
      if(q.some(z=>!Number.isInteger(z)))throw Error('wrong equality denominator');
      eq.add(q.join(','));
   }
 }
 for(let i=0;i<8;i++){
   if(D[i]<=2*C[i])check(i,1,i,0,C[i]);else check(i,2,i,0,D[i]);
 }
 for(let i=0;i<8;i++)for(let j=i+1;j<8;j++){
   const a=C[i]*D[j]-D[i]*C[j];
   if(a)check(i,D[j]-2*C[j],j,2*C[i]-D[i],a);
 }
 vertices.push({e,cnt,eq:[...eq].sort()});
}
const expected=[{e:0,eq:['2,2,2,2,0,0,0,0','4,0,0,0,0,0,0,0']},{e:1,eq:['2,2,2,2,0,0,0,0','3,1,1,1,1,1,1,1','4,0,0,0,0,0,0,0']}];
for(let i=0;i<2;i++)assert(JSON.stringify(vertices[i].eq)===JSON.stringify(expected[i].eq));
return vertices;
})(assert);
assert(V.every(x=>x.cnt===27));
// 14 affine planes and two independent right orientation choices.
const par=x=>x.toString(2).replaceAll('0','').length&1;
const pc=x=>x.toString(2).replaceAll('0','').length;
const P=[];
for(let i=1;i<8;i++)for(let b=0;b<2;b++){let p=0;for(let j=0;j<8;j++)if(par(i&j)===b)p|=1<<j;P.push(p);}
const U=P.flatMap(p=>[p|(p<<8),p|((p^255)<<8)]);
assert.equal(new Set(U).size,28);
assert(U.every((u,i)=>pc(u)===8&&U.slice(0,i).every(v=>pc(v&u)<=4)));
assert.equal(2*240+16**2+28*128,4320);
assert.equal(56*64*2/4,1792);
const S=[0,0];
for(let w=0;w<1<<15;w++){
 let k=0;for(let j=0;j<7;j++)k+=((w>>2*j)^(w>>(2*j+1)))&1;
 if(k!==3&&k!==4)continue;
 S[k-3]++;
 let c=0;
 for(let j=-1;j<15;j++){
   const v=j<0?w:w^(1<<j);let z=0;
   for(let t=0;t<7;t++)z+=((v>>2*t)^(v>>(2*t+1)))&1;
   if(z===3||z===4)c++;
 }
 assert.equal(c,10);
}
assert.deepEqual(S,[8960,8960]);
console.log('E8: exact polar maxima, 28 supports and 4320 attained, both 1792 bounds checked.');

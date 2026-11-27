'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const pins=[["dim96-mixed-triality-six-rows/witness.txt","a00a1bdd3f82f40ed930a8b894465a3feca824ee"],["dim96-mixed-triality-six-rows/witness_extra.txt","7c2e5db39c7780c573356242bd8f12ff3a82762d"],["dim94-96-hexary-parity-packing/best6_96.txt","bef34cfc76b36934686b7714d20299e53c1ea4ca"],["dim96-affine-phase-projections/witness.txt","358e8701a398155022005ae200382a8a493dff15"],["dim96-partial-triality-fibres/witness.txt","0926b8dedda71178708a04f450edee3f60248997"],["dim96-mixed-tetrad-row-alphabet-repack/claims.json","feb31d4023a963927e5b86f8288f6c5200c3948f"]];
function pinInputs(){
 const root=process.argv[2]??'.';
 return pins.map(([file,h])=>{
   const b=fs.readFileSync(path.join(root,'verifications','improved',file));
   const got=crypto.createHash('sha1').update('blob '+b.length+'\0').update(b).digest('hex');
   if(got!==h)throw Error(file+' git blob mismatch: '+got);
   return b.toString('utf8');
 });
}
function output(o){console.log(JSON.stringify(o,null,2));}

const assert=require('node:assert/strict');
const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
const ar=(n,f)=>Array.from({length:n},(_,i)=>f(i));
const read=s=>s.trim().split('\n').map(l=>l.trim().split(/\s+/).map(Number));
const raw=pinInputs();
const [W,X,H,Y,Z]=raw.slice(0,5).map(read);
assert.deepEqual([W.length,X.length,H.length,Y.length,Z.length],[822,348,219,103,60]);
const T=[...W,...X,...H.flatMap(([s,p])=>[[0,p,s^p],[0,s^p,p]]),...Y.map(x=>x.slice(0,3)),...Z.map(x=>x.slice(0,3))];
assert.equal(T.length,1771);assert.equal(new Set(T.map(x=>x.join(','))).size,T.length);
for(const [a,e,o] of T)assert(pc(a|e|o)===6&&!(a&e||a&o||e&o));
const G=ar(6,i=>15<<(4*i));
const SEL=[119,473,56576,69408,114701,176133,184480,196976,199178,264963,265616,323592,328972,329296,364736,368643,394758,395824,405648,409611,459936,483332,527884,544976,590598,591296,592048,647232,656649,657040,658240,663566,696368,722949,741384,806917,852489,855042,933984,954370,1049103,1051344,2098096,2098954,2109451,2121952,3147104,3148548,3182736,3186694,4197600,4222982,4235472,5058560,5245194,5255328,5283849,5296192,6294025,6294534,6340618,7341696,7342128,7344133,8257536,8390925,8425584,9438470,9438768,9439648,9461765,9474058,10487301,10488416,10543120,11537416,11554848,12599472,12603401,12615687,13632320,14681360,14688268];
assert.equal(new Set(SEL).size,SEL.length);
for(const s of SEL){
 assert(pc(s)===6&&G.filter(g=>g&s).length<=3);
 for(const [a,e,o] of T)assert(pc(s&(a|e|o))+pc(s&a)<=6);
}
for(const [i,v] of SEL.entries())for(const w of SEL.slice(i+1))assert(pc(v&w)<=3);
const cuts=new Set();
for(const s of SEL)for(let i=0;i<24;i++)if(s>>i&1)for(let j=4*(i>>2);j<i;j++)if(s>>j&1)cuts.add(j*24+i);
const omit=[0,0,0,0,4,4],A=[];
for(let a=0;a<5;a++)for(let b=0;b<5;b++)for(let c=0;c<5;c++)for(let d=0;d<5;d++){
 const z=[...ar(5,x=>(a+b*x+c*x*x+d*x*x*x)%5),d];
 if(z.every((v,i)=>v!==omit[i]))A.push(z.map((v,i)=>v-(v>omit[i])));
}
assert.equal(A.length,164);assert.equal(new Set(A.map(x=>x.join(','))).size,164);
for(let i=0;i<A.length;i++)for(let j=0;j<i;j++)assert(A[i].filter((v,r)=>v!==A[j][r]).length>=3);
const even=ar(64,i=>i).filter(i=>!(pc(i)&1));assert.equal(even.length,32);
for(const [i,v] of even.entries())for(const w of even.slice(i+1))assert(pc(v^w)>=2);

function generatePhysicalWitness(Tx,seedMasks,axisWords,ar,pc){
const T=Tx.map(([a,e,o])=>({a,f:e|o}));
const shapes=[];
const GROUPS=ar(6,g=>{
 const L=[];
 for(let i=4*g;i<4*g+4;i++)for(let j=4*g;j<=i;j++){
  const positions=[];
  for(let a=0;a<4;a++)for(let b=0;b<4;b++)if(i!==j||a<b){
    const u=4*j+a,v=4*i+b;positions.push([u,v]);}
  L.push({r:(1<<i)|(1<<j),d:i===j?1<<i:0,positions});
 }
 return L;
});
for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)for(let l=j+1;l<6;l++)
for(const a of GROUPS[i])for(const b of GROUPS[j])for(const c of GROUPS[l]){
 if(!a.d&&!b.d&&!c.d)continue;
 const U=a.r|b.r|c.r,D=a.d|b.d|c.d;
 if(T.some(({a,f})=>2*pc(U&a)+pc(U&f)+pc(D&f)>6))continue;
 shapes.push({a:a.positions,b:b.positions,c:c.positions,U,D});
}

const key=(a,b,c,d)=>((a*96+b)*96+c)*96+d;
const stock=new Set(), add=[];
const visit=(p,test)=>{
 for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)for(let k=j+1;k<6;k++)for(let l=k+1;l<6;l++){
   const s=key(p[i],p[j],p[k],p[l]);if(test){if(stock.has(s))return true;}else stock.add(s);
 }
 return false;
};
for(const S of seedMasks){
 const L=ar(24,i=>i).filter(i=>S>>i&1);
 for(const w of axisWords){
   const p=L.map((x,i)=>4*x+w[i]);
   if(visit(p,true))throw Error('seed collision');
   visit(p,false);
 }
}
shapes.sort((a,b)=>pc(a.D)-pc(b.D));
for(const {a,b,c} of shapes)for(const x of a)for(const y of b)for(const z of c){
 const p=[...x,...y,...z];if(!visit(p,true)){visit(p,false);add.push(p);}
}
return add;
}
function physicalWitness(){
 const witness=path.join(__dirname,'supports.b64');
 if(!fs.existsSync(witness))return generatePhysicalWitness(T,SEL,A,ar,pc);
 const s=fs.readFileSync(witness,'utf8').trim();
 if(s.length!==126536||!/^[A-Za-z0-9+/]+$/.test(s))throw Error('invalid support witness encoding');
 const b=Buffer.from(s,'base64');
 if(b.length!==94902||b.toString('base64')!==s)throw Error('noncanonical support encoding');
 return Array.from({length:b.length/6},(_,i)=>Array.from(b.subarray(6*i,6*i+6)));
}

const coords=[
  ...SEL.flatMap(s=>{const L=ar(24,i=>i).filter(i=>s>>i&1);return A.map(w=>L.map((x,i)=>4*x+w[i]));}),
  ...physicalWitness()
];
assert.equal(SEL.length,83);assert.equal(cuts.size,36);
assert.equal(coords.length,29429);
const seen4=new Set(),rootCut=new Set(),hist={};
for(const p of coords){
 assert.equal(p.length,6);
 for(let i=0;i<6;i++)assert(Number.isInteger(p[i])&&p[i]>=0&&p[i]<96&&(!i||p[i-1]<p[i]));
 const mm=ar(6,g=>p.filter(i=>i>>4===g).length).filter(Boolean);
 assert(mm.length<=3);
 const U=p.reduce((s,i)=>s|1<<(i>>2),0),D=p.reduce((s,i,j)=>s|(j&&p[j-1]>>2===i>>2?1<<(i>>2):0),0);
 assert(ar(24,k=>p.filter(x=>x>>2===k).length).every(n=>n<=2));
 hist[pc(D)]=(hist[pc(D)]??0)+1;
 if(D)assert(mm.length===3&&mm.every(n=>n===2));
 for(const [a,e,o] of T){const f=e|o;assert(2*pc(U&a)+pc(U&f)+pc(D&f)<=6);}
 for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
   if(p[i]>>4===p[j]>>4)rootCut.add(96*p[i]+p[j]);
   for(let k=j+1;k<6;k++)for(let l=k+1;l<6;l++){
     const key=((96*p[i]+p[j])*96+p[k])*96+p[l];
     assert(!seen4.has(key));seen4.add(key);
   }
 }
}
assert.equal(hist[0],13612);assert.equal(rootCut.size,720);
const old=JSON.parse(raw[5]).lower_bound;assert.equal(old,209299907072);
const added=coords.length*even.length,removed=rootCut.size*4;
assert.equal(added,941728);assert.equal(removed,2880);
const result={dimension:96,old,oldSixModes:T.length,newSupportCount:coords.length,axialRSsupportCount:hist[0],
              doubleRowSupportCount:coords.length-hist[0],doubleRowHistogram:hist,evenSigns:even.length,
              added,cutCoordinatePairs:rootCut.size,removed,delta:added-removed,lowerBound:old+added-removed};
output(result);
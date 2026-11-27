'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const pins=[["dim94-95-mixed-triality/witness.txt","5c45682ae2ea39b43ab4525d30ecaa8f6df855cd"],["dim94-96-hexary-parity-packing/best6_94.txt","5cfd7cb4793ff637bca3f5b03ef2c6cd852bb3b5"],["dim94-96-post-rs-partial-fibres/witness94_95.txt","90688c853cf8d6272dca658700bdcd7249a05ae1"],["dim94-95-matching-row-append/witness.txt","517698a0364c917753553dff7fde9e70f9820396"],["dim94-95-mixed-tetrad-axis-sign-fibres/claims.json","da251dd30e6e21182f3246f39ca7b8dd13d78315"]];
function pinInputs(){
 const root=process.argv[2]??'.';
 return pins.map(([file,h])=>{
   const b=fs.readFileSync(path.join(root,'verifications','improved',file));
   const got=crypto.createHash('sha1').update('blob '+b.length+'\0').update(b).digest('hex');
   assert.equal(got,h,file+' Git blob');return b.toString('utf8');
 });
}
function output(o){console.log(JSON.stringify(o,null,2));}

const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24;};
const ar=(n,f)=>Array.from({length:n},(_,i)=>f(i));
const read=s=>s.trim().split('\n').map(l=>l.trim().split(/\s+/).map(Number));
const raw=pinInputs();
const [W,X,Y,Z]=raw.slice(0,4).map(read);
assert.deepEqual([W.length,X.length,Y.length,Z.length],[473,152,148,85]);
const T=[...W,...X.flatMap(([s,p])=>[[0,p<<2,(s^p)<<2],[0,(s^p)<<2,p<<2]]),...Y.map(v=>v.slice(0,3))];
assert.equal(T.length,925);assert.equal(new Set(T.map(x=>x.join(','))).size,925);
for(const [a,e,o] of T)assert((a|e|o)<(1<<24)&&!((a|e|o)&3)&&pc(a|e|o)===6&&!(a&e||a&o||e&o));
for(const [s,t] of Z)assert(!((s)&3)&&[4,6,7].includes(pc(s))&&t>=0&&t<3**pc(s));
const seed=[119812,529856,1316884,1573076,1836044,1837376,4201496,4276672,12583116,12632088,12670976];
const gen=()=>{

const group=c=>c<14?0:1+Math.floor((c-14)/16), row=c=>2+(c>>2);
const admissible=(U,D)=>T.every(([a,e,o])=>{
 const F=e|o;return 2*pc(U&a)+pc((U|D)&F)+pc(D&F)<=6;
})&&Z.every(([s])=>(pc(U&s)+pc(D&s))**2<=3*pc(s));
const omit=[0,0,0,0,4,4],A=[];
for(let a=0;a<5;a++)for(let b=0;b<5;b++)for(let c=0;c<5;c++)for(let d=0;d<5;d++){
 const z=[...ar(5,x=>(a+b*x+c*x*x+d*x*x*x)%5),d];if(z.every((v,i)=>v!==omit[i]))A.push(z.map((v,i)=>v-(v>omit[i])));}
const G=ar(6,i=>ar(22,r=>r+2).filter(r=>group(4*(r-2))===i&&group(4*(r-2)+3)===i));
const banks=ar(6,i=>{
 const candidates=[];
 for(const [k,r] of G[i].entries())for(const s of G[i].slice(0,k+1)){
   const coords=[];
   for(let a=0;a<4;a++)for(let b=0;b<4;b++)if(r>s||a<b)coords.push([4*(s-2)+a,4*(r-2)+b]);
   candidates.push({U:(1<<r)|(1<<s),D:r===s?1<<r:0,coords});
 }
 return candidates;
});
const shapes=[];
for(let g=0;g<6;g++)for(let h=g+1;h<6;h++)for(let j=h+1;j<6;j++)
 for(const a of banks[g])for(const b of banks[h])for(const c of banks[j]){
 const U=a.U|b.U|c.U,D=a.D|b.D|c.D;
 if(D&&admissible(U,D))shapes.push({a:a.coords,b:b.coords,c:c.coords,D});
}
const seen4=new Set();
const occupiesFour=(p,check)=>{
 for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)for(let k=j+1;k<6;k++)for(let l=k+1;l<6;l++){
   const key=((p[i]*88+p[j])*88+p[k])*88+p[l];
   if(check){if(seen4.has(key))return true;}else seen4.add(key);
 }
 return false;
};
const points=[];
for(const S of seed){const L=ar(22,i=>i+2).filter(i=>S>>i&1);
 for(const a of A){const p=L.map((r,k)=>4*(r-2)+a[k]);if(occupiesFour(p,true))throw Error("seed conflict");occupiesFour(p,false);points.push(p);}}
shapes.sort((a,b)=>pc(a.D)-pc(b.D));
for(const {a,b,c} of shapes)for(const x of a)for(const y of b)for(const z of c){
 const p=[...x,...y,...z];if(!occupiesFour(p,true)){occupiesFour(p,false);points.push(p);}
}
return {shapeCount:shapes.length,points};

};
const {shapeCount,points}=gen();
assert.equal(shapeCount,248);assert.equal(points.length,3624);
const group=c=>c<14?0:1+Math.floor((c-14)/16), row=c=>2+(c>>2);
const seen4=new Set(),cut=new Set(),hist={};
for(const p of points){
 assert.equal(p.length,6);
 for(let i=0;i<6;i++)assert(Number.isInteger(p[i])&&p[i]>=0&&p[i]<88&&(!i||p[i-1]<p[i]));
 assert(ar(22,r=>p.filter(c=>row(c)===r+2).length).every(x=>x<=2));
 assert(new Set(p.map(group)).size<=3);
 const U=p.reduce((s,c)=>s|1<<row(c),0);
 const D=p.reduce((s,c,i)=>s|(i&&row(p[i-1])===row(c)?1<<row(c):0),0);
 hist[pc(D)]=(hist[pc(D)]??0)+1;
 for(const [a,e,o] of T){const F=e|o;assert(2*pc(U&a)+pc(U&F)+pc(D&F)<=6);}
 for(const [s] of Z){const h=pc(s&U)+pc(s&D);assert(h*h<=3*pc(s));}
 for(let j=0;j<6;j++)for(let i=0;i<j;i++){
   if(group(p[i])===group(p[j]))cut.add(88*p[i]+p[j]);
   for(let k=j+1;k<6;k++)for(let l=k+1;l<6;l++){
     const key=((88*p[i]+p[j])*88+p[k])*88+p[l];assert(!seen4.has(key));seen4.add(key);
   }
 }
}
assert.equal(hist[0],1804);assert.equal(cut.size,358);
const codes=ar(64,i=>i).filter(i=>pc(i)%2===0);assert.equal(codes.length,32);
for(let i=0;i<codes.length;i++)for(let j=0;j<i;j++)assert(pc(codes[i]^codes[j])>=2);
const prior=JSON.parse(raw[4]).claims;
const old=prior.map(x=>x.lower_bound);assert.deepEqual(old,[74147456772,117097130108]);
const added=points.length*32,removed=cut.size*4,delta=added-removed;
assert.equal(delta,114536);
output({old,shapeCount,hist,supports:points.length,rootPairs:cut.size,added,removed,delta,new:old.map(n=>n+delta)});

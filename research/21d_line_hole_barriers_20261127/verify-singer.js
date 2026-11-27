"use strict";
const rot=x=>((x<<1)&((1<<21)-1))|(x>>>20);
function mul(a,b){let o=0;for(let i=0;i<6;i++){if(b&1)o^=a;b>>=1;a<<=1;if(a&64)a^=67;}return o;}
function pow(a,e){let b=1;while(e){if(e&1)b=mul(b,a);a=mul(a,a);e>>=1;}return b;}
if(new Set(Array.from({length:63},(_,i)=>pow(2,i))).size!==63)throw Error('primitive field seed');
const pc=x=>{let k=0;while(x){k++;x&=x-1;}return k;};
let sig=0;for(let i=0;i<21;i++){const a=pow(2,i);if(!(a^pow(a,4)^pow(a,16)))sig|=1<<i;}
const g=[];for(let i=0,r=sig;i<21;i++,r=rot(r))g.push(r);
if(new Set(g).size!==21||g.some(x=>pc(x)!==5)||g.some((x,i)=>g.slice(i+1).some(y=>pc(x&y)!==1)))throw Error('projective lines');
let C=new Set([0]);for(const x of g)for(const y of [...C])C.add(x^y);if(C.size!==1024)throw Error('line span');
const pop=x=>{let c=0;while(x){c++;x&=x-1n;}return c;};
function prove(step){
 const shift=x=>{for(let i=0;i<step;i++)x=rot(x);return x;};
 const orbits=[], cls=new Map;
 for(const x of C){if(cls.has(x))continue;let y=x,o=[];do{o.push(y);cls.set(y,orbits.length);y=shift(y);}while(y!==x);orbits.push(o);}
 const n=orbits.length,B=orbits.map((_,i)=>1n<<BigInt(i)),ALL=(1n<<BigInt(n))-1n;
 const adj=orbits.map(([x])=>g.reduce((a,v)=>a|B[cls.get(x^v)],0n)), close=adj.map((a,i)=>a|B[i]), cost=orbits.map(a=>a.length);
 let F0=0n;for(let j=0;j<n;j++)if(!(adj[j]&B[j]))F0|=B[j];let nodes=0;
 function f(U,F,w){
  nodes++;if(!U)return true;
  let best=n+1,list=0n;
  for(let j=0;j<n;j++)if(U&B[j]){
   const z=close[j]&F;if(!z)return false;const k=pop(z);if(k<best){best=k;list=z;}
  }
  const remaining=63-w,total=pop(U);let opts=[],passes=false;
  for(let j=0;j<n;j++)if((F&B[j])&&cost[j]<=remaining){const c=pop(close[j]&U);if(c){opts.push([j,c]);if(c*remaining>=cost[j]*total)passes=true;}}
  if(!passes)return false;
  opts.sort((a,b)=>b[1]*cost[a[0]]-a[1]*cost[b[0]]||a[0]-b[0]);
  let FF=F;for(const [j] of opts)if(list&B[j]){
   if(f(U&~close[j],FF&~close[j],w+cost[j]))return true;FF&=~B[j];
  }return false;
 }
 if(f(ALL,F0,0))throw Error('found an unexpectedly small pattern; do not claim exclusion');
 return {step,orbits:n,nodes};
}
console.log({sig,scope:[prove(1),prove(3)]});

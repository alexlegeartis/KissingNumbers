'use strict';
const assert=require('node:assert/strict');
const pc=(n)=>{n=n-((n>>>1)&0x55555555);n=(n&0x33333333)+((n>>>2)&0x33333333);return (((n+(n>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
const bl=(b)=>31-Math.clz32(b);
const R=[...new Set(Array.from({length:46},(_,j)=>(j+1)**2%47))].reduce((s,i)=>s|(1n<<BigInt(i)),0n);
const mask=(1n<<47n)-1n, B=Array(48).fill(0n);
for(let j=0;j<47;j++){
 let x=((R<<BigInt(j))|(R>>BigInt(47-j)))&mask;x|=1n<<47n;
 for(let i=47;i>=0;i--)if((x>>BigInt(i))&1n){if(B[i])x^=B[i];else{B[i]=x;break;}}
}
assert.equal(B.filter(Boolean).length,24);
const pref=Array(8).fill(0n),N=Array(40).fill(0n);
for(let x of B.filter(Boolean)){
 for(let j=7;j>=0;j--)if((x>>BigInt(j))&1n){if(pref[j])x^=pref[j];else{pref[j]=x;x=0n;break;}}
 if(x){x>>=8n;for(let j=39;j>=0;j--)if((x>>BigInt(j))&1n){if(N[j])x^=N[j];else{N[j]=x;break;}}}
}
assert(pref.every(Boolean));assert.equal(N.filter(Boolean).length,16);
const syn=(z)=>{for(let j=39;j>=0;j--)if((z>>BigInt(j))&1n&&N[j])z^=N[j];let s=0,k=0;for(let j=0;j<40;j++)if(!N[j]){s|=Number((z>>BigInt(j))&1n)<<k;k++;}return s;};
const C=Array.from({length:40},(_,j)=>syn(1n<<BigInt(j))),T=[0];
for(let i=0;i<8;i++){let x=1n<<BigInt(i);for(let j=7;j>=0;j--)if((x>>BigInt(j))&1n)x^=pref[j];T.push(syn(x>>8n));}
assert(C.slice(0,24).every((v,i)=>v===1<<i));assert.equal(new Set(T).size,9);
const masks=[0x55555555,0x33333333,0x0f0f0f0f,0x00ff00ff], L=1<<19;
let x=new Uint32Array(L),y=new Uint32Array(L);for(const v of T)x[v>>>5]|=1<<(v&31);
const out={};for(let r=1;r<=9;r++){
 y.set(x);
 for(const e of C){const i=e>>>5,k=e&31;
  if(!k){for(let a=0;a<L;a++)y[a^i]|=x[a];continue;}
  for(let a=0;a<L;a++){let w=x[a];if(k&1)w=((w&masks[0])<<1)|((w>>>1)&masks[0]);if(k&2)w=((w&masks[1])<<2)|((w>>>2)&masks[1]);if(k&4)w=((w&masks[2])<<4)|((w>>>4)&masks[2]);if(k&8)w=((w&masks[3])<<8)|((w>>>8)&masks[3]);if(k&16)w=(w<<16)|(w>>>16);y[a^i]|=w;}
 }
 [x,y]=[y,x];if(r>=7){let bad=[];for(let i=0;i<L;i++){let w=~x[i];while(w){const b=bl(w&-w);bad.push(32*i+b);w^=1<<b;}}out[r]=bad;}
}
assert.deepEqual([out[7].length,out[8].length,out[9].length],[1414,28,0]);
const lim=65536,rg=new Int32Array(lim),wt=new Uint8Array(lim);
for(let u=1;u<lim;u++){let b=u&-u;rg[u]=rg[u^b]^C[24+bl(b)];wt[u]=wt[u^b]+1;}
function calc(q,deep){
 const c=new Int32Array(40);let m=0,total=0;
 outer:for(let j=0;j<lim;j++){
  const u=(40503*j+q)&65535,base=q^rg[u],w0=wt[u];if(w0>9)continue;
  for(let k=deep?0:1;k<9;k++){let e=base^T[k],w=pc(e)+w0;if(deep?w!==9:w>9)continue;
   if(deep)assert(k!==0);m++;total+=w;let f=u;
   while(e){let b=e&-e;c[bl(b)]++;e^=b;}while(f){let b=f&-f;c[24+bl(b)]++;f^=b;}
   if(!deep&&m===96)break outer;
  }
 }
 const a=c.reduce((s,d)=>s+(40*d-total)**2,0),den=16000*m*m;
 assert(m>=(deep?397:72)&&total<=9*m&&a*(deep?900:25)<den);
 return {m,a,den};
}
const bests=[];
for(const deep of [false,true]){let lo=Infinity,hi=0,best={q:0,a:0,den:1,m:0};
 for(const q of out[deep?8:7]){const r=calc(q,deep);lo=Math.min(lo,r.m);hi=Math.max(hi,r.m);if(BigInt(r.a)*BigInt(best.den)>BigInt(best.a)*BigInt(r.den))best={q,...r};}
 bests.push({lo,hi,best});
}
assert.deepEqual(bests.map(v=>[v.lo,v.hi,v.best.q,v.best.m,v.best.a]),[[72,96,16614820,80,1324800],[397,424,16623012,405,2605400]]);
const U=666000000n-649n**2n*3360n,V=50000000n+649n**2n*160n;
assert(U<0n&&V>0n&&41n*V**2n>U**2n);
assert(41n*1000000n>80n*649n**2n);
const w=649n*9919n-500n*10081n;
assert(w>0n&&w**2n>(1000000n-649n**2n)*1800n**2n);
assert(11n*9919n-4n*1800n>10n*10081n);
assert(33n*247n-2n*96n>30n*265n);
console.log('verified syndrome coverage, two centroid gates and angular constants',bests);

#!/usr/bin/env node
import assert from "node:assert/strict";
const C="{\n  \"source\": \"Yang et al., arXiv:2609.35051v1, Appendix B, Tables 11 and 12\",\n  \"source_url\": \"https://arxiv.org/html/2609.35051v1\",\n  \"license\": \"CC BY 4.0\",\n  \"length\": 38,\n  \"generators_hex\": [\n    \"2000C0DE0B\", \"1000AA229D\", \"8009F5CD6\", \"4004FAE6B\",\n    \"200ED9AAD\", \"100BC80CE\", \"805E4067\", \"40E56DAB\",\n    \"20B8FB4D\", \"1096303E\", \"84B181F\", \"4EFC197\", \"2BDAD53\", \"1949B31\"\n  ],\n  \"representatives_hex\": [\n    \"0\", \"C6E94\", \"DF6ABF\", \"27B441\", \"D3CE66\", \"41BBBB\",\n    \"BC1E1F\", \"52727\", \"9951C8\", \"7001EB\", \"3F6969\"\n  ]\n}\n";
const z=(()=>{const f=JSON.parse(C),G=f.generators_hex.map(v=>BigInt("0x"+v)),R=f.representatives_hex.map(v=>BigInt("0x"+v));
const b=Array(36).fill(0n),tail=[0n,0n],u=(1n<<36n)-1n;
const put=z=>{for(let i=35;i>=0;i--)if(z>>BigInt(i)&1n){if(b[i])z^=b[i];else{b[i]=z;return;}}throw Error("dependent");};
for(let x of G){let z=x,keep=true;for(let i=0;i<2;i++)if(z>>BigInt(36+i)&1n){
 if(tail[i])z^=tail[i];else{tail[i]=z;keep=false;break;}
}if(keep)put(z);}
if(!tail[0]||!tail[1]||R.some(x=>x>>36n))throw Error("tail");
if(tail[0]>>37n&1n)tail[0]^=tail[1];
const X=b.flatMap((v,i)=>v?[]:[i]);if(X.length!==24)throw Error("rank");
const red=x=>{let z=x&u;for(let i=35;i>=0;i--)if(b[i]&&(z>>BigInt(i)&1n))z^=b[i];return X.reduce((r,k,i)=>r|(Number(z>>BigInt(k)&1n)<<i),0);};
const D=Array.from({length:36},(_,i)=>red(1n<<BigInt(i))),S=R.filter((_,i)=>i!==7).map(red),h0=red(tail[0]),h1=red(tail[1]);if(!h0||!h1||h0===h1)throw Error("tail shift");
const dist=new Uint8Array(1<<24).fill(255),queue=new Uint32Array(dist.length);
let s=0,t=0;for(let r of S){dist[r]=0;queue[t++]=r;}
const hist=[S.length];
while(s<t){const r=queue[s++],d=dist[r]+1;for(let j=0;j<D.length;j++){const a=r^D[j];if(dist[a]===255){dist[a]=d;queue[t++]=a;hist[d]=(hist[d]||0)+1;}}}
if(t!==1<<24)throw Error("coverage");
const off=[0,h0,h1,h0^h1],T=[[9,9,9,9],[9,10,8,9],[10,9,6,8],[8,9,9,10],[6,8,10,9]], Z=[0,1,4,-1,-4];
const H=[0n];for(let z of b)if(z)H.push(...H.map(x=>x^z));
const lift=q=>X.reduce((v,k,i)=>q>>i&1?v|(1n<<BigInt(k)):v,0n);
const code=H.map(x=>[Number(x&0xffffffffn)|0,Number(x>>32n)|0]);
const pop=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
if(X.some((v,i)=>v!==i))throw Error("free coordinates");
const histDiff=new Map();
const MIN=(a,b)=>{const z=a^b;if(histDiff.has(z))return histDiff.get(z);let d=36;for(let [u,v] of code){let r=pop(u^z)+pop(v);d=Math.min(d,r);}histDiff.set(z,d);return d;};
const W=[], tally=[0,0,0,0,0];
for(let q=0;q<dist.length;q++){
 let dd=[dist[q],dist[q^h0],dist[q^h1],dist[q^(h0^h1)]];
 for(let j=0;j<5;j++)if(dd.every((v,i)=>v>=T[j][i])){tally[j]++;if(W.length<20000)W.push([q,Z[j],dd]);}
}
if(W.length!==tally.reduce((s,p)=>s+p,0))return {tally,hist,skipped:true,off,S};
const edge=(a,b)=>{const [q,z]=a,[r,y]=b;const d=MIN(q,r), dot=144-8*d+z*y;return dot<=0||4*dot*dot<=(144+z*z)*(144+y*y);};
let best=[];
if(W.length<=200){
 const E=W.map(x=>W.map(y=>edge(x,y)));
 function search(P,Q){if(Q.length>best.length)best=Q;if(best.length>=5)return true;
  while(P.length+Q.length>best.length){let h=P.pop(),E1=P.filter(x=>E[h][x]);if(search(E1,[...Q,h]))return true;}
  return false;}
 search([...W.keys()],[]);
 return {tally,hist,best:best.map(i=>W[i]),W,off,S};
}
return {tally,hist,W:W.slice(0,30),off,S};})();
assert.deepEqual(z.tally,[4,39,39,39,39]);
assert.equal(z.W.length,160);
assert.deepEqual(z.hist,[10,360,6300,71400,589050,3121552,7309882,5194516,483328,779,38,1]);
assert.equal(z.best.length,4);
console.log(JSON.stringify({status:"PASS",hist:z.hist,bankTally:z.tally,compatibilityClique:z.best,possibleGain:0}));
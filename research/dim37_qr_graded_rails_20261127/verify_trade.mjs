#!/usr/bin/env node
import assert from "node:assert/strict";
const C="{\n  \"source\": \"Yang et al., arXiv:2609.35051v1, Appendix B, Tables 11 and 12\",\n  \"source_url\": \"https://arxiv.org/html/2609.35051v1\",\n  \"license\": \"CC BY 4.0\",\n  \"length\": 38,\n  \"generators_hex\": [\n    \"2000C0DE0B\", \"1000AA229D\", \"8009F5CD6\", \"4004FAE6B\",\n    \"200ED9AAD\", \"100BC80CE\", \"805E4067\", \"40E56DAB\",\n    \"20B8FB4D\", \"1096303E\", \"84B181F\", \"4EFC197\", \"2BDAD53\", \"1949B31\"\n  ],\n  \"representatives_hex\": [\n    \"0\", \"C6E94\", \"DF6ABF\", \"27B441\", \"D3CE66\", \"41BBBB\",\n    \"BC1E1F\", \"52727\", \"9951C8\", \"7001EB\", \"3F6969\"\n  ]\n}\n";
const Z=(()=>{
const f=JSON.parse(C),G=f.generators_hex.map(v=>BigInt("0x"+v)),R=f.representatives_hex.map(v=>BigInt("0x"+v));
const b=Array(36).fill(0n),tail=[0n,0n],u=(1n<<36n)-1n;
const put=z=>{for(let i=35;i>=0;i--)if(z>>BigInt(i)&1n){if(b[i])z^=b[i];else{b[i]=z;return;}}throw Error("dependent");};
for(let x of G){let z=x,keep=true;for(let i=0;i<2;i++)if(z>>BigInt(36+i)&1n){
 if(tail[i])z^=tail[i];else{tail[i]=z;keep=false;break;}
}if(keep)put(z);}
if(!tail[0]||!tail[1]||R.some(x=>x>>36n))throw Error("tail");
if(tail[0]>>37n&1n)tail[0]^=tail[1];
const X=b.flatMap((v,i)=>v?[]:[i]);if(X.length!==24)throw Error("rank");
const red=x=>{let z=x&u;for(let i=35;i>=0;i--)if(b[i]&&(z>>BigInt(i)&1n))z^=b[i];return X.reduce((r,k,i)=>r|(Number(z>>BigInt(k)&1n)<<i),0);};
const D=Array.from({length:36},(_,i)=>red(1n<<BigInt(i))),S=R.map(red),h0=red(tail[0]),h1=red(tail[1]);if(!h0||!h1||h0===h1)throw Error("tail shift");

const H=[0n];for(let z of b)if(z)H.push(...H.map(x=>x^z));
if(H.length!==4096)throw Error("H0");
const off=[0,h0,h1,h0^h1];
const lift=q=>X.reduce((v,k,i)=>q>>i&1?v|(1n<<BigInt(k)):v,0n);
const MIN=(a,b)=>{const z=lift(a^b);let d=36;for(let h of H){let x=h^z,r=0;for(;x;x&=x-1n)r++;d=Math.min(d,r);}return d;};
const rows=[[12518842,4],[8379313,-4],[14022956,4],[1386279,-4]], old=[4,1,-4,-1],conf=[];
const fit=(d,z,m)=>{let x=144-8*d+z*m;return x<=0||4*x*x<=(144+z*z)*(144+m*m);};
for(let [q,z] of rows){
 let hits=[];for(let e=0;e<4;e++)for(let i=0;i<11;i++){let d=MIN(q,S[i]^off[e]);if(!fit(d,z,old[e]))hits.push([e,i,d]);}
 conf.push(hits);
}
const pairs=rows.map(([q],i)=>rows.map(([r],j)=>{let d=MIN(q,r);if(i!==j&&!fit(d,rows[i][1],rows[j][1]))throw Error("new cross");return d;}));
const minH0=Math.min(...H.slice(1).map(z=>{let n=0;while(z){z&=z-1n;n++;}return n;}));
return {rows,conf,pairs,minH0};
})();
assert.deepEqual(Z.conf,[[[1,7,1]],[[3,7,1]],[[2,7,1]],[[0,7,1]]]);
assert.deepEqual(Z.pairs,[[0,11,10,11],[11,0,11,10],[10,11,0,11],[11,10,11,0]]);
assert(Z.minH0>=12);
assert(4*32**2<=160*32 && 4*16**2<=160*32);
assert(2*1536**2<=2304**2);
assert(4*4**2<=160);
console.log(JSON.stringify({status:"PASS",...Z,removed:4*4096,added:4*4096}));
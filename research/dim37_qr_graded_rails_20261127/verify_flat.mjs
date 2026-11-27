#!/usr/bin/env node
import assert from "node:assert/strict";
import fs from "node:fs";
import {createHash} from "node:crypto";
// Optional argv[2] is the parent read-only worktree.
const C="{\n  \"source\": \"Yang et al., arXiv:2609.35051v1, Appendix B, Tables 11 and 12\",\n  \"source_url\": \"https://arxiv.org/html/2609.35051v1\",\n  \"license\": \"CC BY 4.0\",\n  \"length\": 38,\n  \"generators_hex\": [\n    \"2000C0DE0B\", \"1000AA229D\", \"8009F5CD6\", \"4004FAE6B\",\n    \"200ED9AAD\", \"100BC80CE\", \"805E4067\", \"40E56DAB\",\n    \"20B8FB4D\", \"1096303E\", \"84B181F\", \"4EFC197\", \"2BDAD53\", \"1949B31\"\n  ],\n  \"representatives_hex\": [\n    \"0\", \"C6E94\", \"DF6ABF\", \"27B441\", \"D3CE66\", \"41BBBB\",\n    \"BC1E1F\", \"52727\", \"9951C8\", \"7001EB\", \"3F6969\"\n  ]\n}\n";
if(process.argv[2]){
 const b=fs.readFileSync(process.argv[2]+"/verifications/improved/dim37-qr-flat-repack/dense_parent38.json");
 assert.equal(createHash("sha1").update("blob "+b.length+"\0").update(b).digest("hex"),"3343690a0bd2ebdf7b0c30fd80bee29b0fe8702a");
 assert.equal(b.toString(),C);
}
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
const D=Array.from({length:36},(_,i)=>red(1n<<BigInt(i))),S=R.map(red),h0=red(tail[0]),h1=red(tail[1]);if(!h0||!h1||h0===h1)throw Error("tail shift");
const dist=new Uint8Array(1<<24).fill(255),queue=new Uint32Array(dist.length);
let s=0,t=0;for(let r of S){dist[r]=0;queue[t++]=r;}
const hist=[S.length];
while(s<t){const r=queue[s++],d=dist[r]+1;for(let j=0;j<D.length;j++){const a=r^D[j];if(dist[a]===255){dist[a]=d;queue[t++]=a;hist[d]=(hist[d]||0)+1;}}}
if(t!==1<<24)throw Error("coverage");
const off=[0,h0,h1,h0^h1];
const pure=[],tops=[];let survivors=0;
for(let q=0;q<dist.length;q++){
 const a=dist[q],b=dist[q^h0],c=dist[q^h1],d=dist[q^h0^h1],m=Math.min(a,b,c,d);
 pure[m]=(pure[m]||0)+1;
 if(a===10)tops.push([q,a,b,c,d]);
 if(a>=10&&b>=9&&c>=6&&d>=8)survivors++;
}
return {hist,pure,tops,survivors,off};})();
assert.deepEqual(z.hist,[11,396,6930,78540,647955,3394353,7403146,4915178,330564,141,2]);
assert.deepEqual(z.pure,[11,396,6930,78540,620655,2321973,1164862,937].map(x=>x*4));
assert.deepEqual(z.off,[0,11149981,12639755,7011478]);
assert.deepEqual(z.tops,[[12518842,10,1,9,10],[14022956,10,9,1,10]]);
assert.equal(z.survivors,0);
console.log(JSON.stringify({status:"PASS",...z}));
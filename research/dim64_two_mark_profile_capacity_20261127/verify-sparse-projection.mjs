#!/usr/bin/env node
import assert from "node:assert/strict";
import fs from "node:fs";
import {createHash} from "node:crypto";
// Optional argv[2] points to read-only worktree of pinned parent.
const cubic="0 1241 18626 36507 41078 73094 86429 87364 121823 128306 138938 152225 153208 186595 194062 197436 198629 216062 232871 238410 273268 274349 287670 297730 304623 341746 359145 359984 384582 389291 406990 424405 425228 449914 455575 469064 470161 483466 493630 501459 529653 539160 563822 564919 578220 606641 618332 624616 642035 642858 672397 683104 689364 706767 707606 726985 735524 760146 761227 774544 790892 802689 831770 832963 850392 869416 879301 896156 909447 910430 935700 944633 962464 975803 976738 987728 998589 1028646 1029887 1047268\n";
if(process.argv[2]){
 const name=process.argv[2]+"/verifications/improved/dim63-67-signed-overlap/cubic80.txt";
 const b=fs.readFileSync(name);
 assert.equal(createHash("sha1").update("blob "+b.length+"\0").update(b).digest("hex"),"5d374d5e49ab72f28d580008868330a47dcfa390");
 assert.equal(b.toString(),cubic);
}
const {out:shell,cubes:mother}=(()=>{function shell(cubic){const pop=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};const MON=[];for(let i=0;i<64;i++){let a=0,b=0;for(let j=0;j<32;j++){if((j&i)===i)a|=1<<j;if(((j+32)&i)===i)b|=1<<j;}MON.push([a,b]);}
  const cs=MON.filter((x,i)=>pop(i)===3),qs=MON.filter((x,i)=>pop(i)===2),as=MON.filter((x,i)=>pop(i)===1);const xor=(X,Y)=>[X[0]^Y[0],X[1]^Y[1]];const from=(S,v)=>S.reduce((a,b,i)=>v>>i&1?xor(a,b):a,[0,0]);const C=cubic.trim().split(/\s+/).map(i=>from(cs,+i)),Q=Array.from({length:32768},(_,i)=>from(qs,i)),A=Array.from({length:64},(_,i)=>from(as,i));const out=[];
  for(let k=0;k<C.length;k++){let [ca,cb]=C[k];for(let j=0;j<Q.length;j++){let qa=Q[j][0]^ca,qb=Q[j][1]^cb;for(let i=0;i<64;i++){let a=qa^A[i][0],b=qb^A[i][1],w=pop(a)+pop(b);if(w===16)out.push([a>>>0,b>>>0]);if(w===48)out.push([(~a)>>>0,(~b)>>>0]);}}}
  out.sort((x,y)=>x[1]-y[1]||x[0]-y[0]);if(out.length!==37996||out.some((x,i)=>i&&x[0]===out[i-1][0]&&x[1]===out[i-1][1]))throw Error("shell count");return {out,cubes:C.map(([a,b])=>[a>>>0,b>>>0])};
}
return shell(cubic);})();
const z=(()=>{const dg=(()=>{function gf(a,b){let t=0;for(let i=0;i<5;i++){if(b&1)t^=a;b>>>=1;a<<=1;if(a&32)a^=37;}return t;}
function tr(a){let t=0;for(let i=0;i<5;i++){t^=a;a=gf(a,a);}if(t>1)throw Error("trace");return t;}
const P=[];for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)P.push([i,j]);const dg=[];
for(let a=0;a<32;a++)for(let b=0;b<32;b++){
 const Q=(x,y)=>tr(gf(a,gf(x,y))^gf(b,gf(gf(x,x),y)^gf(x,gf(y,y))));
 const d=Array.from({length:5},(_,i)=>Q(1<<i,1<<i));let u=0;for(let [k,[i,j]] of P.entries())u|=(j<5?(Q(1<<i,1<<j)^(d[i]*d[j])):d[i])<<k;
 dg.push(u);
}
const l=Array(15).fill(0);for(let z of dg){let e=z^dg[0];for(let i=14;i>=0;i--)if(e>>i&1){
 if(l[i])e^=l[i];else{l[i]=e;break;}
}}
if(new Set(dg).size!==1024||l.filter(Boolean).length!==15)throw Error('DG');return dg;})();
const M=[];for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)M.push((1<<i)|(1<<j));
const mask=Array.from({length:64},(_,i)=>1n<<BigInt(i));
function proj(r, P){let v=0;for(let i=0;i<P.length;i++){let j=P[i];v|=((j<32?r[0]>>>j:r[1]>>>(j-32))&1)<<i;}return v;}
const mon=t=>mask.reduce((v,x,i)=>(i&t)===t?v^x:v,0n);
const U=M.map(mon);
const DI=dg.map(x=>U.reduce((v,p,i)=>x>>i&1?v^p:v,0n));
const AI=[0,...[0,1,2,3,4,5].map(i=>1<<i)].map(mon);
const bit=r=>[Number(r&0xffffffffn),Number(r>>32n)];
const [D,A]=[DI.map(bit),AI.map(bit)];
const counts=Array(4096).fill(0),cert=[];let seen=0,hist={};
for(let id=0;id<shell.length;id++){
 const [a,b]=shell[id],m=(BigInt(b)<<32n)|BigInt(a);
 const P=mask.flatMap((x,i)=>x&m?[i]:[]);
 if(P.length!==16)throw Error("support");
 const B=Array(16).fill(0);let rank=0;for(let a0 of A){let d=proj(a0,P);for(let k=15;k>=0;k--)if(d>>k&1){if(B[k])d^=B[k];else{B[k]=d;rank++;break;}}}
 if(rank<7)continue;
 const T=new Set();for(let q of D){let d=proj(q,P);for(let k=15;k>=0;k--)if(B[k]&&(d>>k&1))d^=B[k];T.add(d);}
 const n=T.size*(1<<rank);hist[n]=(hist[n]||0)+1;seen++;
 if(n===32768){const p=x=>{x^=x>>>16;x^=x>>>8;x^=x>>>4;return (0x6996>>>(x&15))&1;};if(![...A,...U.map(bit)].every(([x,y])=>p(a&x)===p(b&y)))throw Error('quadratic parity'); if(!mother.every(([x,y])=>p((a&x)>>>0)===p((b&y)>>>0)))throw Error('cubic parity');cert.push(id);for(let i=0;i<16;i++)for(let j=i+1;j<16;j++)counts[P[i]*64+P[j]]++;}
}
let z=[];for(let i=0;i<64;i++)for(let j=i+1;j<64;j++)z.push([counts[i*64+j],i,j]);
z.sort((a,b)=>a[0]-b[0]);
return {min:z.slice(0,20),max:z.slice(-4),seen,hist,cert};})();
assert.equal(shell.length,37996);
assert.equal(z.seen,35392);
assert.deepEqual(z.hist,{24576:64,28672:2176,30720:1792,32768:31360});
assert.equal(z.cert.length,31360);
assert.equal(z.min[0][0],1792);assert.equal(z.max.at(-1)[0],1944);
assert.equal(1792/8,224);assert.equal(224*2048-131072,327680);
console.log(JSON.stringify({status:"PASS",originalSupports:shell.length,rankSeven:z.seen,universal:z.cert.length,minIncidence:z.min[0][0],maxIncidence:z.max.at(-1)[0]}));
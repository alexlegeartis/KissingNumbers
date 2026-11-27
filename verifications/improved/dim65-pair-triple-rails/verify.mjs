#!/usr/bin/env node
import assert from "node:assert/strict";
import fs from "node:fs";
import {createHash} from "node:crypto";
// Optional argv[2] verifies pinned parent Git blobs in a read-only worktree.
const cfg={"i":1,"removed":[12,21,36,42],"m1":[21,42],"m2":[7,36],"a":8,"b":13,"sig":-1,"keep":[[4,-1],[5,1],[16,-1],[19,-1],[33,-1],[34,-1],[38,1],[39,-1],[44,-1],[44,1],[53,-1],[53,1],[54,1],[60,1],[70,-1],[75,1],[77,-1],[79,-1],[89,-1],[89,1],[91,1],[94,-1],[104,-1],[105,1],[107,-1],[107,1],[112,1],[114,-1],[114,1],[119,-1]]};
const assignment="21\n0 1 15 54\n12 21 36 42\n14 16 50 57\n15 4 17 28\n17 12 13 38\n18 3 52 58\n19 8 44 46\n31 20 35 45\n32 6 29 48\n35 26 37 59\n36 0 24 51\n39 14 22 49\n43 11 47 60\n54 9 56 63\n62 5 33 41\n63 19 40 43\n66 23 30 39\n70 10 55 61\n75 34 53 62\n77 18 27 32\n79 2 25 31\n1 3 4 5 16 17 19 24 25 26 27 32 33 34 35 38 39 44 51 53 54 60 61 65 70 72 73 75 77 78 79 88 89 91 92 94 104 105 107 108 112 114 115 117 119 120\n";
const heads="1 0 4225 15\n1 1 4225 15\n2 0 4268 15\n2 1 4654 15\n3 0 1577 15\n3 1 1577 15\n4 0 5795 15\n4 1 5795 15\n5 0 1157 15\n5 1 1157 15\n6 0 1647 15\n6 1 1647 15\n7 0 5837 15\n7 1 96 15\n8 0 5828 15\n8 1 1674 15\n9 0 4230 15\n9 1 4230 15\n10 0 194 15\n10 1 194 15\n11 0 68 15\n11 1 68 15\n12 0 4136 15\n12 1 4136 15\n13 0 1103 15\n13 1 1103 15\n14 0 513 15\n14 1 513 15\n15 0 5347 15\n15 1 5347 15\n16 0 4775 15\n16 1 4775 15\n17 0 1263 15\n17 1 5186 15\n18 0 5710 15\n18 1 5710 15\n19 0 518 15\n19 1 518 15\n20 0 1612 15\n20 1 4642 15\n21 0 1069 15\n21 1 1069 15\n22 0 4111 15\n22 1 4111 15\n23 0 4175 15\n23 1 4324 15\n24 0 551 15\n24 1 1548 15\n25 0 173 15\n25 1 5803 15\n26 0 1187 15\n26 1 1187 15\n27 0 1575 15\n27 1 1575 15\n28 0 673 15\n28 1 4841 15\n29 0 5350 15\n29 1 5350 15\n30 0 1220 15\n30 1 1220 15\n31 0 4777 15\n31 1 4777 15\n32 0 1132 15\n32 1 1132 15\n33 0 4292 15\n33 1 4292 15\n34 0 229 15\n34 1 4778 15\n35 0 1184 15\n35 1 1184 15\n36 0 5224 15\n36 1 5224 15\n37 0 5320 15\n37 1 5320 15\n38 0 1057 15\n38 1 1573 15\n39 0 105 15\n39 1 105 15\n40 0 79 15\n40 1 79 15\n41 0 687 15\n41 1 687 15\n42 0 238 15\n42 1 238 15\n43 0 5166 15\n43 1 5166 15\n44 0 1066 15\n44 1 5638 15\n45 0 138 15\n45 1 1607 15\n46 0 719 15\n46 1 719 15\n47 0 4138 15\n47 1 4802 15\n48 0 37 15\n48 1 37 15\n49 0 4161 15\n49 1 1 15\n50 0 1258 15\n50 1 1258 15\n51 0 1122 15\n51 1 5805 15\n52 0 1195 15\n52 1 1195 15\n53 0 4097 15\n53 1 4097 15\n54 0 1056 15\n54 1 1056 15\n55 0 555 15\n55 1 555 15\n56 0 5131 15\n56 1 5131 15\n57 0 5230 15\n57 1 552 15\n58 0 1030 15\n58 1 1030 15\n59 0 1539 15\n59 1 4234 15\n60 0 527 15\n60 1 527 15\n61 0 4676 15\n61 1 4676 15\n62 0 1255 15\n62 1 1255 15\n63 0 5731 15\n63 1 524 15\n";
const labels="0 1241 18626 36507 41078 73094 86429 87364 121823 128306 138938 152225 153208 186595 194062 197436 198629 216062 232871 238410 273268 274349 287670 297730 304623 341746 359145 359984 384582 389291 406990 424405 425228 449914 455575 469064 470161 483466 493630 501459 529653 539160 563822 564919 578220 606641 618332 624616 642035 642858 672397 683104 689364 706767 707606 726985 735524 760146 761227 774544 790892 802689 831770 832963 850392 869416 879301 896156 909447 910430 935700 944633 962464 975803 976738 987728 998589 1028646 1029887 1047268\n";
assert.equal(cfg.sig,-1);
const lines=assignment.trim().split("\n");
const old=lines.slice(1,22).map(s=>s.split(/\s+/).map(Number));
assert.equal(Number(lines[0]),21);
assert.deepEqual(old[cfg.i],cfg.removed);
const groups=[...old.filter((_,i)=>i!==cfg.i),[cfg.a,...cfg.m1],[cfg.b,...cfg.m2]];
assert.equal(groups.length,22);
assert.equal(new Set(groups.map(x=>x[0])).size,22);
assert(groups.every(x=>x[0]>=0&&x[0]<80));
assert.deepEqual(groups.flatMap(x=>x.slice(1)).sort((a,b)=>a-b),Array.from({length:64},(_,i)=>i));
if(process.argv[2]){
 const repo=process.argv[2],folder="verifications/improved/";
 const check=(name,sha)=>{
  const bytes=fs.readFileSync(repo+"/"+folder+name);
  assert.equal(createHash("sha1").update("blob "+bytes.length+"\0").update(bytes).digest("hex"),sha);
  return bytes.toString();
 };
 assert.equal(check("dim65-68-dg-diagonal-triples/triple_assignment.txt","b59d97220036ed6900db59c9837349212e4ff880"),assignment);
 assert.equal(check("dim65-68-complete-hyperplane-atlas/hyperplane_heads.txt","fe8ba58f9864d1e57ab07e80082f642e12e8a181"),heads);
 assert.equal(check("dim63-67-signed-overlap/cubic80.txt","5d374d5e49ab72f28d580008868330a47dcfa390"),labels);
 assert.equal(JSON.parse(fs.readFileSync(repo+"/"+folder+"dim65-68-dg-diagonal-triples/claims.json","utf8")).lower_bounds[65],424177492);
 check("dim65-dg-three-site-donors/verify_gates.py","9e6532e37ebdf15dbbbf1aa929e76d3dbcb60ad7");
}
const z=(()=>{const D=1000n,S=D*D;
// Angular skips are safe after independently relabeling the four-site packing.
// A leave link has min degree >=1 (63 minus twice the number of blocks through a pair).
// A maximal matching of size s in a 64-vertex link with m edges satisfies m>=64-s.
for (const [blocks, expected] of [[10741, 31], [10726, 29]]) {
   const uncovered = (65*64*63)/6 - 4*blocks;
   const edges = Math.floor(3*uncovered/65);
   assert.equal(63 % 2, 1);
   assert.equal(64-edges, expected);
   assert(expected >= 2); // two disjoint marked pairs
}
const base=335544320+40414*2048+10741*16+2;
const old=[
  [-5645n,1250n,1250n,0n],
  [1400n,-5500n,1280n,1070n],
  [1400n,1280n,-5500n,-1070n],
];
const group=sg=>[
  [-5679n,1250n,-sg*420n],
  [1250n,-5679n,-sg*420n],
  [2780n,2780n,sg*4290n],
];
const ps=[old,group(1n),group(-1n)];
const sum=x=>x.reduce((a,b)=>a+b,0n);
const ab=x=>x<0n?-x:x;
const norm=p=>BigInt(65-p.length)*S+sum(p.map(x=>x*x));
const delta=p=>sum(p.slice(0,-1).map(x=>ab(x-D)));
const count={};
function gate(k,t,n,m){if(t>0n&&4n*t*t>n*m)throw Error(k+" "+[t,n,m]);count[k]=(count[k]||0)+1;}
for(let a=0;a<ps.length;a++){
 for(let i=0;i<3;i++){
  let p=ps[a][i], n=norm(p), r=p.length-1; if(n>96n*S)throw Error("norm");
  gate("dense",(32n*D+delta(p))*D,n,64n*S);
  gate("tail pole",ab(p.at(-1))*D,n,S);
  for(let b=0;b<ps.length;b++)for(let j=0;j<3;j++){
   const q=ps[b][j],m=norm(q),tt=p.at(-1)*q.at(-1);
   if(a<b||(a===0&&b===0))gate("disjoint",32n*S+(delta(p)+delta(q))*D+tt,n,m);
   if(a!==b)continue;
   if(i!==j)gate("same owner",BigInt(64-r)*S+sum(p.map((x,k)=>x*q[k])),n,m);
   for(let e=0;e<(1<<r);e++){
    let t=16n*S+tt;for(let k=0;k<r;k++)t+=((e>>k&1)?1n:-1n)*(p[k]*q[k]-S);
    gate("distinct DG owners",t,n,m);
   }
  }
  for(let e=0;e<(1<<r);e++){
    let t=16n*D,z=0;for(let k=0;k<r;k++)if(e>>k&1){z++;t+=ab(p[k])-D;if(p[k]<0n)t-=2n*D;}
    gate("zero-tail weight 16",t*D,n,16n*S);
    for(let h=0;h<2;h++){
      if(a!==0&&i===2&&z===2&&h===1)continue; // the leave-link relabeling excludes this triple for each new pair
      let v=BigInt(4-z-h)*D+(h?ab(p.at(-1)):0n);
      for(let k=0;k<r;k++)if(e>>k&1)v+=ab(p[k]);
      gate("four-position",v*D,n,4n*S);
    }
  }
 }
}
const total=base+22*2*131072;
if(total!==424251218)throw Error("count");
return {count,base,total,prior:424177492,gain:total-424177492,norms:ps.map(arr=>arr.map(p=>+norm(p).toString()))};})();
assert.equal(z.base+21*2*131072+46*2*2048+2,z.prior);
assert.equal(z.prior+2*131072-46*2*2048-2,z.total);
const coeff=labels.trim().split(/\s+/).map(Number);
assert.equal(coeff.length,80); assert.equal(new Set(coeff).size,80);
const pop6=x=>{let r=0;for(;x;x&=x-1)r++;return r;};
const mon3=Array.from({length:64},(_,i)=>i).filter(i=>pop6(i)===3).map(i=>Array.from({length:64},(_,k)=>k).reduce((x,j)=>(j&i)===i?x|(1n<<BigInt(j)):x,0n));
assert.equal(mon3.length,20);
const cus=[cfg.a,cfg.b].map(c=>mon3.reduce((x,v,j)=>coeff[c]>>j&1?x^v:x,0n));
const cubes=cus.map(x=>[Number(x&0xffffffffn),Number(x>>32n)]);
const V=(()=>{const P=assignment.trim().split("\n");const old=P.slice(1,22).map(s=>s.split(/\s+/).map(Number));
const H=P[22].split(/\s+/).map(Number);
const head=heads.trim().split("\n").map(s=>s.split(/\s+/).map(Number));
const pop=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24;};
const mf=i=>{let u=0;for(let j=0;j<32;j++)if((i&j)===i)u|=1<<j;return u;};
const ml=Array.from({length:32},(_,i)=>[i,mf(i)]);
const mon2=ml.filter(([i])=>pop(i)<=2).map(([,a])=>a),mon3=ml.filter(([i])=>pop(i)>2).map(([,a])=>a);
if(mon2.length!==16||mon3.length!==16)throw Error("count");
const make=(src,n)=>src.reduce((v,x,i)=>n>>i&1?v^x:v,0);
const Q=Int32Array.from({length:65536},(_,i)=>make(mon2,i));
if(new Set(Q).size!==65536)throw Error("quadratics");
const cs=cubes;
const prof=sg=>[[-5679,1250,-sg*420],[1250,-5679,-sg*420],[2780,2780,sg*4290]];
const nr=p=>62000000+p.reduce((v,x)=>v+x*x,0);
function full(){
 const ans=[];
 for(let h of H){
  let [u,w,selection,k]=head[h];if(k!==15)throw Error("allocation");
  const I=Array.from({length:64},(_,i)=>i).filter(i=>(pop(u&i)&1)===w);
  const lcb=31-Math.clz32(u&-u);
  if(I.length!==32||I.some((s,i)=>((s&((1<<lcb)-1))|((s>>>(lcb+1))<<lcb))!==i))throw Error("lcb order");
  const base=make(mon3,selection),a=[];
  for(let g=0;g<cubes.length;g++){
   let [l,m]=cs[g],p=base;
   for(let k=0;k<32;k++){let j=I[k];if((j<32?l>>>j:m>>>(j-32))&1)p^=1<<k;}
   const low=[];
   for(let q of Q){
    const z=p^q,d=pop(z);if(d<4||(d&1))throw Error("inherited gate");
    if(d<8)low.push([z,d]);
   }
   a.push(low);
  }
  ans.push([I,a]);
 }
 return ans;
}
function allow(full,marks,cubic,sg){
 const P=prof(sg);
 return full.map(([I,U])=>{
  const a=I.indexOf(marks[0]),b=I.indexOf(marks[1]);
  return [-1,1].map(sh=>U[cubic].every(([m,d])=>P.every(p=>{
   const e=a<0?0:m>>>a&1?-1:1,f=b<0?0:m>>>b&1?-1:1;
   const cap=16000+Math.abs(p[0]-1000)+Math.abs(p[1]-1000)+2*Math.abs(p[2]);
   const n=nr(p); if(cap*cap>9*n)throw Error("far");
   const v=32000-2000*d+e*(p[0]-1000)+f*(p[1]-1000)+2*sh*p[2];
   if(!Number.isSafeInteger(v*v)||!Number.isSafeInteger(9*n))throw Error("exact");
   return v<=0||v*v<=9*n;
  })));
 });
}
return {old,H,full,allow};})();
const F=V.full();
const u=V.allow(F,cfg.m1,0,cfg.sig),v=V.allow(F,cfg.m2,1,-cfg.sig);
const retained=V.H.flatMap((h,i)=>[-1,1].flatMap((s,j)=>u[i][j]&&v[i][j]?[[h,s]]:[]));
assert.deepEqual(retained,cfg.keep);
assert.equal(retained.length,30);
const result=z.total+retained.length*2048;
assert.equal(result,424312658);
assert.equal(result-(10741-10726)*16,424312418);
console.log(JSON.stringify({status:"PASS",old:z.prior,new:result,gain:result-z.prior,atlas:retained.length,quadraticComparisons:V.H.length*2*65536,geometry:z.count}));

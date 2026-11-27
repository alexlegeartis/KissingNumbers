
const count=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24};
const mbase=(a,b)=>Array.from({length:64},(_,mask)=>mask).filter(mask=>{const p=count(mask);return p>=a&&p<=b});
const pairs=mbase(2,2),cube=mbase(3,3),aff=mbase(0,1);
const halves=mask=>{let u=0,v=0;for(let j=0;j<32;j++){if((j&mask)===mask)u|=1<<j;if(((j+32)&mask)===mask)v|=1<<j;}return [u,v]};
const mat=masks=>{const N=1<<masks.length,L=new Int32Array(N),H=new Int32Array(N);let j=0;for(const m of masks){const [l,h]=halves(m);let n=1<<j++;for(let a=0;a<n;a++){L[n+a]=L[a]^l;H[n+a]=H[a]^h}}return [L,H]};
const Q=mat(pairs), V=mat(aff), U=mat(cube);
const B=[0,1,20,21].map(Number);const fixedmask=(1<<0)|(1<<20);
const choices=Array.from({length:128},(_,i)=>i).filter(i=>!(V[0][i]&fixedmask));
const mul=(a,b)=>{let s=0;for(let i=0;i<5;i++){if(b&1)s^=a;b>>=1;a<<=1;if(a&32)a^=37}return s};
const tr=a=>{let r=0;for(let i=0;i<5;i++){r^=a;a=mul(a,a)}if(r>1)throw Error("tr");return r};
const dg=[];for(let a=0;a<32;a++)for(let b=0;b<32;b++){
 const bil=(x,y)=>tr(mul(a,mul(x,y))^mul(b,mul(mul(x,x),y)^mul(x,mul(y,y))));
 const diag=Array.from({length:5},(_,i)=>bil(1<<i,1<<i));let c=0;
 for(let j=0;j<15;j++){const p=pairs[j],i=31-Math.clz32(p&-p),k=31-Math.clz32(p);const q=k===5?diag[i]:bil(1<<i,1<<k)^(diag[i]*diag[k]);c|=q<<j;}
 dg.push(c);
}
const by=[[],[],[],[]];for(let i=0;i<128;i++)by[(V[0][i]&1)|((V[0][i]>>>19)&2)].push(i);
function min16(id) {
const l=U[0][id],h=U[1][id],f=[];
for(let q=0;q<32768;q++){
 const a=Q[0][q]^l,b=Q[1][q]^h;
 // a's values at marked sites require the affine matching prefix, not 00 itself.
 for(const i of by[(a&1)|((a>>>19)&2)]){let x=a^V[0][i],y=b^V[1][i];if(!(x&fixedmask)&&count(x)+count(y)===16)f.push([q,x,y]);}
}
return f;
}

const assert=(v,m)=>{if(!v)throw Error(JSON.stringify(m));};
const C="0 1241 18626 36507 41078 73094 86429 87364 121823 128306 138938 152225 153208 186595 194062 197436 198629 216062 232871 238410 273268 274349 287670 297730 304623 341746 359145 359984 384582 389291 406990 424405 425228 449914 455575 469064 470161 483466 493630 501459 529653 539160 563822 564919 578220 606641 618332 624616 642035 642858 672397 683104 689364 706767 707606 726985 735524 760146 761227 774544 790892 802689 831770 832963 850392 869416 879301 896156 909447 910430 935700 944633 962464 975803 976738 987728 998589 1028646 1029887 1047268".split(/\s+/).map(Number);
const TC=[0,359984,153208,469064,760146,976738,642858,831770], O=[0,1241,18626,36507,41078,73094,86429,87364,121823,128306];
const W=[[1,0,0,[743,3886,7685,7154,5731,412]],[2,1,0,[8186,1594,6557,5340,141,515]],[3,2,0,[5798,2705,6705,1535,2836,205]],[4,3,0,[5964,915,7403,500,654,1137]],[5,4,0,[5293,7512,6124,2499,6597,258]],[6,5,0,[7091,560,2233,4068,1215,206]],[7,6,0,[939,6003,3696,5265,4724,314]],[8,7,0,[3546,2989,646,5950,6976,105]],[9,8,0,[7351,7586,8102,4650,5341,250]],[10,9,0,[55,418,1260,5872,5248,2180]],[11,0,1,[1021,4178,3311,3962,5405,267]],[12,1,1,[3578,6514,187,1401,5456,609]],[13,2,1,[379,1416,4188,4608,7157,181]],[14,3,1,[4881,1875,5763,5607,3548,253]],[15,4,16,[3882,6798,994,5902,1083,112]],[16,5,1,[872,5728,863,7252,5850,1032]],[17,6,1,[1482,1272,2247,1804,3021,4210]],[18,7,1,[129,8047,6607,7174,6192,2156]],[19,8,1,[5445,5560,6357,1015,2619,70]],[21,9,4,[3246,1011,399,5300,4368,87]],[22,0,3,[7800,2565,2795,2773,6473,791]],[23,1,3,[3724,6974,6349,5927,3246,1025]],[24,2,2,[7179,157,88,5721,1653,803]],[25,3,2,[6102,7263,993,2327,7120,177]],[26,4,2,[872,2628,1574,5049,5493,113]],[27,5,8,[817,5949,733,8164,5702,1070]],[28,6,66,[2106,3292,3790,2361,5978,239]],[29,7,32,[6113,6502,5883,3880,3291,1102]],[30,8,16,[5113,3740,833,4492,7009,191]],[31,9,2,[1274,3215,6213,1778,6340,293]],[32,0,6,[3742,6106,1103,157,3578,77]],[33,1,8,[6004,2359,7617,5452,7607,1030]],[34,2,16,[1701,5249,5116,3669,630,131]],[35,3,4,[5184,6699,2615,2226,2844,55]],[36,4,8,[1800,1070,2548,746,5915,135]],[37,5,11,[7344,7769,2571,1461,3184,112]],[38,6,16,[7303,7572,6472,1850,6537,2065]],[39,7,3,[147,1578,7513,1861,5921,597]],[40,8,3,[7456,3563,5761,3601,4505,179]],[41,9,6,[917,1892,1017,1931,7652,2351]],[42,0,9,[1007,4383,1337,4595,6997,318]],[43,1,15,[4829,219,5226,2466,3119,349]],[44,2,21,[1018,5230,3623,7864,2588,260]],[45,3,5,[282,6255,4044,7521,2020,253]],[46,4,1,[7343,1565,1459,7632,5847,141]],[47,5,16,[4413,1596,967,7067,308,178]],[48,6,2,[6802,5687,7868,7957,5247,237]],[49,7,24,[5622,3009,2670,1000,3480,197]],[50,8,21,[6794,6591,2743,779,7203,212]],[51,9,13,[1625,1029,3181,314,5977,159]],[52,0,10,[4818,4986,3828,579,2710,187]],[53,1,33,[70,2358,5646,197,7481,1064]],[54,2,3,[5047,7566,3663,518,3488,1176]],[55,3,9,[5488,5852,5719,7385,5437,1060]],[56,4,6,[1574,3771,665,2438,7140,178]],[57,5,9,[2345,36,3705,7995,4054,537]],[58,6,19,[6800,6514,3280,8027,4350,308]],[59,7,2,[3668,645,1494,7822,8188,83]],[60,8,19,[869,3525,6794,4081,2948,201]],[61,9,5,[1341,7716,3301,3854,4896,147]],[62,0,65,[4167,1198,5053,2018,7626,37]],[63,1,29,[5598,705,1225,5284,6580,83]]];
const span=bas=>{let S=[0];for(const b of bas)S=S.concat(S.map(x=>x^b));return S;};
assert(C.length===80&&new Set(C).size===80,"parent labels");
assert(JSON.stringify(span([359984,153208,760146]))===JSON.stringify(TC),"8-cubic subgroup");
const macro=O.flatMap(p=>TC.map(c=>c^p));
assert(macro.length===80&&new Set(macro).size===80&&C.every(x=>macro.includes(x)),"ten disjoint cubic translates");
const cf=C.map(c=>[U[0][c],U[1][c]]);
const RM=mbase(0,2).map(halves);
const wt=([l,h])=>count(l)+count(h);
const ort=([l,h],[u,v])=>!((count(l&u)+count(h&v))&1);
assert(RM.every(a=>wt(a)%4===0)&&cf.every(a=>wt(a)%4===0),"mother generators doubly even");
assert(RM.every(a=>RM.every(b=>ort(a,b)))&&cf.every(a=>RM.every(b=>ort(a,b))&&cf.every(b=>ort(a,b))),"mother degree<=2 and cubics self-orthogonal");
const D=new Set(dg);
assert(D.size===1024&&D.has(0),"DG chart");
const K=dg.filter(x=>dg.every(y=>D.has(x^y)));
assert(K.length===32&&new Set(K).size===32&&K.every(a=>K.every(b=>K.includes(a^b))),"32 additive stabilizers");
const pivot=new Int32Array(15);
for(let x of K)for(let j=14;j>=0&&x;j--)if(x>>j&1){if(pivot[j])x^=pivot[j];else{pivot[j]=x;break}}
const free=Array.from({length:15},(_,j)=>j).filter(j=>!pivot[j]);
assert(JSON.stringify(free)===JSON.stringify([0,1,2,3,4,10,11,12,13,14]),"fixed coordinate convention");
function red(q){for(let i=14;i>=0;i--)if(pivot[i]&&(q>>i&1))q^=pivot[i];let t=0;for(let i=0;i<free.length;i++)t|=((q>>free[i])&1)<<i;return t;}
assert(K.every(k=>red(k)===0),"kernel K");
assert(choices.length===32&&by.every(x=>x.length===32),"rank two on marks");
const F=TC.map(min16);
assert(F[0].length===1457,"RM2 exact min16 words with marks fixed");
const seen=new Set(), sites=new Set(), marks=[0,20], counts=[], badcounts=[];
for(const [j,p,o,B] of W){
 assert(Number.isInteger(j)&&j>=0&&j<64&&!marks.includes(j)&&!sites.has(j),"ordinary unique site");sites.add(j);
 assert(Number.isInteger(p)&&p>=0&&p<O.length&&Number.isInteger(o)&&o>=0&&o<8192,"macro and translation");
 assert(B.every(b=>Number.isInteger(b)&&b>=0&&b<8192)&&B.length===6,"selector basis");
 const S=span(B);
 assert(S.length===64&&new Set(S).size===64,"selector rank six");
 const bad=new Set();
 for(let i=0;i<TC.length;i++)for(const [q,l,h]of F[i])if(!(((j<32?l:h)>>>(j%32))&1))bad.add((i<<10)|red(q));
 assert(!bad.has(0)&&S.every(x=>!bad.has(x)),["no weight16 when marks and spike agree",j]);
 for(const x of S){
   const t=x^o, u=t>>>10, q=t&1023, cc=O[p]^TC[u], key=cc*1024+q;
   assert(!seen.has(key),"disjoint owner tiles");
   seen.add(key);
 }
 counts.push(S.length);badcounts.push(bad.size);
}
assert(sites.size===62&&seen.size===3968,"all 62 sites; all 3968 different tiles");
const gates=[];
const gate=(name,d,a,b)=>{assert(d<=0||4n*d*d<=a*b,[name,String(d)]);gates.push(name)};
gate("new-new same spike with t=16 and opposite source spike signs",15n,77n,77n);
gate("new-new same spike with t>=20",37n,77n,77n);
gate("new-new different spike, parent t>=16",36n,77n,77n);
gate("old dense, parent t>=16 minus single marked difference",35n,77n,64n);
gate("old unmarked sparse, one opposite parity",17n,77n,16n);
gate("old single-mark sparse, one opposite body parity",16n,77n,16n);
gate("old body affine plane",7n,77n,4n);
gate("old head affine plane",5n,77n,4n);
gate("old body pole",4n,77n,1n);
gate("old tail pole",0n,77n,1n);
const prev=232191236, old=230159620, oldSpikes=62*32*1024, add=seen.size*K.length*choices.length;
assert(old+oldSpikes===prev&&add===4063232&&old+add===234222852,"exact ledger");
console.log(JSON.stringify({status:"PASS",prev,retained:old,oldSpikes,replacement:add,count:old+add,gain:old+add-prev,macroTileCounts:O.map(p=>[...seen].filter(x=>TC.includes((x>>>10)^p)).length),badCounts:[Math.min(...badcounts),Math.max(...badcounts)],min16Counts:F.map(x=>x.length),angleGates:gates}));
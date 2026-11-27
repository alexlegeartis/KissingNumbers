'use strict';
// 33D exact seven-high continuous paid-price atlas. Node v20+.
const assert=require('node:assert/strict'),fs=require('node:fs'),crypto=require('node:crypto');
const pc=x=>{x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return Math.imul((x+(x>>>4))&0xf0f0f0f,0x01010101)>>>24};
const alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
function base64(s){let x=0,k=0,out=[];for(const c of s){const v=alphabet.indexOf(c);if(v<0)continue;x=x<<6|v;k+=6;if(k>=8){k-=8;out.push(x>>>k&255);}}return Uint8Array.from(out);}
function inf(by,n){let at=0,k=0,b=0,w=0;const O=new Uint8Array(n);
const rd=(j)=>{while(k<j){if(at>=by.length)throw Error('EOF');b|=by[at++]<<k;k+=8;}const x=b&((1<<j)-1);b>>>=j;k-=j;return x};
const R=(L)=>{const N=new Int32Array(16),F=new Int32Array(16),T=Array.from({length:16},()=>new Map());for(const x of L)if(x)N[x]++;for(let i=1;i<16;i++)F[i]=(F[i-1]+N[i-1])<<1;for(let i=0;i<L.length;i++)if(L[i])T[L[i]].set(F[L[i]]++,i);return T;};
const sym=T=>{let x=0;for(let i=1;i<16;i++){x=x<<1|rd(1);if(T[i].has(x))return T[i].get(x);}throw Error('code');};
const LB=[3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258],LE=LB.map((x,i)=>i<8||i===28?0:(i>>2)-1);
const DB=[1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577],DE=DB.map((x,i)=>i<4?0:(i>>1)-1);
let done=0;while(!done){done=rd(1);const flag=rd(2);let L,D;
 if(flag===0){rd(k&7);const m=rd(16);if((m^rd(16))!==65535)throw Error('stored');for(let i=0;i<m;i++)O[w++]=rd(8);continue;}
 if(flag===1){L=R(Array.from({length:288},(_,i)=>i<144?8:i<256?9:i<280?7:8));D=R(Array(32).fill(5));}
 if(flag===2){const m=rd(5)+257,n=rd(5)+1,g=rd(4)+4,ord=[16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15],u=Array(19).fill(0);for(let i=0;i<g;i++)u[ord[i]]=rd(3);
  const T=R(u),S=[];while(S.length<m+n){const x=sym(T);if(x<16)S.push(x);else if(x===16)S.push(...Array(rd(2)+3).fill(S.at(-1)));else if(x===17)S.push(...Array(rd(3)+3).fill(0));else S.push(...Array(rd(7)+11).fill(0));}
  L=R(S.slice(0,m));D=R(S.slice(m));}
 if(!L)throw Error('invalid block');
 for(;;){const x=sym(L);if(x<256)O[w++]=x;else if(x===256)break;else {const z=LB[x-257]+rd(LE[x-257]),q=sym(D),s=DB[q]+rd(DE[q]);if(w<s)throw Error('bad distance');for(let y=0;y<z;y++){O[w]=O[w-s];w++;}}}
}
if(w!==n)throw Error('length');return O;}
function decode(arr){
 const by=Uint8Array.from(arr),V=new DataView(by.buffer);const rd=(j)=>V.getUint32(j,true),rh=(j)=>V.getUint16(j,true);
 let e=by.length-22;while(e>=0&&rd(e)!==0x06054b50)e--;if(e<0)throw Error('missing dir');
 const n=rh(e+10);let j=rd(e+16),Z={};for(let i=0;i<n;i++){
 if(rd(j)!==0x02014b50)throw Error('central dir');
 const sz=rd(j+20),out=rd(j+24),m=rh(j+28),x=rd(j+42),name=String.fromCharCode(...by.slice(j+46,j+46+m)),off=x+30+rh(x+26)+rh(x+28);
 const T=by.slice(off,off+sz),U=rh(j+10)===8?inf(T,out):T;
 if(U[0]!==147||U.length!==out)throw Error('NPY');
 const v=new DataView(U.buffer,U.byteOffset,U.byteLength),q=U[6]===1?10:12,end=q+(q===10?v.getUint16(8,true):v.getUint32(8,true));const head=String.fromCharCode(...U.slice(q,end)),D=head.match(/'descr': '([^']+)/)[1],shape=head.match(/'shape': \(([^)]*)/)[1].split(',').map(x=>x.trim()).filter(Boolean).map(Number);
 const cls={'<i2':Int16Array,'<i4':Int32Array,'<i8':BigInt64Array,'|i1':Int8Array,'<u2':Uint16Array,'<u4':Uint32Array,'<u8':BigUint64Array,'|u1':Uint8Array}[D];
 if(!cls)throw Error(D);
 const data=Array.from(new cls(U.buffer,U.byteOffset+end,(U.length-end)/cls.BYTES_PER_ELEMENT),Number);
 if(data.length!==shape.reduce((a,b)=>a*b,1))throw Error('shape');
 Z[name.replace('.npy','')]={shape,data};j+=46+m+rh(j+30)+rh(j+32);
 }
 return Z;
}


function makeNineData(D){
const P=Array.from({length:33},(_,i)=>1n<<BigInt(i)), B=D.supports.data.map(BigInt),A=P.map((p,i)=>i), C=new Set();
assert.equal(B.length,1837);assert(B.every(m=>pcB(m)===8));assert(D.parity.data.every(x=>x===0));
for(const m of B){const V=P.filter(p=>m&p);for(let a=0;a<8;a++)for(let b=a+1;b<8;b++)
 for(let c=b+1;c<8;c++)for(let d=c+1;d<8;d++)for(let e=d+1;e<8;e++)C.add(V[a]|V[b]|V[c]|V[d]|V[e]);}
let U=[];for(let a=0;a<29;a++)for(let b=a+1;b<30;b++)for(let c=b+1;c<31;c++)for(let d=c+1;d<32;d++)for(let e=d+1;e<33;e++){
 const m=P[a]|P[b]|P[c]|P[d]|P[e];if(!C.has(m))U.push(m);}
const counts=[U.length];
for(let k=6;k<=9;k++){
 const prev=new Set(U),V=new Set();
 for(const m of U)for(const p of P)if(p>m&&P.every(q=>!(m&q)||prev.has((m^q)|p)))V.add(m|p);
 U=[...V];counts.push(U.length);
}
assert.deepEqual(counts,[134464,79706,3044,197,18]);
assert(U.every(m=>m&P[32]));
assert.deepEqual(D.sparse.shape,[4348,33]);
const S=U.map(m=>A.filter(i=>m&P[i]));return S;
}
function pcB(x){let n=0;for(;x;x&=x-1n)n++;return n}

function makeCharts(D,S){
 const A=Array.from({length:33},(_,i)=>i),P=A.map(i=>1n<<BigInt(i)),B=D.supports.data.map(BigInt),VV=new Map();
 for(const s of S)for(const j of A)if(!s.includes(j)){const X=[...s,j].sort((a,b)=>a-b);VV.set(X.join(','),X);}
 assert.equal(VV.size,420);
 const R=[],H={};
 for(const X of VV.values()){const support=X.reduce((m,i)=>m|P[i],0n),F=[...new Set(B.map(s=>s&support).filter(s=>pcB(s)===5))];
  for(let a=0;a<10;a++)for(let b=a+1;b<10;b++)for(let c=b+1;c<10;c++){
   const V=[X[a],X[b],X[c]],m=P[X[a]]|P[X[b]]|P[X[c]],Q=new Set();
   if(F.some(f=>pcB(f&m)<2))continue;
   for(const f of F){const q=f&m;if(pcB(q)===2)Q.add(Number(q));}
   if(V.includes(32))Q.add(Number(m^P[32]));
   if(Q.size>=3)continue;
   const e=Q.size;if(!H[e])H[e]=[0,0];H[e][V.includes(32)?1:0]++;
   R.push({H:X.filter(i=>!V.includes(i)),L:V,con:[...Q],support:X}); 
  }
 }
 assert.deepEqual(H,{0:[10,0],1:[110,20],2:[23,8]});
 return R;
}
function rows(D,k){return Array.from({length:k},(_,i)=>D.sparse.data.slice(33*i,33*i+33));}

function atlas3(){
const vert=[
['lo','Z','q','T','hi','p'],
['lo','Z','q T','hi','p'],
['lo','Z','T','q','hi','p'],
['lo','Z T','q','hi','p'],
['lo','T','Z','q','hi','p'],
['lo T','hi p q Z'],
['T','lo','p','hi','q','Z'],
['T','p lo','hi','q','Z'],
['T','p','lo','hi','q','Z'],
['T','p','lo hi','q','Z']
];
const perms=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
const seen=new Set(),A=[];
for(let a=0;a<vert.length;a++){
 const V=vert[a],ids=Object.fromEntries(V.flatMap((v,i)=>v.split(' ').map(x=>[x,2*i])));
 for(let g=ids.lo;g<=ids.hi;g++){
 const lhs=x=>Math.sign(g-ids[x]);
 const pair=[[0,lhs('p'),-lhs('q')],[lhs('p'),0,1],[-lhs('q'),1,0]];
 function arm(n,t){if(t<0)return -arm(n.map(z=>-z),1);
 const pos=[],neg=[];n.forEach((x,i)=>{if(x>0)pos.push(i);if(x<0)neg.push(i);});
 if(t===1){
   if(pos.length===0||pos.length===1)return -1;
   if(pos.length===2)return neg.length===0?pair[pos[0]][pos[1]]:neg[0]===0?lhs('Z'):-1;
   return 1;
 }
 if(!pos.length)return neg.length?-1:0;
 if(!neg.length)return 1;
 if(pos.length===2)return neg[0]===2?lhs('T'):1;
 if(neg.length===2)return pos[0]===2?-lhs('T'):-1;
 const i=pos[0],j=neg[0],d=Math.sign(i-j);
 if((i===0&&j===1||i===1&&j===0)&&lhs('lo')===0)return 0;
 if((i===1&&j===2||i===2&&j===1)&&lhs('hi')===0)return 0;
 if(a===9)return 0;
 return d;
 }
 for(const p of perms){const table=[],pairs=[];
 for(let i=0;i<3;i++)for(let j=i+1;j<3;j++)pairs.push(pair[p[i]][p[j]]);
 for(let k=0;k<27;k++){
  const q=[k%3-1,(k/3|0)%3-1,(k/9|0)-1],n=[0,0,0];for(let i=0;i<3;i++)n[p[i]]=q[i];
  for(let t=-1;t<=1;t++)table.push(arm(n,t));
 }
 const key=table.join(',');if(seen.has(key))continue;seen.add(key);
 A.push({pair:pairs,table,repr:[a,g,p]});
 }
 }
}
return A;
}
function fullAtlas(D,S,C){
 const rank=v=>{const B=Array(10).fill(0);for(let x of v)while(x){const j=31-Math.clz32(x);if(B[j])x^=B[j];else{B[j]=x;break;}}return B.filter(x=>x!==0).length};
 const VV=new Map();
 for(const I of S)for(let i=0;i<33;i++)if(!I.includes(i)){
  const X=[...I,i].sort((a,b)=>a-b);VV.set(X.join(','),X);
 }
 assert.equal(VV.size,420);
 assert(D.dense_positions.data.every((x,i)=>x===i));
 for(const I of VV.values()){
  const P=I.filter(i=>i!==32);assert.equal(P.length,9);
  assert.equal(rank(D.dense_generator.data.map(v=>P.reduce((m,i,j)=>m|((v>>>i&1)?1<<j:0),0))),9);
 }
 const atl=atlas3();assert.equal(atl.length,127);
 const X=rows(D,2216).map(r=>r.map(x=>x/2));
 assert(X.every(r=>r.every(x=>[-1,0,1].includes(x))&&r.filter(x=>x!==0).length===8));
 const dist={};
 for(const {H,L,con} of C){
 const pp=[[L[0],L[1]],[L[0],L[2]],[L[1],L[2]]].map(([a,b])=>Number(1n<<BigInt(a)|1n<<BigInt(b)));
 const Z=atl.filter(A=>con.every(z=>A.pair[pp.indexOf(z)]<=0));assert(Z.length);
 const Y=X.flatMap(r=>{const h=H.map(j=>r[j]),l=L.map(j=>r[j]),a=h.reduce((s,z)=>s+Math.abs(z),0);return a>=3?[{h,l}]:[];});
 const comps=[];
 for(let v=0;v<8;v++){
   const plug=Array.from({length:27},(_,k)=>{let q=[k%3-1,(k/3|0)%3-1,(k/9|0)-1];for(let i=0;i<3;i++)if(v>>>i&1)q[i]=-q[i];return q[0]+1+3*(q[1]+1)+9*(q[2]+1);});
   for(const A of Z){const ids=[];for(let k=0;k<27;k++)for(let t=0;t<3;t++)if(A.table[3*plug[k]+t]>0)ids.push(k*3+t);comps.push(ids);}
 }
 for(let h=0;h<128;h++){
   let must=0;const B=new Int16Array(81);
   for(const {h:w,l} of Y){
    const a=w.reduce((s,z,j)=>s+((h>>>j&1)?-z:z),0),k=l[0]+1+3*(l[1]+1)+9*(l[2]+1);
    if(a>=6)must++;
    else if(a>=3)B[k*3+5-a]++;
   }
   let opt=9999;
   for(const ids of comps)opt=Math.min(opt,must+ids.reduce((sum,i)=>sum+B[i],0));
   dist[opt]=(dist[opt]||0)+1;
 }
 }
 assert.deepEqual(dist,{6:21,7:137,8:252,9:524,10:1046,11:1791,12:2033,13:2625,14:2163,15:2864,16:2061,17:2079,18:1457,19:1196,20:717,21:332,22:326,23:105,24:131,25:8,26:20});
 return dist;
}

function main(buffer){
 assert(Buffer.isBuffer(buffer));
 const sha=crypto.createHash('sha1').update('blob '+buffer.length+'\0').update(buffer).digest('hex');
 assert.equal(sha,'a282e197b71dc04ee5295ebfbeafa223fd51cac0');


 const D=decode(buffer), S=makeNineData(D), R=makeCharts(D,S);
 console.log(JSON.stringify({atlas:atlas3().length,charts:R.length,prices:fullAtlas(D,S,R)}));
}
main(fs.readFileSync(process.argv[2]));
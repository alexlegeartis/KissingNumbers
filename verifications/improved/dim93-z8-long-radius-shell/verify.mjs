// Exact standalone geometry replay; source and minimum premise in README.
const old=[[0,0,2,0,2,0,0],[0,4,2,0,-2,0,0],[0,32,-2,0,2,0,0],[0,36,-2,0,-2,0,0],[0,1,2,0,0,2,0],[0,5,2,0,0,-2,0],[0,33,-2,0,0,2,0],[0,37,-2,0,0,-2,0],[0,2,2,0,0,0,2],[0,6,2,0,0,0,-2],[0,34,-2,0,0,0,2],[0,38,-2,0,0,0,-2],[0,8,0,2,2,0,0],[0,12,0,2,-2,0,0],[0,40,0,-2,2,0,0],[0,44,0,-2,-2,0,0],[0,9,0,2,0,2,0],[0,13,0,2,0,-2,0],[0,41,0,-2,0,2,0],[0,45,0,-2,0,-2,0],[0,10,0,2,0,0,2],[0,14,0,2,0,0,-2],[0,42,0,-2,0,0,2],[0,46,0,-2,0,0,-2],[2,60,0,0,0,0,1],[2,47,0,0,2,0,-2],[2,63,0,0,2,0,2],[1,34,0,0,2,-2,0],[2,62,0,0,-2,0,2],[2,61,0,0,0,-2,2],[2,42,0,0,0,2,-2],[1,2,0,0,-1,0,0],[1,58,0,0,1,0,0],[2,36,2,2,0,0,0],[1,42,0,0,0,-1,0],[1,26,0,0,-2,2,0],[1,50,0,0,0,1,0],[2,46,0,0,0,-2,-2],[1,18,0,0,-2,-2,0],[2,4,2,-2,0,0,0],[2,44,0,0,0,0,-1],[2,52,-2,-2,0,0,0],[2,20,-2,2,0,0,0],[2,41,0,0,-2,0,-2],[2,57,0,0,0,2,2],[1,10,0,0,2,2,0]];
const A=Array.from({length:8},(_,k)=>Array.from({length:4},(_,x)=>1-2*(((k>>2)^((k&1)*(x>>1))^(((k>>1)&1)*(x&1)))&1)));
const B=A.map(a=>a.map((v,j)=>j?v:-v));
const C=A.map((a,k)=>[0,1,2,3].map(j=>j===k%4?(k<4?2:-2):0));
const T=[A,B,C],dot=(a,b)=>a.reduce((x,y,i)=>x+y*b[i],0),bf=(c,p,q)=>48-dot(T[c][p>>3],T[c][q>>3])-dot(T[c][p&7],T[c][q&7]);
const r=old.filter(x=>dot(x.slice(2),x.slice(2))===8);
const G=[[7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,0,0,0,0,0,0,0,3],[0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,0,0,0,0,0,0,3],[0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,0,0,0,0,0,3],[0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,0,0,0,0,3],[0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,0,0,0,3],[0,0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,0,0,3],[0,0,0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,0,3],[0,0,0,0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,0,3],[0,0,0,0,0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,0,3],[0,0,0,0,0,0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,0,3],[0,0,0,0,0,0,0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,0,3],[0,0,0,0,0,0,0,0,0,0,0,7,1,4,4,2,7,3,3,4,7,2,1,3]];
const U=G.slice(2), M22=(1<<22)-1;
const pc=x=>{x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0xf0f0f0f)*0x1010101)>>>24;};
const word=k=>{let z=new Uint8Array(24);for(let i=0;i<10;i++){let a=k>>>(2*i)&3;"a";for(let j=2;j<24;j++)z[j]=(z[j]+a*U[i][j])&7;}return z;};
const BM=U.map(row=>row.slice(2).reduce((a,x,i)=>a|((x&1)<<i),0));
const Q=[0];for(let b of BM)Q.push(...Q.map(q=>q^b));
const et=(x)=>{let w=0;for(let i=0;i<10;i++)w|=(x>>i&1)<<(2*i);return w;};
function distZ(z){let a=0,b=0;for(let j=2;j<24;j++){let x=z[j];if(x&3)a|=1<<(j-2);else if(x&4)b|=1<<(j-2);}
let wb=pc(a);if(wb>=16)return 2*wb;let m=5, mm=M22^a;
for(const t of Q){let d=pc((t^b)&mm);if(d<m)m=d;if(!m)break;}return 2*wb+4*m;}

const W=[0,925866,138681,778062,125631,582428,262875,669053,843288,490391,493044,894580,636966,77810,890680,357816,831563,792586,602376,17125,101086,235865,951219,399933,293269], Z=[0,1044849,104145,201198],
  sh=[0,732955,114648,366284,890643,625560,753704,1007067], pref=[8,16,24,32], signBases=[0,12,4,8];
const check=(s,m)=>{if(!s)throw Error(m);};
const norm=a=>dot(a,a);
check(G.length===12&&G.every((r,i)=>r.length===24&&r[i]===7&&r.slice(0,i).every(v=>!v)),'source rank');
check(Q.length===1024&&new Set(Q).size===1024&&Math.min(...Q.slice(1).map(pc))===8,'4G0 size and min');
check(old.length===46&&r.length===40&&old.filter(x=>norm(x.slice(2))===1).length===6,'retained dense count');
for(let i=0;i<46;i++){
 const a=old[i], u=a.slice(2), n=norm(u);
 check(n<=8&&2*(40+n)<=88+n,'old internally');
 for(let j=0;j<i;j++){const b=old[j],v=b.slice(2),k=a[0]===b[0]?bf(a[0],a[1],b[1]):44;
  const d=k+dot(u,v);check(d<=0||4*d*d<=(88+n)*(88+norm(v)),'old-old Gram');}
}
const ks=25, amp=14, den=5, R2=88*den*den+5*amp*amp;
const plus=(q,p)=>((q^p)^(((q&p)&0x55555)<<1));
const af=(p,z)=>{const x=7*(p>>3)&7,y=7*((p&7)-x*G[0][1])&7;const h=word(z);return h.map((a,i)=>(a+x*G[0][i]+y*G[1][i])&7);};
function mind(a,b,limit){
 let A=0,B=0;for(let j=2;j<24;j++){let v=(a[j]-b[j])&7;if(v&3)A|=1<<(j-2);else if(v&4)B|=1<<(j-2);}
 let ar=2*pc(A),M=M22^A;
 if(ar>=limit)return true;
 for(let q of Q)if(ar+4*pc((B^q)&M)<limit)return false;
 return true;
}
const wh=W.map(word);
check(W.length===ks&&new Set(W).size===ks,'literal quotient word count');
for(let i=0;i<ks;i++)for(let j=0;j<i;j++)check(mind(wh[i],wh[j],32),'intra-head exact coset distance');
for(let i=0;i<4;i++)for(let j=0;j<i;j++)if(pc(signBases[i]^signBases[j])===1)
 for(let a of W)for(let b of W)
  check(mind(af(pref[i],plus(a,Z[i])),af(pref[j],plus(b,Z[j])),24),'cross square edge exact coset distance');
const H=[];
for(let i=0;i<4;i++)for(let o=0;o<8;o++){
 const masks=[1,2,16], s=masks.reduce((a,m,j)=>a^(o>>j&1?m:0),signBases[i]);
 const v=Array.from({length:5},(_,j)=>(s>>j&1?-amp:amp));
 const cands=W.map(x=>plus(plus(x,sh[o]),Z[i]));
 check(cands.every(x=>x>=0&&x<1<<20),'coset label');
 const z=af(pref[i],cands[0]);check(8*z[0]+z[1]===pref[i],'label projection');
 check(!old.some(r=>r[0]===1&&r[1]===pref[i]),'old/new prefix disjoint');
 H.push({c:1,p:pref[i],s,v,cands,group:i,o});
}
check(H.length===32&&new Set(H.map(x=>x.s)).size===32,'all 32 signs');
for(let i=0;i<4;i++){let C=H.filter(x=>x.group===i).flatMap(x=>x.cands);check(C.length===200&&new Set(C).size===200,'common translate disjointness');}
for(let i=0;i<H.length;i++){
 let a=H[i];for(let h of old){let v=h.slice(2).map(x=>den*x), b=h[0]===1?bf(1,h[1],a.p):44,d=b*den*den+dot(a.v,v);check(d<=0||4*d*d<=R2*(88*den*den+norm(v)),'new/old dense Gram');}
 for(let j=0;j<=i;j++){
  let b=H[j], d=pc(a.s^b.s);
  const need=Math.ceil((R2-4*amp*amp*d)/(4*den*den));
  if(i===j)check(need===32,'one long coset stock distances');
  else if(a.group===b.group) check(need<=24&&b.cands.every(x=>!a.cands.includes(x)),'same-prefix signs');
  else if(d===1)check(a.o===b.o&&pc(signBases[a.group]^signBases[b.group])===1&&need===24,'four edge gate covers translate');
  else{
   const prefixLoss= (8-(dot(T[1][a.p>>3],T[1][b.p>>3])+dot(T[1][a.p&7],T[1][b.p&7])))/2;
   check(prefixLoss<=4&&need<=24-prefixLoss,'non-edge inherited distance');
  }
 }
}
check(R2===3180&&4*(2*amp)**2<=2*R2&&4*(amp+2*den)**2<=2*R2,'pure and mixed D93 roots');
for(const w of [10,11])check(4*(4*w*den)**2<=8*w*R2,'retained Golay');
for(const w of [5,6,7,8])check(4*8*(w*den)**2<=4*w*R2,'retained matching');
check(4*(22*den)**2<=28*R2&&4*(16*den)**2<=16*R2&&4*(4*den)**2<=8*R2,'all other sparse stock');
const inherited=46*8**10+153534224+1800, newCount=H.length*W.length*Q.length;
check(inherited===49545659928&&newCount===819200&&inherited+newCount===49546479128,'full ledger');
console.log(JSON.stringify({dimension:93,old:inherited,retainedDense:46*8**10,retainedBodyOnly:153534224,retainedRoots:1800,heads:H.length,cosets:W.length*H.length,cosetSize:Q.length,added:newCount,lowerBound:inherited+newCount,edgePairs:4*W.length**2}));

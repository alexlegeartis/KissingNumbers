
/* 38D fixed planar-QR aperture, rigorous arbitrary head no-append through radius 93/20.
   No native C++; this standalone Node script rebuilds quotient distances in pure JS.
   Expect ~80 MiB scratch (before JS runtime overhead) and large finite traversal.
   Inherits source and sparse geometry in lower57_p48_and_qr_planar_20261126.md. */
"use strict";
const parent={
  "source": "Yang et al., arXiv:2609.35051v1, Appendix B, Tables 11 and 12",
  "source_url": "https://arxiv.org/html/2609.35051v1",
  "license": "CC BY 4.0",
  "length": 39,
  "generators_hex": [
    "4000152727", "2000C0DE0B", "1000AA229D", "8009F5CD6",
    "4004FAE6B", "200ED9AAD", "100BC80CE", "805E4067",
    "40E56DAB", "20B8FB4D", "1096303E", "84B181F",
    "4EFC197", "2BDAD53", "1949B31"
  ],
  "representatives_hex": [
    "0", "ADCFB", "F7795F", "740E29", "940DFF",
    "3B66AB", "6CD301", "4709D4", "98A926", "D23688"
  ]
}
;
const build=(function(parent,marks){
const G=parent.generators_hex.map(x=>BigInt("0x"+x)), B=parent.representatives_hex.map(x=>BigInt("0x"+x));
const split=x=>[Number(x&0xffffffffn)>>>0,Number(x>>32n)>>>0];
const bodyidx=Array.from({length:39},(_,i)=>i).filter(i=>!marks.includes(i));
const pstart=G.length,max=1<<pstart;
const pre=x=>marks.reduce((p,j,i)=>p+(Number(x>>BigInt(j)&1n)<<i),0);
const body=x=>bodyidx.reduce((p,j,i)=>p+((x>>BigInt(j)&1n)<<BigInt(i)),0n);
const H=[0n]; for(let g of G)H.push(...H.slice().map(x=>x^g));
const K=H.filter(x=>pre(x)===0).map(body);
if(K.length!==4096)throw Error("bad shortened rank");
const U=Array.from({length:8},(_,p)=>{let s=marks.map((_,i)=>p>>i&1?-2:2),t=s.reduce((a,b)=>a+b,0);return s.map(x=>3*x-t);});
const H0=H.slice(0).sort((a,b)=>a<b?-1:a>b?1:0);
const P=[];for(let i=0;i<B.length;i++)for(let p=0;p<8;p++){
 if((p===0||p===7)&&i!==0)continue;
 let ph=p===0?4:p===7?3:p; // marks[2]=31
 let x=H0.find(x=>pre(x^B[i])===p);
 P.push({s:body(x^B[i]),ph,i,p});
}
const E=new Map(); for(let x of K){let y=x;while(y){let j=y.toString(2).length-1;if(E.has(j))y^=E.get(j);else{E.set(j,y);break;}}}
const F=Array.from({length:36},(_,i)=>i).filter(i=>!E.has(i));
const elim=[...E.entries()].sort((a,b)=>b[0]-a[0]);
function synd(x){let q=x;for(let [i,g]of elim)if(q>>BigInt(i)&1n)q^=g;return F.reduce((s,j,i)=>s+(Number(q>>BigInt(j)&1n)<<i),0);}
const A=Array.from({length:36},(_,i)=>synd(1n<<BigInt(i)));
function pc(x){x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0F0F0F0F)*0x01010101)>>>24;}
const pairs=K.map(split);
function md(u){let [a,b]=split(u),x=36;for(let [c,d] of pairs)x=Math.min(x,pc(a^c)+pc(b^d));return x;}
const norms=U.map(y=>1296+y.reduce((s,x)=>s+x*x,0));
function needed(u,v){
 const d=U[u].reduce((s,x,j)=>s+x*U[v][j],0);
 for(let i=0;i<=12;i++){let x=1296-72*i+d;if(x<=0||4*x*x<=norms[u]*norms[v])return i;}
 throw Error('undefined');
}
let minimum=40,pairChecks=0,same={};
for(let i=0;i<P.length;i++)for(let j=0;j<i;j++){
 let d=md(P[i].s^P[j].s);let t=needed(P[i].ph,P[j].ph);
 if(d<t)throw Error(JSON.stringify([i,j,d,t,P[i].ph,P[j].ph]));
 minimum=Math.min(minimum,d);pairChecks++;
}
const Kmin=K.filter(x=>x).reduce((s,x)=>{let [a,b]=split(x);return Math.min(s,pc(a)+pc(b));},40);
const ps=P.map(x=>({...x,s:String(x.s),q:synd(x.s)}));
const dist=Array.from({length:8},(_,i)=>Array.from({length:8},(_,j)=>needed(i,j)));
return {K:K.map(String),U,A,F,dist,ps,pairChecks,minimum,Kmin,marks};
});
const h=build(parent,[0,1,31]);
const scan=(function(h){
const Q=1<<24, D=new Uint8Array(Q);D.fill(255); D[0]=0;let t=1,queue=new Uint32Array(Q);
for(let i=0;i<t;i++){
 let x=queue[i],n=D[x]+1;
 for(let j=0;j<h.A.length;j++){let y=x^h.A[j];if(D[y]===255){D[y]=n;queue[t++]=y;}}
}
if(t!==Q)throw Error("rank");
const B=h.U.slice(1,7).map(x=>x.map(z=>20*z)),norm=h=>h.reduce((s,x)=>s+x*x,0),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const n0=518400,cap=77841;
function need(a,b){let c=n0+dot(a,b),A=n0+norm(a),B=n0+norm(b);for(let d=6;d<=11;d++){let x=c-d*28800;if(x<=0||4*x*x<=A*B)return d;}throw Error("angle");}
let menu=new Map();
for(let a=-230;a<=230;a++)for(let b=-230;b<=230;b++){
 let p=[a,b,-a-b];if(norm(p)>cap)continue;
 let r=B.map(x=>need(p,x)),k=r.join();if(!menu.has(k)||norm(menu.get(k).p)<norm(p))menu.set(k,{p,d:r});
}
const ord=[];for(let i=0;i<62;i++)for(let k=1;k<=6;k++){const v=h.ps.filter(p=>p.ph===k)[i];if(v)ord.push(v);}
const pos=[],surv=[];let cut=0,hi=0,relN=0,menuN=0,total=0;const S=new Uint32Array(4**6);
function outer(L,U){
 if(L.some((x,i)=>x>U[i])||L[0]+L[1]+L[2]>0||U[0]+U[1]+U[2]<0)return false;
 for(let i=0;i<27;i++){
  let act=[i%3,(i/3|0)%3,i/9|0],n=0,s=0,f=0;
  for(let j=0;j<3;j++){if(act[j]===1)n++;else{let x=act[j]===0?L[j]:U[j];s+=x;f+=x*x;}}
  if(n===0){if(s===0&&100*f<=2802276)return true;}else if(act.every((x,k)=>x!==1||(-s>=L[k]*n&&-s<=U[k]*n))&&100*(f*n+s*s)<=2802276*n)return true;
 }
 return false;
}
let V=Array.from({length:4**6},(_,i)=>Array.from({length:6},(_,j)=>7+(i>>(2*j)&3)));
const O=V.map(m=>outer([Math.max(-137,575-72*m[0]),Math.max(-137,575-72*m[1]),Math.max(-137,575-72*m[3])],[Math.min(137,72*m[5]-575),Math.min(137,72*m[4]-575),Math.min(137,72*m[2]-575)]));
const F=V.map(m=>[...menu.values()].find(o=>o.d.every((x,i)=>x<=m[i])));
for(let q=0;q<Q;q++){
 let m=[10,10,10,10,10,10],good=true;
 for(let j=0;j<ord.length;j++){
  const u=ord[j],x=D[q^u.q];
  if(x<7){good=false;break;}
  if(x<m[u.ph-1])m[u.ph-1]=x;
 }
 if(!good)continue;hi++;surv.push([q,...m]);
 let i=0;for(let k=0;k<6;k++)i|=(m[k]-7)<<(2*k);S[i]++;
 if(O[i]){relN++;if(pos.length<500)pos.push({q,m,positive:F[i]?.p});}
 if(F[i])menuN++;
}
return {Q,menu:[...menu.values()],hi,relN,menuN,pos,good:Array.from(S).reduce((s,v)=>s+v,0),surv,patterns:V.filter((x,i)=>S[i]&&O[i]).map(m=>[m,S[V.findIndex(x=>x.join()===m.join())]])};
});
function need(x,s){if(!x)throw Error(s);}
// For s≤93/20, actual half-dot threshold at d≤6 is less than 161/2,
// while raw dot≥96−(93/20)(10/3)=161/2.
need(166*155<161*161,"cross distance cannot be below seven");
// Upper E in z=3 head coordinates, E=1/2*sqrt(1392*(1296+||z||²)).
need(721*721*400>=348*(1296*400+77841),"E upper for head radius");
need(3*137*137*100>2*2802276,"coordinate cap after multiply by twelve");
const Q=scan(h);
need(Q.Q===16777216&&Q.hi===880&&Q.good===880&&Q.relN===18&&Q.pos.length===18,"exhaustive outer profile counts");
const a=JSON.stringify(Q.patterns);
need(a==="[[[9,9,8,9,8,8],9],[[8,8,9,8,9,9],9]]","two actual residual profiles");
for(let [v,n]of Q.patterns){
 const allminus=[0,1,3].every(i=>v[i]===8),allplus=[5,4,2].every(i=>v[i]===8);
 need(allminus||allplus,"all three demand signs the same");
}
// A residual profile gives z_i≤b−2 for ALL i, or z_i≥2−b for ALL i.
// Their sum is zero, so b≥2; E≤721 gives b≤25/12. Thus |z_i|≤1/6.
// But the defining norm q(b)=96+48b+(12/29)b²≥192.
need(3<192*36,"impossible residual norm");
const ret={syndromes:Q.Q,coarse_distance_screen:Q.hi,slack_outer:Q.relN,outer_tuples:Q.patterns,exact_append_through_radius_93_over_20:0,scope:"fixed 62 old dense fibres plus one complete affine body fibre, any continuous planar head radius≤93/20; not general continuous exchanges or unrestricted 38D"};
if(typeof report==="function")report(ret);else console.log(JSON.stringify(ret,null,2));


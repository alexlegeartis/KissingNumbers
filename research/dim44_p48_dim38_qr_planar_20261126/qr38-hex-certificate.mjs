
/* 38D planar-hexagon transfer, Lower57. Pure JS/Node exact verifier.
   The Lin octads and retained sparse stars/roots are unchanged from
   KissingNumbers 26998bee... verifications/improved/dim38-qr-octad-switches.
   Optional native witness is injected only as coordinate data, never as a solver premise. */
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
const result={deleted:[],added:[]};
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
const M=build(parent,[0,1,31]),B=M.U.slice(1,7).map(z=>z.map(x=>2*x));
function need(b,s){if(!b)throw Error(s);}
const norm=h=>h.reduce((s,x)=>s+x*x,0),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
function md(s){let z=36;for(const e of M.K){let x=s^BigInt(e),n=0;while(x){x&=x-1n;n++;}if(n<z)z=n;}return z;}
function unpack(q){need(Number.isInteger(q)&&q>=0&&q<2**24,"syndrome");return M.F.reduce((x,j,i)=>x|(BigInt(q>>i&1)<<BigInt(j)),0n);}
need(M.ps.length===62&&M.pairChecks===1891&&M.K.length===4096&&M.Kmin===12,"inherited geometry");
need(new Set(result.deleted).size===result.deleted.length&&result.deleted.every(i=>Number.isInteger(i)&&i>=0&&i<62),"distinct deleted slots");
const T=M.ps.map(p=>[BigInt(p.s),B[p.ph-1]]).filter((_,i)=>!result.deleted.includes(i));
for(let [q,...h] of result.added){
 need(h.length===3&&h.every(Number.isInteger)&&h.reduce((s,x)=>s+x,0)===0&&norm(h)<=384,"safe rational head");
 const q0=unpack(q),n=5184+norm(h);
 for(let [z,k] of T){
  let d=md(q0^z),c=5184-288*d+dot(h,k);
  need(c<=0||4*c*c<=n*(5184+norm(k)),"all new cross pairs");
 }
 T.push([q0,h]);
}
// In coordinate scale with dense magnitudes ±2 and body length36, any head radius r≤√(32/3)<10/3.
// Uniform signed octad bound: 16 (no plane), 14+r (one plane coordinate),
// or 12+√2 r (two plane coordinates); compare √(288+2r²) to each.
need(3*32<=100,"coarse radius cap");
need(16*16<=288&&9*92-9*28*10/3+100>0&&2*100<=9*36,"octad inequalities");
need((4*3+20)**2<9*288&&(8*3+20)**2<9*288,"root and star coarse bounds");
const sparse=395264+2572+288;
const ret={dimension:38,marks:[0,1,31],kernel:[36,12,12],body_head_fibres:T.length,count:sparse+4096*T.length,sparse,deleted:result.deleted,added:result.added,plane:"3 coordinates of sum zero, scaled by 1/6, hence two physical head coordinates"};
if(typeof report==="function")report(ret);else console.log(JSON.stringify(ret,null,2));


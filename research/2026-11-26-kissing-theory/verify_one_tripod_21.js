// Node.js >=16. No dependencies. Bounds enclose real magnitudes with denominator 32768;
// this is not a grid restriction. See the companion proof for justification of ext() and rounding.
"use strict";

const M=32768,P=M*M,H=M/2,sign=(n,i)=>n>>i&1?1:-1;
function rt(n){if(n<=0)return 0;let s=Math.floor(Math.sqrt(n));while((s+1)*(s+1)<=n)s++;while(s*s>n)s--;return s;}
function ac(n){let t=rt(n);return t*t<n?t+1:t;}
const T=Math.ceil(ac(10089*P)/2),mask=Array.from({length:32},(_,m)=>m).filter(m=>(m&17)!=17),vs=[[1,2,2],[2,1,-2],[2,-2,1]],cuts=Array.from({length:32},(_,m)=>vs.flatMap(r=>[1,-1].map(e=>[50,...r.map(t=>24*t*e),49].map((t,i)=>t*sign(m,i)))));
function make(occ){
 const db=Array(32),es=[];
 for(let m of occ)db[m]={l:Array(5).fill(0),h:Array(5).fill(H)};
 for(let m of occ)for(let i=0;i<5;i++)if(db[m^1<<i]){
  db[m].l[i]=rt(P/8);if(m<(m^1<<i))es.push([m,m^1<<i,i]);}
 for(let m of occ)if((m&17)==1||(m&17)==16){const sg=[1,2,3].map(i=>sign(m,i));
 if(vs.some(v=>v.every((t,i)=>Math.sign(t)==sg[i])||v.every((t,i)=>Math.sign(t)==-sg[i]))){let j=(m&1)?4:0;db[m].l[j]=Math.max(db[m].l[j],Math.floor(3*M/10));}}
 return {occ,db,es};
}
function prop(W){
 const {occ,db,es}=W;let any=true;
 for(let round=0;round<100&&any;round++){any=false;
 for(let m of occ){const {l,h}=db[m];let L=l.reduce((s,x)=>s+x*x,0),U=h.reduce((s,x)=>s+x*x,0);if(L>P||U<P)return false;
 for(let i=0;i<5;i++){let x=rt(P-U+h[i]*h[i]),y=ac(P-L+l[i]*l[i]);if(x>l[i]){l[i]=x;any=true;}if(y<h[i]){h[i]=y;any=true;}if(l[i]>h[i])return false;}
 for(let a of cuts[m]){let s=a.reduce((t,x,i)=>t+x*(x>=0?l[i]:h[i]),0);if(s>T)return false;
 for(let i=0;i<5;i++)if(a[i]>0){let z=l[i]+Math.ceil((T-s)/a[i]);if(z<h[i]){h[i]=z;any=true;}}else{let z=h[i]+Math.floor((T-s)/a[i]);if(z>l[i]){l[i]=z;any=true;}}}
 }
 for(let [m,n,i] of es){let a=db[m],b=db[n],ds=0,dx=0;
 for(let j=0;j<5;j++)if(j!=i){ds+=Math.max(a.h[j]-b.l[j],b.h[j]-a.l[j])**2;dx+=a.l[j]*b.l[j];}
 if((a.h[i]+b.h[i])**2+ds<P||dx-a.h[i]*b.h[i]>P/2)return false;let s=rt(P-ds);
 let v=Math.max(s-b.h[i],b.h[i]?Math.floor((dx-P/2)/b.h[i]):0);if(v>a.l[i]){a.l[i]=v;any=true;}
 v=Math.max(s-a.h[i],a.h[i]?Math.floor((dx-P/2)/a.h[i]):0);if(v>b.l[i]){b.l[i]=v;any=true;}
 }
 }return true;
}
function ext(l,h,axes,need){
 const x=axes.map(i=>l[i]),V=[];let z=need-x.reduce((s,t)=>s+t*t,0);
 if(z<=0)return [x];let dd=axes.map(i=>h[i]*h[i]-l[i]*l[i]);if(dd.reduce((s,t)=>s+t,0)<z)return [];
 for(let j=0;j<axes.length;j++)for(let w=0;w<1<<axes.length;w++)if(!(w&1<<j)){let v=z,ar=x.slice();for(let k=0;k<axes.length;k++)if(w&1<<k){v-=dd[k];ar[k]=h[axes[k]];}if(v>=0&&v<=dd[j]){ar[j]=rt(x[j]*x[j]+v);V.push(ar);}}
 return V;
}
function ps(m,n){let df=m^n,pa=[],ne=[];for(let i=0;i<5;i++)((df>>i)&1?ne:pa).push(i);return [pa,ne];}
function cb(m,b,W,ns){let {l,h}=b;for(let i=0;i<5;i++)if(l[i]>h[i])return false;
 if(l.reduce((s,x)=>s+x*x,0)>P||h.reduce((s,x)=>s+x*x,0)<P)return false;
 for(let a of cuts[m]){let pa=[],neg=0,q=P;for(let i=0;i<5;i++)if(a[i]>0)pa.push(i);else{neg+=a[i]*h[i];q-=h[i]*h[i];}
 let shapes=ext(l,h,pa,q);if(!shapes.some(v=>v.reduce((s,t,k)=>s+t*a[pa[k]],neg)<=T))return false;}
 for(let [n,pa,ne] of ns){let c=W.db[n],th=P/2;for(let k of ne)th+=h[k]*c.h[k];let q=P,qq=P;for(let k of ne){q-=h[k]*h[k];qq-=c.h[k]*c.h[k];}
 let A=ext(l,h,pa,q),B=ext(c.l,c.h,pa,qq);let ok=false;for(let x of A){for(let y of B){let s=0;for(let i=0;i<x.length;i++)s+=x[i]*y[i];if(s<=th){ok=true;break;}}if(ok)break;}if(!ok)return false;}
 return true;
}
function cook(W,step,lim=30){
 if(!prop(W))return false;let {occ,db}=W;
 const ns=occ.map(m=>occ.map(n=>[n,...ps(m,n)]).filter(([n,p,v])=>v.length>0&&v.length<=2));
 for(let round=0;round<lim;round++){let changed=false;
 for(let j=0;j<occ.length;j++){let m=occ[j],b=db[m];if(!cb(m,b,W,ns[j]))return false;
 for(let k=0;k<5;k++)for(let d of [-1,1])while(b.h[k]-b.l[k]>=step){let x={l:b.l.slice(),h:b.h.slice()};if(d==-1)x.h[k]=x.l[k]+step;else x.l[k]=x.h[k]-step;if(cb(m,x,W,ns[j]))break;if(d==-1)b.l[k]+=step;else b.h[k]-=step;changed=true;}
 }
 if(!changed)break;if(!prop(W))return false;
 }return true;
}
function copy(W){return {occ:W.occ,es:W.es,db:W.db.map(x=>({l:x.l.slice(),h:x.h.slice()}))};}
function dfs(W,depth=0,lim=18,step=256,stat=[]){
 stat[depth]=(stat[depth]||0)+1;if(!cook(W,step,20))return 0;if(depth>=lim)return 1;
 let ch,best=-1;for(let m of W.occ){let b=W.db[m],deg=W.es.reduce((s,[x,y])=>s+(x==m||y==m),0);for(let i=0;i<5;i++){let dx=b.h[i]-b.l[i];if(dx<256)continue;let r=dx*(1+deg*.12)*((m&17)!=0&&(i==0||i==4)?1.8:1)*(i==0||i==4?1.15:1);if(r>best){best=r;ch=[m,i,Math.floor((b.l[i]+b.h[i])/2)];}}}
 if(!ch)return 1;let [m,i,b]=ch,W2=copy(W);W.db[m].h[i]=b;W2.db[m].l[i]=b;return dfs(W,depth+1,lim,step,stat)||dfs(W2,depth+1,lim,step,stat);
}

const perms=[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]];function sym(m,p,n){let ss=m^2,t=ss&17;for(let j=1;j<=3;j++){let s=!!(ss&(1<<p[j-1]));if(s!=n)t+=1<<j;}return t^2;}function canon(ch){let rep=Number.MAX_SAFE_INTEGER;for(let p of perms)for(let n of [0,1]){let ar=ch.map(m=>sym(m,p,n)).sort((a,b)=>a-b);rep=Math.min(rep,ar[0]*1024+ar[1]*32+ar[2]);}return rep;}function quotient(chs){let M=new Map;for(let ch of chs){let k=canon(ch);if(!M.has(k))M.set(k,ch);}return [...M.values()];}

function vacancies(){let a=[];for(let i=0;i<24;i++)for(let j=0;j<i;j++)for(let k=0;k<j;k++)a.push([mask[i],mask[j],mask[k]]);return quotient(a);}

function run(){
 const reps=vacancies();if(reps.length!==246)throw Error("symmetry");
 let totals=[0,0,0,0],details=[];
 for(const ch of reps){let w=make(mask.filter(m=>!ch.includes(m)));if(!prop(w))continue;totals[0]++;if(!cook(w,400))continue;totals[1]++;if(!cook(w,75))continue;totals[2]++;let tree=[];if(dfs(w,0,22,256,tree))throw Error("UNRESOLVED "+ch);let nodes=tree.reduce((x,y)=>x+y,0);totals[3]+=nodes;details.push([ch,nodes]);}
 if(totals.join()!=[71,17,8,656].join())throw Error("unexpected totals: "+totals.join());
 console.log(JSON.stringify({missingTriples:2024,reps:reps.length,totals,details,status:"all 21-row configurations excluded"}));
}run();

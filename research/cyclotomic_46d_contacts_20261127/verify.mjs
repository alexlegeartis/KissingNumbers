
const n=47,qs=new Set(Array.from({length:23},(_,i)=>((i+1)**2)%n));
const pal=Array.from({length:n},(_,i)=>[...Array.from({length:n},(_,j)=>j===i||qs.has((j-i+n)%n)?1:2),2]);
function basis(V){
 const T=V.map(a=>a.slice());let p=0;
 for(let i=0;i<T[0].length;i++){
  let j=T.findIndex((a,j)=>j>=p&&a[i]);if(j<0)continue;
  [T[p],T[j]]=[T[j],T[p]];if(T[p][i]===2)T[p]=T[p].map(x=>(2*x)%3);
  for(let j=0;j<T.length;j++)if(j!==p&&T[j][i]){const k=T[j][i];for(let l=i;l<T[0].length;l++)T[j][l]=(T[j][l]+6-k*T[p][l])%3;}
  p++;
 }return T.slice(0,p);
}
const C=basis(pal);
if(C.length!==24||!pal.every(a=>pal.every(b=>a.reduce((s,x,i)=>s+x*b[i],0)%3===0)))throw "bad code";

function stats(x){let u=0,v=0;for(let i=0;i<47;i++){u+=x[i];v+=x[i]*x[i];}return [u,47*v-u*u];}
function projected(x,y){let u=0,v=0,p=0,r=0,t=0;for(let i=0;i<47;i++){p+=x[i]*y[i];u+=x[i];v+=y[i];r+=x[i]*x[i];t+=y[i]*y[i];}return [47*p-u*v,47*r-u*u,47*t-v*v];}
function cert(x,y){const [p,r,t]=projected(x,y);return p>0&&4*p*p>r*t;}

const HITS=[[[0,2,0,2,-2,0,2,0,2,2,0,0,2,0,2,0,0,0,2,2,-2,-2,2,2,2,0,-2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,-2,0,2,0,2,4],[0,0,-2,2,0,0,2,0,-2,0,-2,0,0,-2,2,0,-2,0,0,0,-2,-2,0,2,0,0,0,-2,-2,0,-2,0,0,0,-2,0,0,-2,0,0,0,0,-2,-2,2,-2,0,-4]],[[0,2,0,2,-2,0,2,0,2,2,0,0,2,0,2,0,0,0,2,2,-2,-2,2,2,2,0,-2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,-2,0,2,0,2,4],[-1,1,-1,1,-1,-1,1,1,1,1,-1,1,1,1,3,-1,-1,1,1,-1,-1,-3,1,-1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,-1,3,-1,1,-5]],[[0,2,0,2,-2,0,2,0,2,2,0,0,2,0,2,0,0,0,2,2,-2,-2,2,2,2,0,-2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,-2,0,2,0,2,4],[-1,-1,1,1,-1,-1,1,-1,1,1,1,1,1,1,1,-1,-3,-1,1,3,-1,-1,1,1,1,1,-1,-3,-1,1,1,-1,-1,1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,1,5]],[[1,1,1,-1,3,-1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,-3,-1,-1,1,-3,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-5],[1,1,-1,1,3,-1,-1,1,1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,1,1,-1,1,-3,-1,1,-1,-1,-1,-3,1,1,-1,1,-1,5]]];
const CN={"x":[0,2,0,2,-2,0,2,0,2,2,0,0,2,0,2,0,0,0,2,2,-2,-2,2,2,2,0,-2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,-2,0,2,0,2,4],"z":[0,2,-2,2,0,2,2,2,0,0,0,-2,0,2,2,0,0,0,0,0,0,-2,4,2,2,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,4]};
const u=[...new Array(47).fill(1),-5];
function demand(q,m){if(!q)throw Error(m);}
function sq(x){return x.reduce((a,b)=>a+b*b,0);}
function profile(x){
 const a=x.slice(0,47).reduce((u,v)=>u+v,0), inf=x[47];
 return [(a-5*inf)/12,(a+7*inf)/12];
}
function member(x){
 if(x.length!==48)return false;
 const p=((x[0]%2)+2)%2;
 if(x.some(z=>!Number.isInteger(z)||((z-p)%2)))return false;
 const b=x.map((t,i)=>(t-(p?u[i]:0))/2);
 return b.reduce((s,v)=>s+v,0)%2===0 &&
 pal.every(row=>b.reduce((s,v,i)=>s+v*row[i],0)%3===0);
}
const expectedProfiles=[[[0,4],[0,-4]],[[0,4],[2,-3]],[[0,4],[-2,3]],[[2,-3],[-2,3]]];
const expectedNorm=[14,12,6,14], results=[];
for(let i=0;i<HITS.length;i++){
 const [x,y]=HITS[i];
 demand(member(x)&&member(y)&&sq(x)===96&&sq(y)===96,"literal members");
 demand(JSON.stringify([profile(x),profile(y)])===JSON.stringify(expectedProfiles[i]),"literal profiles");
 const d=sq(x.map((z,j)=>z-y[j]))/12;
 demand(d===expectedNorm[i]&&cert(x,y),"true orthogonal projection conflict");
 results.push({profiles:[profile(x),profile(y)],normDifference:d,projectedGram:projected(x,y)});
}
const c=CN.x.map((v,j)=>v-CN.z[j]);
demand(member(CN.x)&&member(CN.z)&&member(c)&&sq(CN.x)===96&&sq(CN.z)===96&&sq(c)===72,"literal center contact");
demand(JSON.stringify([profile(CN.x),profile(CN.z),profile(c)])==="[[0,4],[0,4],[0,0]]","central profile");
demand(cert(CN.x,c)&&CN.x.slice(0,47).every(v=>Math.abs(v)<=2),"flat A to center");
console.log({status:"PASS",results,flatAActualCentralBlocker:projected(CN.x,c)});

const W=[0,2,0,0,2,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,2,0,0,-2,0,0,0,0,-2,0,0,2,0,2,2,0,0,0,0,2,0,0,2,4];
demand(member(W)&&sq(W)===72&&profile(W).join()==="0,4","literal norm-six rail W");
const palettes=new Map;
for(const a of qs)for(let b=0;b<47;b++){
  const r=new Array(48).fill(0);
  for(let j=0;j<47;j++)r[(a*j+b)%47]=W[j];r[47]=W[47];
  demand(member(r)&&sq(r)===72&&profile(r).join()==="0,4","orbit rail");
  palettes.set(r.join(),r);
}
demand(palettes.size===1081,"full norm-six profile orbit");
const z=[...palettes.values()].filter(t=>{
 const d=CN.x.map((v,j)=>v-t[j]);
 return sq(d)===72;
});
demand(z.length===10&&z.every(t=>member(CN.x.map((v,j)=>v-t[j]))&&cert(CN.x,CN.x.map((v,j)=>v-t[j]))),"literal ten strong blockers");
console.log({strongCentralBlockersOnLiteralFlatA:10,fullA6Orbit:palettes.size});

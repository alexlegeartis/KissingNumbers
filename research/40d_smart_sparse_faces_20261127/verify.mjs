// Node 20+: node verify.mjs a40.8.8.4267H
import {readFileSync} from "node:fs";
import {createHash} from "node:crypto";
const input=readFileSync(process.argv[2]);
if(createHash("sha256").update(input).digest("hex")!=="0566ebf4da8b4997c9a782ea537ba3f35125d9c53201fa0cd94551a803ac47e2")throw Error("source SHA-256");
const S=input.toString("ascii").split(/\r?\n/).map(l=>l.trim()).filter(l=>/^[0-9a-fA-F]{10}$/.test(l));
const json=(x)=>JSON.stringify(x);
const need=(v,l)=>{if(!v)throw Error(l)};
const eq=(a,b,l)=>need(json(a)===json(b),l+": "+json(a)+" ≠ "+json(b));
function combos(n,k,f){let X=[];function r(j){if(X.length===k){f(X);return;}for(;j<=n-k+X.length;j++){X.push(j);r(j+1);X.pop();}}r(0);}
const range=n=>Array.from({length:n},(_,i)=>i);
const B=S.map(x=>Array.from(x.match(/../g)).reverse().flatMap((t,j)=>range(8).filter(i=>parseInt(t,16)>>i&1).map(i=>8*j+i)));
need(B.length===4267&&new Set(S).size===4267&&B.every(x=>x.length===8),"source cardinalities");
const Q=range(40).map(n=>range(6).map(k=>{let p=1;for(let i=1;i<=k;i++)p=p*(n-i+1)/i;return Math.round(p);}));
const key=a=>a.reduce((s,x,i)=>s+Q[x][i+1],0);
const F=new Uint16Array(658008), old=new Set(B.map(a=>a.join(",")));
for(const [i,b] of B.entries())combos(8,5,x=>{const k=key(x.map(j=>b[j]));need(!F[k],"two Smart octads intersect above 4");F[k]=i+1;});
const order=(a,b)=>a-b, cov=a=>F[key(a.sort(order))], overlap=(a,b)=>a.filter(x=>b.includes(x)).length;
const per=range(9).map(()=>[0,0,0]), H=[], donors=[], own=range(4267).map(()=>[]), pairs=new Map();
function traverse(a,start,u,v) {
  const k=a.length;
  per[k][u?(v?2:1):0]++;
  if(k===7&&!u)H.push(a.slice());
  if(k===8){if(!old.has(a.join(",")))donors.push([a.slice(),u-1,v-1]);return;}
  for(let i=start;i<40;i++){let p=u,q=v;const off=Q[i][5];let ok=true;
    outer:for(let j=0;j<k;j++)for(let l=j+1;l<k;l++)for(let m=l+1;m<k;m++)for(let n=m+1;n<k;n++){
      const z=F[a[j]+Q[a[l]][2]+Q[a[m]][3]+Q[a[n]][4]+off];
      if(z&&z!==p&&z!==q){if(!p)p=z;else if(!q)q=z;else{ok=false;break outer;}}
    }
    if(ok){a.push(i);traverse(a,i+1,p,q);a.pop();}
  }
}
traverse([],0,0,0);
eq(per,[[1,0,0],[40,0,0],[780,0,0],[9880,0,0],[91390,0,0],[419056,238952,0],[306005,1065808,1246909],[714,86618,277768],[0,4283,648]],"exact owner census");
const hist=[0,0,0];for(const t of H)for(let j=0;j<40;j++)if(!t.includes(j)){
 let common=t.slice(), count=0;
 combos(7,4,u=>{let v=u.map(i=>t[i]);if(cov([...v,j])){count++;common=common.filter(k=>v.includes(k));}});
 need(count>0&&common.length<=2,"7-hole extension ranks");hist[common.length]++;
}
eq(hist,[22074,1416,72],"covered-quintuple tail intersections");
let gaps={7:0,8:0};
for(let block=0;block<5;block++)combos(8,4,T=>{
 const head=T.map(t=>8*block+t),Bs=B.filter(b=>overlap(head,b)===4&&b.filter(i=>(i>>3)===block).length===4);
 let U=new Set();for(const b of Bs)for(const i of b)if((i>>3)!==block){need(!U.has(i),"tail frames disjoint");U.add(i);}
 need(Bs.length===7||Bs.length===8,"group-tetrad incidence");gaps[Bs.length]++;
 if(Bs.length===7){
  const missing=range(40).filter(i=>(i>>3)!==block&&!U.has(i));eq(missing.length,4,"one vacant tetrad");
  for(const i of missing)need(B.some(b=>b.includes(i)&&overlap(head,b)>=3),"vacancy wall at each site");
 }
});
eq(gaps,{7:14,8:336},"Smart five-byte orthoplexes");
for(const [a,u,v] of donors){if(v<0)own[u].push(a);else{const p=u<v?[u,v]:[v,u];const k=p.join(",");if(!pairs.has(k))pairs.set(k,{p,a:[]});pairs.get(k).a.push(a);}}
eq(donors.length,664,"new support inventory");eq(own.filter(x=>x.length===1).length,16,"single owners");
need(!own.some(x=>x.length>1),"one owner never two carriers");
let charged=0,rare=[];
for(const {p:[u,v],a} of pairs.values()){
 const pool=[...a,...own[u],...own[v]];
 need(pool.every((x,i)=>pool.slice(0,i).every(y=>overlap(x,y)>=5)),"all two-owner candidate carriers conflict");
 function faces(b){let L=[];combos(8,5,X=>{let w=X.map(i=>b[i]);L.push(pool.reduce((s,x,i)=>s|(overlap(w,x)===5?1<<i:0),0));});return L;}
 const all=(1<<pool.length)-1, ff=faces(B[u]),gg=faces(B[v]);
 if(ff.some(f=>gg.some(g=>(f|g)===all)))charged++;
 else rare.push([u,v]);
}
eq(pairs.size,252,"two-owner blocks");eq(charged,244,"two fixed 5-face charges");eq(rare.length,8,"explicit unresolved 2-owner signed exceptions");
function det(A){
 if(A.length===1)return A[0][0];
 if(A.length===2)return A[0][0]*A[1][1]-A[0][1]*A[1][0];
 return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}
function raymax(n,M,L){
 const rows=[range(n).map(i=>i?2:1),range(n).map(i=>M.filter(j=>j<=i+1).length)];
 const b=[200,400];if(L!==null){rows.push(range(n).map(i=>i+1));b.push(L);}
 let topN=0,topD=1;
 for(let k=1;k<=b.length;k++)combos(n,k,J=>combos(b.length,k,I=>{
  const A=I.map(i=>J.map(j=>rows[i][j]));let D=det(A);if(!D)return;
  const V=J.map((j,h)=>det(A.map((a,i)=>a.map((v,l)=>h===l?b[I[i]]:v)))*(D>0?1:-1));D=Math.abs(D);
  if(V.some(v=>v<0)||rows.some((r,i)=>V.reduce((s,x,h)=>s+x*r[J[h]],0)>b[i]*D))return;
  const seq=range(n).map(i=>V.reduce((s,v,h)=>s+(J[h]>=i?v:0),0)), N=seq.reduce((s,v)=>s+v*v,0);
  const d=10000*D*D;
  need(Number.isSafeInteger(N)&&Number.isSafeInteger(d),"checked exact number arithmetic");
  if(BigInt(N)*BigInt(topD)>BigInt(topN)*BigInt(d)){topN=N;topD=d;}
 }));
 const gcd=(a,b)=>b?gcd(b,a%b):a,g=gcd(topN,topD);return [topN/g,topD/g];
}
const walls=[[3,4,5,6,7],[3,4,6,7,8]];
eq(walls.map(m=>raymax(11,m,null)),[[196,25],[196,25]],"support eleven radius");
eq(walls.map(m=>raymax(40,m,957)),[[79769,10000],[239747,30000]],"arbitrary support radius with ℓ1 cap 9.57");
console.log(JSON.stringify({sealed:true,source_supports:B.length,owners:per[8],seven_holes:H.length,common_quartet_hist:hist,group_fibres:gaps,charged_two_owner_packets:charged,unresolved_signed_packets:rare,maximum_support11:"196/25",maximum_l1_957_100:"239747/30000"}));

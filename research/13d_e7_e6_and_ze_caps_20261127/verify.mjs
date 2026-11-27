// Node 20+: node verify.mjs
function buildE7E6() {
 const n=7,L=Array.from({length:n},(_,i)=>i), H=[0,1,3], need=(x,l)=>{if(!x)throw Error(l)};
 const ip=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0), neg=v=>v.map(x=>-x), ref=v=>v.find(x=>x)>0?v:neg(v);
 const T=L.map(j=>Array.from(H,i=>(j+i)%7).sort((a,b)=>a-b)), Q=T.map(t=>L.filter(i=>!t.includes(i)));
 const E=L.map(i=>L.map(j=>i===j?2:0));
 for(const q of Q)for(let s=0;s<8;s++){let x=L.map(()=>0);for(let j=0;j<q.length;j++)x[q[j]]=j>0&&s>>(j-1)&1?-1:1; E.push(x);}
 need(E.length===63,"E7 line count");
 for(let i=0;i<E.length;i++)for(let j=0;j<i;j++)need(Math.abs(ip(E[i],E[j]))<=2,"E7 min angle");
 const t=L.map(i=>H.includes(i)?1:0), e=E.filter(x=>ip(x,t)===0);
 need(e.length===36,"E6 root lines");
 function framePartition(V,k){
   const n=V.length, A=V.map((x,i)=>V.flatMap((y,j)=>i<j&&ip(x,y)===0?[j]:[]));
   const frames=[], inc=Array.from({length:n},()=>[]);
   function list(a,W){
     if(a.length===k){const ix=frames.length;frames.push(a.slice());for(const j of a)inc[j].push(ix);return;}
     for(const j of W)list([...a,j],A[j].filter(h=>W.includes(h)));
   }
   list([],Array.from({length:n},(_,i)=>i));
   const mask=frames.map(s=>s.reduce((m,i)=>m|(1n<<BigInt(i)),0n));let ex=null;
   function solve(M,st){
    if(!M){ex=st;return;}
    let opt=null;for(let i=0;i<n;i++)if((M>>BigInt(i))&1n){
      let row=inc[i].filter(j=>(M&mask[j])===mask[j]);if(!row.length)return;
      if(!opt||row.length<opt.length)opt=row;
    }
    for(const j of opt){solve(M^mask[j],[...st,j]);if(ex)return;}
   }
   solve((1n<<BigInt(n))-1n,[]);
   return {numFrames:frames.length,ix:ex?.map(x=>frames[x])??null,F:ex?.map(x=>frames[x].flatMap(i=>[V[i],neg(V[i])]))??null};
 }
 const a=framePartition(E.slice(7),7),b=framePartition(e,4);
 const e7Frames=[E.slice(0,7).flatMap(v=>[v,neg(v)]),...(a.F??[])];
 const E7=E.flatMap(v=>[v,neg(v)]), E6=e.flatMap(v=>[v,neg(v)]);
 const Y=E7.filter(r=>ip(r,t)).map(r=>r.map((x,i)=>3*x-ip(r,t)*t[i]));
 need(Y.length===54&&Y.every(y=>ip(y,y)===24),"E6* shells in rational numerator");
 need(Y.every((a,i)=>Y.slice(0,i).every(b=>ip(a,b)<=12)),"E6* normalized mutual angles");
 return {ip,neg,t,e,E,E6,E7,e7Frames,a,b,Y};
}
function buildPure7(E) {
 const ip=E.ip, N=7,L=Array.from({length:N},(_,i)=>i), H=[0,1,3], F=L.map(j=>Array.from(H,i=>(j+i)%7).sort((a,b)=>a-b));
 const A=[],all=[];
 for(let i=0;i<7;i++)for(let j=i+1;j<7;j++)for(let a of [-1,1])for(let b of [-1,1]){
   let x=L.map(()=>0);x[i]=a;x[j]=b;A.push({x,n:2,t:"D7"});}
 for(const f of F)for(let s=0;s<8;s++){let x=L.map(()=>0);f.forEach((i,j)=>x[i]=s>>j&1?-1:1);A.push({x,n:3,t:"Fano"});}
 for(let i=0;i<7;i++){
   const P=F.filter(f=>f.includes(i)).map(f=>f.filter(j=>j!==i));
   if(P.length!==3||new Set(P.flat()).size!==6)throw Error("Fano star");
   for(let k=0;k<3;k++){const Q=P.filter((_,j)=>j!==k).flat();
     for(let s=0;s<32;s++){let x=L.map(()=>0);x[i]=s&16?-2:2;Q.forEach((i,j)=>x[i]=s>>j&1?-1:1);A.push({x,n:8,t:"Spike",i,k});}
   }
 }
 if(A.length!==812)throw Error("scope");
 for(const z of A){
  if(ip(z.x,z.x)!==z.n)throw Error("norm");
  for(const y of E.E7){
   const d=ip(y,z.x);
   if(d>0&&d*d>2*z.n)throw Error("polar"); // d sqrt(8/n) <= 4.
  }
 }
 const C=Array.from({length:A.length},()=>[]);
 for(let i=0;i<A.length;i++)for(let j=0;j<i;j++){
  const d=ip(A[i].x,A[j].x);
  if(d>0&&4*d*d>A[i].n*A[j].n){C[i].push(j);C[j].push(i);}
 }
 return {A,C};
}
function loadCertificateL80() {
 const n=7, range=n=>Array.from({length:n},(_,i)=>i), Z=range(n), F=Z.map(i=>[i,(i+1)%7,(i+3)%7]), key=a=>a.slice().sort((x,y)=>x-y).join(",");
 const fs=new Set(F.map(key)), perms=[];
 function perm(a,U){
  if(!U.length){perms.push(a);return;}
  for(const x of U)perm([...a,x],U.filter(y=>y!==x));
 }
 perm([],Z);const sym=perms.filter(p=>F.every(f=>fs.has(key(f.map(i=>p[i])))));
 if(sym.length!==168)throw Error("Fano automorphism order");
 const done=new Set(), reps=[];for(const p of perms){
   const k=p.join(",");if(done.has(k))continue;reps.push(p); for(const s of sym)done.add(p.map(i=>s[i]).join(","));
 }
 if(reps.length!==30||done.size!==5040)throw Error("order cover");
 const cg=(a,b)=>b?cg(b,a%b):(a>=0n?a:-a);
 const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0n);
 const norm=a=>dot(a,a);
 const sqrt=(a)=>{let h=1n<<BigInt((a.toString(2).length+1)>>1);for(;;){let j=(h+a/h)>>1n;if(j>=h)return h;h=j;}};
 const sh=(x,d,m)=>x.map(v=>v*m).concat([d*m]);
 function outer(r,cutInput=null,maxSteps=300){
   const inv=Z.map(x=>r.indexOf(x)), f=F.map(s=>s.map(i=>inv[i]).sort((a,b)=>a-b)), Q=f.map(p=>Z.filter(i=>!p.includes(i)));
   const S=f.flatMap(p=>[Z.map(i=>BigInt(p.includes(i)?1:0)),Z.map(i=>BigInt(i===p[2]?-1:p.includes(i)?1:0))]);
   const walls=[{a:Z.map(i=>i===0?1n:0n), b:2n},...Z.map(j=>({a:Z.map(i=>i===j?-1n:i===j+1?1n:0n),b:0n}))];
   let P=range(n+1).map(k=>({x:Z.map(i=>i<k?2n:0n),d:1n,m:0n}));
   for(const v of P)for(const [i,w] of walls.entries())if(dot(v.x,w.a)===v.d*w.b)v.m|=1n<<BigInt(i);
   function clip(w){
    const bit=1n<<BigInt(walls.length);walls.push(w);
    const T=P.map(v=>dot(v.x,w.a)-v.d*w.b);
    if(!T.some(a=>a>0n)){
       for(const [i,v] of P.entries())if(T[i]===0n)v.m|=bit;
       return;
    }
    const K=[];for(let i=0;i<P.length;i++)if(T[i]<=0n)K.push({...P[i],m:P[i].m|(T[i]===0n?bit:0n)});
    const seen=new Set();
    for(let i=0;i<P.length;i++)if(T[i]<0n)for(let j=0;j<P.length;j++)if(T[j]>0n){
       const c=P[i].m&P[j].m;
       if(P.some((v,k)=>k!==i&&k!==j&&(v.m&c)===c))continue;
       const h=P[i],q=P[j], a=T[j],b=-T[i];
       let x=h.x.map((v,k)=>a*v+b*q.x[k]),d=a*h.d+b*q.d,g=x.reduce(cg,d);
       x=x.map(v=>v/g);d/=g;
       const l=x.join(",")+"/"+d;if(seen.has(l))continue;seen.add(l);
       K.push({x,d,m:c|bit});
    }
    P=K;
   }
   for(const Qr of Q)clip({a:Z.map(i=>Qr.includes(i)?1n:0n),b:4n});
   const cuts=[], cNum=5n,cDen=8n;
   function lowcut(v){
     const forms=S.map(s=>dot(v.x,s)), chosen=forms.flatMap((z,k)=>z>0n&&z*z>12n*v.d*v.d?[k]:[]);
     const h=BigInt(chosen.length);
     if(!h||cDen*chosen.reduce((t,j)=>t+forms[j]**2n-12n*v.d*v.d,0n)<=12n*cNum*v.d*v.d)return null;
     // Weights follow the actual rational slopes; reduce by their gcd for concise exact replay.
     let u=chosen.map(j=>forms[j]),g=u.reduce(cg);u=u.map(a=>a/g);
     const d=10000n, target=12n*(h*cDen+cNum)*norm(u)*d*d;
     let b=sqrt(target/cDen);if(b*b*cDen<target)b++;
     const a=Z.map(i=>d*chosen.reduce((z,j,k)=>z+S[j][i]*u[k],0n));
     if(dot(a,v.x)<=b*v.d)throw Error("non-separating tangent");
     return {chosen,u,d,b,a};
   }
   if(cutInput){for(const rec of cutInput){
     const {chosen}=rec,u=rec.u.map(BigInt),d=BigInt(rec.d),b=BigInt(rec.b);
     if(new Set(chosen).size!==chosen.length||chosen.some(j=>j<0||j>=S.length)||u.some(x=>x<=0n)||u.length!==chosen.length||d<=0n||b<=0n||b*b*cDen<12n*(BigInt(chosen.length)*cDen+cNum)*norm(u)*d*d)throw Error("Cauchy tangent not certified");
     clip({a:Z.map(i=>d*chosen.reduce((z,j,k)=>z+S[j][i]*u[k],0n)),b});
    }
   }else for(let i=0;i<maxSteps;i++){
    let v=null;
    for(const k of P)if(norm(k.x)>=8n*k.d*k.d&&(!v||norm(k.x)*v.d*v.d>norm(v.x)*k.d*k.d))v=k;
    if(!v)break;
    const C=lowcut(v);if(!C)return {bad:[v.x.map(String),String(v.d)],r,cuts};
    cuts.push({chosen:C.chosen,u:C.u.map(String),d:String(C.d),b:String(C.b)});
    clip(C);
   }
   let topN=0n,topD=1n;for(const v of P){const x=norm(v.x),q=v.d*v.d;if(x*topD>topN*q){topN=x;topD=q;}}
   const g=cg(topN,topD);return {r,cuts,verts:P.length,max:[String(topN/g),String(topD/g)],closed:topN<8n*topD};
 }
 return {reps,outer};
}
const receipt=[{"r":[0,1,2,3,4,5,6],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":28,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,3,4,6,5],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":36,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,3,5,4,6],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,3,5,6,4],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,3,6,4,5],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"},{"chosen":[0,12],"u":["1","1"],"d":"10000","b":"79373"}],"verts":49,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,3,6,5,4],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"},{"chosen":[0,12],"u":["1","1"],"d":"10000","b":"79373"}],"verts":37,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,4,3,5,6],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":28,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,4,3,6,5],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":32,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,4,5,3,6],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,4,5,6,3],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,4,6,3,5],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":32,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,4,6,5,3],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":28,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,5,3,4,6],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["12602898161","1600000000"],"closed":true},{"r":[0,1,2,5,3,6,4],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["12602898161","1600000000"],"closed":true},{"r":[0,1,2,5,4,3,6],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["12602898161","1600000000"],"closed":true},{"r":[0,1,2,5,4,6,3],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["12602898161","1600000000"],"closed":true},{"r":[0,1,2,5,6,3,4],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["12602898161","1600000000"],"closed":true},{"r":[0,1,2,5,6,4,3],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["12602898161","1600000000"],"closed":true},{"r":[0,1,2,6,3,4,5],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"},{"chosen":[0,12],"u":["1","1"],"d":"10000","b":"79373"}],"verts":43,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,6,3,5,4],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"},{"chosen":[0,12],"u":["1","1"],"d":"10000","b":"79373"}],"verts":37,"max":["787853129","100000000"],"closed":true},{"r":[0,1,2,6,4,3,5],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":32,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,6,4,5,3],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":28,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,6,5,3,4],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,2,6,5,4,3],"cuts":[{"chosen":[0,1],"u":["1","1"],"d":"10000","b":"79373"}],"verts":24,"max":["3151019387","400000000"],"closed":true},{"r":[0,1,3,2,4,5,6],"cuts":[{"chosen":[0],"u":["1"],"d":"10000","b":"44159"},{"chosen":[0,1],"u":["44159","35841"],"d":"10000","b":"3192017054"}],"verts":44,"max":["1027520673009327277","128457728100000000"],"closed":true},{"r":[0,1,3,2,4,6,5],"cuts":[{"chosen":[0],"u":["1"],"d":"10000","b":"44159"},{"chosen":[0,1],"u":["44159","35841"],"d":"10000","b":"3192017054"}],"verts":52,"max":["1027520673009327277","128457728100000000"],"closed":true},{"r":[0,1,3,2,5,4,6],"cuts":[{"chosen":[0],"u":["1"],"d":"10000","b":"44159"},{"chosen":[0,1],"u":["44159","35841"],"d":"10000","b":"3192017054"}],"verts":40,"max":["342355820654151587","42819242700000000"],"closed":true},{"r":[0,1,3,2,5,6,4],"cuts":[{"chosen":[0],"u":["1"],"d":"10000","b":"44159"},{"chosen":[0,1],"u":["44159","35841"],"d":"10000","b":"3192017054"}],"verts":40,"max":["342355820654151587","42819242700000000"],"closed":true},{"r":[0,1,3,2,6,4,5],"cuts":[{"chosen":[0],"u":["1"],"d":"10000","b":"44159"},{"chosen":[0,1],"u":["44159","35841"],"d":"10000","b":"3192017054"},{"chosen":[0,12],"u":["44159","35841"],"d":"10000","b":"3192017054"}],"verts":72,"max":["1027520673009327277","128457728100000000"],"closed":true},{"r":[0,1,3,2,6,5,4],"cuts":[{"chosen":[0],"u":["1"],"d":"10000","b":"44159"},{"chosen":[0,1],"u":["44159","35841"],"d":"10000","b":"3192017054"},{"chosen":[0,12],"u":["44159","35841"],"d":"10000","b":"3192017054"}],"verts":60,"max":["1027520673009327277","128457728100000000"],"closed":true}];
const F=loadCertificateL80(), E=buildE7E6();
const verified=receipt.map((v,i)=>{
 if(JSON.stringify(v.r)!==JSON.stringify(F.reps[i]))throw Error("case mismatch");
 const w=F.outer(v.r,v.cuts);
 if(!w.closed||JSON.stringify(w.max)!==JSON.stringify(v.max)||w.verts!==v.verts)throw Error("wrong final cap certificate");
 return w;
});
const X=buildPure7(E),ip=E.ip;
const D=X.A.filter(v=>v.t==="Fano").map(v=>v.x), d=X.A.filter(v=>v.t!=="Fano").map(v=>v.x.map(x=>v.t==="D7"?2*x:x).join(","));
const Mid=new Set();
for(let i=0;i<D.length;i++)for(let j=0;j<i;j++){
 const q=ip(D[i],D[j]);if(![1,-1,-3].includes(q))throw Error("E7* code");
 if(q===1)Mid.add(D[i].map((x,k)=>x+D[j][k]).join(","));
}
if(D.length!==56||Mid.size!==756||d.length!==756||!d.every(x=>Mid.has(x)))throw Error("literal edge midpoint family");
const zero=Array(7).fill(0);
const list=[];
for(let i=0;i<9;i++)for(const a of E.e7Frames[i])for(const b of E.b.F[i])list.push({x:[...a,...b],n:8});
for(const q of X.A.filter(v=>v.t==="D7"))list.push({x:[...q.x,...zero],n:2});
for(const y of E.Y)list.push({x:[...zero,...y],n:24});
if(list.length!==1146)throw Error("PackingStar replay size");
for(let i=0;i<list.length;i++){
 const a=list[i];if(ip(a.x,a.x)!==a.n)throw Error("replay norm");
 for(let j=0;j<i;j++){const b=list[j],t=ip(a.x,b.x);if(t>0&&4*t*t>a.n*b.n)throw Error("replay cross geometry");}
}
const max=verified.reduce((m,x)=>BigInt(x.max[0])*BigInt(m.max[1])>BigInt(m.max[0])*BigInt(x.max[1])?x:m);
console.log(JSON.stringify({orderedRepresentatives:verified.length,cuts:receipt.reduce((z,r)=>z+r.cuts.length,0),
 maximumSquaredNorm:max.max,strictUnitLoadLower:"5/8",pureE7CapBound:89,edgeAndVertexCatalogueSharpBound:84,
 E7frames:E.e7Frames.length,E6frames:E.b.F.length,fixedMixedStock:1008,exactPackingStarCount:1146,capOnlyFamilyCeiling:1151}));

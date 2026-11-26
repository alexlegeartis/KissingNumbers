import fs from "node:fs";
import crypto from "node:crypto";
import assert from "node:assert/strict";
const bytes = fs.readFileSync(process.argv[2]);
assert.equal(crypto.createHash("sha256").update(bytes).digest("hex"),"5018425664146ee2ace418021537b8b13863fc57c7dc64376ab6f9e811aa5162");
const T = JSON.parse(bytes).template_numerators_and_denominators;
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const pop=x=>{let n=0;while(x){x&=x-1;n++}return n};
assert(T.every(t=>t.slice(0,5).every(v=>Math.abs(v)<=78) && Math.abs(t[5])<=159));
// Products in the Cramer numerators are bounded by 3*315*2*78^2;
// det<=3*78*2*78^2; squared-norm comparison is only on cube vertices,
// so each component in the numerator is then <=abs(det).
const budgets=[0,...T.map((_,i)=>1<<i),...T.flatMap((_,i)=>T.slice(i+1).map((_,k)=>(1<<i)|(1<<(i+1+k))))];
assert.equal(budgets.length,352);
let hist={},vertices=0;
for(let u=0;u<5;u++) for(let v=u+1;v<5;v++) for(const su of [-1,1]) for(const sv of [-1,1]){
  const ii=[0,1,2,3,4].filter(i=>i!==u&&i!==v);
  const A=T.map(t=>ii.map(i=>t[i])), b=T.map(t=>t[5]-su*t[u]-sv*t[v]);
  for(let k=0;k<3;k++){let e=[0,0,0];e[k]=1;A.push(e,e.map(x=>-x),e);b.push(1,1,0)}
  const masks=new Map();
  for(let i=0;i<A.length;i++) for(let j=i+1;j<A.length;j++){
    const ij=cross(A[i],A[j]);
    for(let k=j+1;k<A.length;k++){
      let det=dot(ij,A[k]);
      if(!det)continue;
      const jk=cross(A[j],A[k]),ki=cross(A[k],A[i]);
      let pt=ij.map((x,f)=>b[i]*jk[f]+b[j]*ki[f]+b[k]*x);
      if(det<0){det=-det;pt=pt.map(x=>-x)}
      if(Math.max(...pt.map(Math.abs))>det)continue;
      const L=25*dot(pt,pt),R=49*det*det;
      assert(Number.isSafeInteger(L)&&Number.isSafeInteger(R));
      if(L<R)continue;
      let mask=0,n=0;
      for(let q=0;q<T.length;q++) if(dot(A[q],pt)>b[q]*det){mask|=1<<q;if(++n>2)break}
      if(n>2)continue;
      vertices++;
      let m=0;
      for(let sig=0;sig<8;sig++)if(pt.every((x,i)=>(sig>>i&1)?x<=0:x>=0))m|=1<<sig;
      masks.set(mask,(masks.get(mask)??0)|m);
    }
  }
  let max=0;
  for(const drop of budgets){
    let sig=0;
    for(const [v,w] of masks)if((v&~drop)===0)sig|=w;
    max=Math.max(max,pop(sig));
  }
  hist[max]=(hist[max]??0)+1;
}
assert.deepEqual(hist,{"2":37,"3":3});
assert.equal(vertices,3027);
console.log("PASS: every hot choice and every <=2 paid old rows leaves at most 3 eligible regular octants.");

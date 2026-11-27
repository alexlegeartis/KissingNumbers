// Lower61, 22D A2 density-torus transfer. Fixture-free integer sign-triangle obstruction.
import assert from "node:assert/strict";

const N=24, G=[0]; const base=((1<<11)|(1<<10)|(1<<6)|(1<<5)|(1<<4)|(1<<2)|1);
function pop(x){let n=0;for(let z=x;z;z&=z-1)n++;return n;}
for(let i=0;i<12;i++){let v=base<<i;v^=(pop(v)&1)<<23;G.push(...G.slice().map(x=>x^v));}
const W={};for(let u of G){let n=pop(u);(W[n]??=[]).push(u);}
if(G.length!==4096||new Set(G).size!==4096||W[8].length!==759||W[12].length!==2576||W[16].length!==759)throw Error("Golay input");

const tab=[{"a":-1,"b":-1,"n":20736,"R":44,"key":"10,10"},{"a":1,"b":2,"n":2816,"R":36,"key":"0,6"},{"a":-1,"b":0,"n":20736,"R":44,"key":"8,2"},{"a":1,"b":-1,"n":2816,"R":36,"key":"6,6"},{"a":2,"b":1,"n":2816,"R":36,"key":"6,0"},{"a":2,"b":0,"n":891,"R":32,"key":"8,8"},{"a":2,"b":2,"n":891,"R":32,"key":"4,4"},{"a":0,"b":-1,"n":20736,"R":44,"key":"2,8"},{"a":0,"b":0,"n":49896,"R":48,"key":"0,0"},{"a":0,"b":1,"n":20736,"R":44,"key":"10,4"},{"a":1,"b":0,"n":20736,"R":44,"key":"4,10"},{"a":1,"b":1,"n":20736,"R":44,"key":"2,2"}];
for(const f of tab)f.X=new Int8Array(f.n*22);
const I=new Map(tab.map((f,i)=>[[f.a,f.b].join(),i])),C=new Int32Array(tab.length),x=new Int8Array(24);
function put(){
 let a=(x[0]+x[1])/2,b=(x[0]+x[2])/2,i=I.get([a,b].join());if(i===undefined)return;
 let j=C[i]++*22,t=tab[i].X;t[j]=x[0]-x[1]-x[2];
 for(let k=3;k<24;k++)t[j+k-2]=x[k];
}
for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(let a of [-4,4])for(let b of [-4,4]){
 x[i]=a;x[j]=b;put();x[i]=x[j]=0;
}
for(let u of W[8]){let S=[];for(let i=0;i<24;i++)if(u>>i&1)S.push(i);for(let z=0;z<256;z++)if(!(pop(z)&1)){for(let i=0;i<8;i++)x[S[i]]=z>>i&1?-2:2;put();}for(let i of S)x[i]=0;}
for(let u of G){for(let i=0;i<24;i++)x[i]=u>>i&1?-1:1;for(let i=0;i<24;i++){x[i]*=-3;put();x[i]/=-3;}}
if(C.join()!==tab.map(x=>x.n).join())throw Error("census "+C);
const mod=x=>((x%12)+12)%12, key=(a,b)=>[mod(4*a-2*b),mod(4*b-2*a)].join(),q=(a,b)=>4*a*a-4*a*b+4*b*b,Q=new Map();
for(let a=-8;a<=8;a++)for(let b=-8;b<=8;b++){let k=key(a,b);if(!Q.has(k)||q(a,b)<Q.get(k))Q.set(k,q(a,b));}
function upper(f,g){const r=f.R+g.R-48+Q.get(key(f.a-g.a,f.b-g.b));return r<=0||r*r<=f.R*g.R;}
function dot(x,y,i,j){let s=0;for(let k=1;k<22;k++)s+=x[i+k]*y[j+k];return 3*s+x[i]*y[j];}

const ids=[I.get("-1,-1"),I.get("-1,0"),I.get("1,2"),I.get("1,-1"),I.get("2,1")],P=[];
for(let i of ids){
 const f=tab[i],m=new Map();
 for(let k=0;k<f.n;k++){let q=f.X[k*22];if(!m.has(q))m.set(q,[]);m.get(q).push(k*22);}
 for(let [q,rows]of [...m].sort((a,b)=>a[0]-b[0]))P.push({i,q,r:2*f.R-q*q,n:rows.length,rows});
}
function cross(i,j){
 const x=P[i],y=P[j],f=tab[x.i],g=tab[y.i],wa=Q.get(key(f.a-g.a,f.b-g.b));
 return 2*(f.R+g.R-48+wa-x.q*y.q);
}
function ip(f,g,a,b){let t=0;for(let k=1;k<22;k++)t+=f.X[a+k]*g.X[b+k];return 6*t;}

const TR=[[1,3,12],[4,6,9]],out=[];
function cover(i,j,w){
 const x=P[i],y=P[j],f=tab[x.i],g=tab[y.i],target=cross(i,j),diag=i===j;
 const ax=x.rows.slice(),by=y.rows.slice(),X=new Set(),Y=diag?X:new Set();
 const v=(a,b)=>ip(f,g,a,b)===target;
 const d=[];
 if(diag&&ax.length<2000){
  for(let h=0;h<ax.length&&d.length<210;h++){
    let a=ax[h];if(X.has(a))continue;
    let made=false;
    for(let j=h+1;j<by.length&&!made;j++){
      let b=by[j];if(X.has(b)||!v(a,b))continue;
      for(let k=j+1;k<by.length;k++){
        let c=by[k];if(X.has(c)||!v(a,c)||!v(b,c))continue;
        d.push([a/22,b/22,c/22]);X.add(a);X.add(b);X.add(c);made=true;break;
      }
    }
  }
 }
 let W=d.length*2, tri=d.length;
 for(let a of ax){
   if(X.has(a))continue;
   for(let b of by){
     if(a===b&&diag||Y.has(b)||!v(a,b))continue;
     X.add(a);Y.add(b);d.push([a/22,b/22]);W++;break;
   }
   if(W>=w)break;
 }
 return {i,j,W,tri,target,d};
}
for(let t of TR)for(let k=0;k<t.length;k++)for(let l=k;l<t.length;l++){
 let i=t[k],j=t[l];out.push(cover(i,j,912));
}

assert.equal(ids.map(i=>tab[i].n).join(),"20736,20736,2816,2816,2816");
assert.equal(ids.reduce((s,i)=>s+tab[i].n,0),49920);
assert.equal(P.length,17);
assert.equal(new Set(TR.flat()).size,6);
for(const p of P)for(const a of p.rows){
 const f=tab[p.i];let s=0;for(let k=1;k<22;k++)s+=f.X[a+k]**2;
 assert.equal(3*s,p.r);assert.equal(f.X[a],p.q);
}
for(const {i,j,W,target,d} of out){
 const f=tab[P[i].i],g=tab[P[j].i],diag=i===j;
 assert.equal(target,diag?2*P[i].r-96:96);
 let X=new Set(),Y=diag?X:new Set(),w=0;
 for(const edge of d){
  assert(edge.length===2||diag&&edge.length===3);
  const a=edge[0]*22,b=edge[1]*22;
  assert(P[i].rows.includes(a)&&P[j].rows.includes(b));
  assert(a!==b||!diag);
  assert(!X.has(a)&&!Y.has(b));X.add(a);Y.add(b);
  assert.equal(ip(f,g,a,b),target);
  w++;
  if(edge.length===3){
   const c=edge[2]*22;assert(P[j].rows.includes(c));assert(!Y.has(c));Y.add(c);
   assert.equal(ip(f,g,a,c),target);assert.equal(ip(f,g,b,c),target);w++;
  }
 }
 assert.equal(W,w);assert(W>=912);
}
let n=49920+(891+891)+2-912*2;
assert.equal(n,49880);
console.log(JSON.stringify({sourceKinds:P.length,triangles:TR,verifiedCovers:out.map(({d,...r})=>r),
  rawCandidateCount:49920,with891BanksAndTwoPoles:n},null,2));

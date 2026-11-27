// Lower61, 22D fixed two nearest A2 projection banks, direct exact residual check.
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
const A=tab[I.get("-1,-1")],B=tab[I.get("-1,0")];
if(!upper(A,B))throw Error("incompatible baseline");
let calls=0;
function survive(f){
 let qs=[A,B].filter(a=>!upper(f,a)),Y=[];
 if(!qs.length)return {c:f.n,all:true,Y:null};
 for(let i=0;i<f.n;i++){
  let fail=false;
  scan:for(let bank of qs){
   let H=0;while((H+1)*(H+1)<=bank.R*f.R)H++;
   for(let j=0;j<bank.n;j++){
    calls++;
    const v=dot(f.X,bank.X,i*22,j*22);
    if(v>H){if(v*v<=bank.R*f.R)throw Error("threshold");fail=true;break scan;}
   }
  }
  if(!fail)Y.push(i);
 }
 return {c:Y.length,all:false,Y};
}
const res=tab.map(f=>f===A||f===B?{c:f.n,anchor:true}:survive(f));

assert.equal(res.map(x=>x.c).join(),"20736,0,20736,0,2816,0,0,0,0,0,0,0");
assert.equal(res.reduce((z,r)=>z+r.c,0),44288);
for(let i=0;i<tab.length;i++){
 const f=tab[i];for(let j=0;j<f.n;j++)assert.equal(dot(f.X,f.X,22*j,22*j),2*f.R);
}
console.log(JSON.stringify({exactDotComparisons:calls,
  fixedCompleteAnchorUpper:44288,survivingFibreCounts:res.map(r=>r.c)},null,2));

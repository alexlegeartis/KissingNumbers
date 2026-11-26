import {readFileSync} from "node:fs";

const H=x=>{x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0xf0f0f0f)*0x1010101)>>>24};
const par=x=>H(x)&1, span=gs=>gs.reduce((a,g)=>a.concat(a.map(x=>x^g)),[0]);
const ebl=Array.from({length:16},(_,i)=>i);
const code=f=>ebl.reduce((a,i)=>a|(f(i)<<i),0);
const lin=[0,1,2,3].map(i=>code(x=>x>>i&1)), mon=[];for(let i=0;i<4;i++)for(let j=0;j<i;j++)mon.push(lin[i]&lin[j]);
const aff=span([65535,...lin]), reps=span(mon), rm2=reps.flatMap(p=>aff.map(a=>p^a)), oct=aff.filter(a=>H(a)==8);
const bmap=reps.map((v,i)=>aff.some(a=>H(a^v)==6)), bis=(a,b)=>par(reps[a]&reps[b]);
const s6=rm2.filter(w=>H(w)==6);
const h0=0,h1=1,hot=(1<<h0)|(1<<h1),LL=aff.filter(x=>(x&hot)==0);
const q= bmap.findIndex(Boolean), S0=aff.map(x=>x^reps[q]).filter(x=>H(x)==6);
const E=reps.map((r,i)=>i).filter(i=>i>0&&bmap[i]&&!bis(i,q));
let F=null;const search=(a,b)=>{if(a.length==6){F=a;return true}for(let i=0;i<b.length;i++)if(search(a.concat([b[i]]),b.slice(i+1).filter(j=>bmap[j^b[i]])))return true;return false;};
search([0],E);if(F===null)throw Error('empty 6');
const Dens=F.flatMap(q=>aff.map(x=>x^reps[q])), flatSet=new Set(Dens);
const w=rm2.find(x=>!(x&hot)&&Dens.every(s=>!flatSet.has(s^x)));
if(w===undefined)throw Error('cannot shift');
const can=(s,m=65535)=>Math.min(...LL.map(l=>s^(l&m)));
const pool=[];
function insert(t,m,s,j,seed){pool.push([t,m,s,j,seed?1:0]);}
for(let m of s6){let bits=ebl.filter(i=>m>>i&1);for(let z=0;z<64;z++){let s=bits.reduce((b,i,k)=>b|((z>>(k)&1)<<i),0);if(s!=can(s,m))continue;if((m&hot)==hot&&!(s&hot))continue;let T=S0.includes(m)&&par(s);insert(0,m,s,-1,T);}}
for(let j of ebl)for(let s of rm2){if(s!=can(s))continue;if((j==h0||j==h1)&&!(s&hot))continue;insert(1,65535,s,j,false);}
for(let s of rm2)if(s==can(s))insert(2,65535,s,-1,flatSet.has(s));
const flip=(p)=>p[2]^w;
const kern=(a,b,l=0)=>{
 let d=a[2]^b[2]^l, u=a[0],v=b[0],m=a[1]&b[1];
 if(u>v)throw Error('sort types');
 let sd=H(m)-2*H(d&m);
 let sign=j=>1-2*(d>>j&1);
 if(u==0&&v==0)return sd>2;
 if(u==0&&v==1)return sd+((a[1]>>b[3]&1)?2*sign(b[3]):0)>4;
 if(u==0&&v==2)return sd>4;
 if(u==1&&v==1)return sd+(a[3]==b[3]?8*sign(a[3]):2*(sign(a[3])+sign(b[3])))>8;
 if(u==1&&v==2)return sd+2*sign(a[3])>8;
 return sd>7;
};
function conflict(a,b){
 if(a[0]>b[0])[a,b]=[b,a];
 if(b[0]==2&&a[0]>=1&&can(a[2]^b[2])==can(w))return true;
 return LL.some(l=>kern(a,b,l));
}
const seeds=pool.flatMap((p,i)=>p[4]?[i]:[]);
for(let i=0;i<seeds.length;i++)for(let j=0;j<i;j++)if(conflict(pool[seeds[i]],pool[seeds[j]]))throw Error('seed clash');

function buildStatus(indices,shift,old=false){
 let V=[];const add=(x,y=0,r=0)=>V.push({x:Int16Array.from(x),y,r:typeof r=='number'?null:Int16Array.from(r)});
 for(let i=0;i<16;i++)for(let j=i+1;j<16;j++)for(let a of [-12,12])for(let b of [-12,12]){if(!old&&i==h0&&j==h1&&a==12&&b==12)continue;let x=Array(16).fill(0);x[i]=a;x[j]=b;add(x);}
 for(let o of oct){let p=ebl.filter(j=>o>>j&1);for(let s=0;s<256;s++)if(par(s)){
   let x=Array(16).fill(0);p.forEach((j,k)=>{x[j]=s>>k&1?-6:6});add(x);
 }}
 for(let s of [-1,1])add(Array(16).fill(0),12*s);
 for(let side of [1,-1]){
  let ws=side>0?0:shift;
  if(old){
   for(let m of S0){for(let s of aff){let p=0;while(!(m>>p&1))p++;let word=(s^(1<<p)^ws)&m;
    add(ebl.map(j=>(m>>j&1)?(word>>j&1?-6:6):0),side*6);
   }}
   for(let s of Dens){let z=s^ws;add(ebl.map(j=>z>>j&1?-4:4),side*4);}
  }else{
   for(let i of indices){const [t,m,s,j]=pool[i];for(let g of LL){let z=s^g^ws;let x=ebl.map(k=>!(m>>k&1)?0:(t==0?6:t==2?4:k==j?9:3)*(z>>k&1?-1:1));add(x,(t==2?4:6)*side);}}
  }
 }
 if(!old)for(let s of [-1,1])add(Array(16).fill(0),s*6,ebl.map(j=>j==h0||j==h1?6:0));
 return V;
}
const sgn=(x,y)=>{if(!y)return Math.sign(x);if(!x)return Math.sign(y);if((x>0)===(y>0))return Math.sign(x);return x*x>3*y*y?Math.sign(x):Math.sign(y);};
function rr(a,b){let x=2*a.y*b.y,y=0;if(!a.r&&!b.r){for(let j=0;j<16;j++)x+=a.x[j]*b.x[j];return [x,0];}
 for(let j=0;j<16;j++){let aa=a.r?a.r[j]:0,bb=b.r?b.r[j]:0;x+=a.x[j]*b.x[j]+3*aa*bb;y+=aa*b.x[j]+bb*a.x[j];}
 return [x,y];
}
function verifyCode(V,n){if(V.length!=n)throw Error('count');
 let eq=0;
 for(let i=0;i<V.length;i++){
   let a=V[i],p=rr(a,a);if(p[0]!=288||p[1]!=0)throw Error('norm '+[i,p]);
   for(let j=0;j<i;j++){let [x,y]=rr(a,V[j]);let s=sgn(x-144,y);if(s>0)throw Error('pair '+[i,j,x,y]);if(s==0)eq++;}
 }return {n,pairs:n*(n-1)/2,contacts:eq};
}

if(pool.length!==7696||seeds.length!==86||w!==34952)throw Error('unexpected generation');
if(process.argv[2]){
 let q=readFileSync(process.argv[2],'utf8').trim().split(/\s+/).map(Number),[m,s,...rows]=q;
 if(m<88||rows.length!==m||!rm2.includes(s)||(s&hot)||new Set(rows).size!==rows.length||rows.some(i=>!pool[i]))throw Error('invalid positive witness input');
 console.log(JSON.stringify(verifyCode(buildStatus(rows.slice(0,88),s),5731)));
}else{
 console.log(JSON.stringify(verifyCode(buildStatus([],w,true),5730)));
 console.log(JSON.stringify(verifyCode(buildStatus(seeds,w),5699)));
}


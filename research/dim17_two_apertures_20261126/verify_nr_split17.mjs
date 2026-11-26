
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

const ds=[7440,7441,7442,7443,7492,7493,7494,7495,7532,7533,7534,7535,7548,7549,7550,7551,7572,7573,7574,7575,7624,7625,7626,7627,7640,7641,7642,7643,7680,7681,7682,7683], D=[[0,1,3,2],[4,5,7,6],[16,17,19,18],[20,21,23,22],[24,25,27,26],[28,29,31,30],[10,11,9,8],[14,15,13,12]];
let C=[];for(let i=0;i<pool.length;i++)if(pool[i][0]!=2){let z=0;for(let j=0;j<32;j++)if(LL.some(g=>kern(pool[i],pool[ds[j]],g)))z=(z|(1<<j))>>>0;
 if(H(z)<=8)C.push([i,z]);}
let N=C.length, GG=new Uint8Array(N*N);for(let i=0;i<N;i++)for(let j=0;j<i;j++){let a=pool[C[j][0]],b=pool[C[i][0]],r=LL.some(g=>kern(a,b,g));if(r)GG[i*N+j]=GG[j*N+i]=1;}
let cert={},max=0,fail=[],top=[];function cv(V,mode){let S=V.slice();if(mode==1)S.reverse();if(mode>=2){let dd=new Int16Array(N);for(let a of V){let k=0,z=a*N;for(let b of V)k+=GG[z+b];dd[a]=k;}S.sort((a,b)=>(mode==2?dd[b]-dd[a]:dd[a]-dd[b])||a-b);}
 let cl=[];for(let a of S){let g=cl.find(k=>k.every(b=>GG[a*N+b]));if(g)g.push(a);else cl.push([a])}return cl.length;}
function walk(r,d,counts){if(r==8){let V=[];for(let i=0;i<N;i++)if(!(C[i][1]&~d))V.push(i);let x=Math.min(cv(V,0),cv(V,1));if(x>=64)x=Math.min(x,cv(V,2),cv(V,3));max=Math.max(max,x);if(x>=64)fail.push([d>>>0,x,V.length]);cert[x]=(cert[x]||0)+1;
 return;}for(let j=0;j<4;j++)if(counts[j]<2){counts[j]++;walk(r+1,d|(1<<D[r][j]),counts);counts[j]--;}}
walk(0,0,[0,0,0,0]);if(N!==3968||max!==38||fail.length||Object.values(cert).reduce((a,b)=>a+b,0)!==2520)throw Error('cover mismatch');console.log(JSON.stringify({N,cert,max,nfail:fail.length}));

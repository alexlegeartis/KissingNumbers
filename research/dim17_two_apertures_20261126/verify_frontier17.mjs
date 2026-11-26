
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

function bad0(a,b){if(a[0]>b[0])[a,b]=[b,a];return LL.some(l=>kern(a,b,l))}
const mapp=new Map(), pbits=[], Rsz=11;
for(let i=0;i<pool.length;i++)if(!pool[i][4]){let a=0n,c=0;for(let j=0;j<seeds.length;j++)if(bad0(pool[i],pool[seeds[j]])){a|=1n<<BigInt(j);c++;}if(c>Rsz)continue;
 if(!mapp.has(a))mapp.set(a,[0,c]);mapp.get(a)[0]++;pbits.push(a.toString());}
const K=[...mapp.keys()], S=new Set(K), Q=K.slice(), count={},best={};
function bc(b){let n=0;while(b){b&=b-1n;n++;}return n}
let all=0, bad=0;
for(let i=0;i<Q.length;i++){
 const a=Q[i],z=bc(a);let c=0;for(let k of K){if((k&a)==k){c+=mapp.get(k)[0];continue;}if(z>=Rsz)continue;const b=k|a;if(!S.has(b)&&bc(b)<=Rsz){S.add(b);Q.push(b);}}
 count[z]=(count[z]||0)+1;best[z]=Math.max(best[z]??-1111,c-z);if(c>=z+2)bad++;
}
console.log(JSON.stringify({bots:pbits.length,init:K.length,count,best,total:Q.length,bad}));
const hi=Array(100).fill(0);
for(let m=0;m<65536;m++)if(H(m)==5&&oct.every(o=>H(o&m)<=4)){let k=S0.reduce((z,o)=>z+[0,0,8,4,10,6][H(o&m)],0);hi[k]++;}
if(Q.length!==27363||bad||hi[60]!==1152||hi[82]!==1440||hi[86]!==96||hi.reduce((s,x)=>s+x,0)!==2688)throw Error('unexpected exact certificate');
console.log(JSON.stringify({high5Supports:2688,high5IndividualCandidates:86016,fullCost:[...hi.entries()].filter(r=>r[1]),scope:"exterior matching and at most eleven old rank-three orbits"}));


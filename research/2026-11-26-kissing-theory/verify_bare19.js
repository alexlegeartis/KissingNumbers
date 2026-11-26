
const R=[[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],P=[];for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)P.push([i,j]);const V=P.flatMap(([i,j])=>[1,-1].map(s=>R.map(row=>s*(row[i]-row[j])/2)));const parities=Array.from({length:8},(_,c)=>P.reduce((a,[i,j],k)=>a|(3*(1-((c>>i^c>>j)&1))<<(2*k)),0));const cf=[0,1,1,0,1,0,0,1];

function boxes(q){
 let E=[];
 for(let i=-4*q;i<=4*q;i+=2)for(let j=-4*q;j<=4*q;j+=2)for(let k=-4*q;k<=4*q;k+=2){
 let t=[i,j,k],u=t.map(x=>Math.abs(x)+1),z=t.map(x=>Math.max(Math.abs(x)-1,0)).sort((a,b)=>b-a);
 let r=u[0]**2+u[1]**2+u[2]**2,N=16*q*q+r,B=2*N;const big=d=>d>0&&d*d>B;
 if(r<2*q*q||big(2*(z[0]+z[1]))||big(2*q+2*z[0]))continue;
 let A=0,ok=true;for(let h=0;h<12;h++){let v=V[h],z=v[0]*i+v[1]*j+v[2]*k-2;if(big(4*q+z)){ok=false;break}if(big(6*q+z))A|=1<<h}for(let z=0;z<6;z++){let a=2*(z%2?-1:1)*t[z>>1]-2;if(big(4*q+a))A|=1<<(12+z)}if(ok)E.push({t,r,A,diag:true,N,
 i,j,k, abs:Math.abs(i)+Math.abs(j)+Math.abs(k)})
 }return E
}
function lb(a,b){let d=8*q0*q0+a.i*b.i+a.j*b.j+a.k*b.k-a.abs-b.abs-3;return d<=0||4*d*d<=a.N*b.N}
function boxstates(B){
 let sets=[new Set(),new Set()],sing=[new Set(),new Set()],Gs=new Map();
 for(let b of B){let k=(b.A<<1)+(b.diag?1:0);if(!Gs.has(k))Gs.set(k,[]);Gs.get(k).push(b);sing[0].add(b.A);if(b.diag)sing[1].add(b.A)}
 let gg=[...Gs.keys()],pairtest=0,miss=0;
 for(let i=0;i<gg.length;i++)for(let j=0;j<=i;j++){let a=gg[i],b=gg[j],s=a>>1|b>>1,odd=a&b&1;if(sets[odd].has(s))continue;
 let L=Gs.get(a),M=Gs.get(b),ok=false;
 out:for(let v of L)for(let w of M){pairtest++;if(lb(v,w)){ok=true;break out}}
 if(ok){sets[0].add(s);if(odd)sets[1].add(s)}else miss++
 }
 const table=([s0,s1])=>[s0,s1].map(S=>{let t=new Uint8Array(262144);for(let s of S)t[s]=1;for(let b=0;b<18;b++)for(let m=0;m<262144;m++)if(m>>b&1)t[m]|=t[m^(1<<b)];return t});
 const sTab=table(sing),pTab=table(sets);
 let s=Array.from({length:262144},(_,m)=>{let score=0,flags=[];for(let c=0;c<8;c++){let g=262143^(m^(parities[c]|((cf[c]?0:63)<<12))),lx=pTab[cf[c]][g]?2:sTab[cf[c]][g]?1:0;score+=lx;flags.push(lx)};return [score,m,flags]});
 return {b:B.length,pairtest,miss,h:s.reduce((h,[a])=>(h[a]=(h[a]||0)+1,h),{}),best:s.reduce((a,b)=>Math.max(a,b[0]),0),winner:s.filter(x=>x[0]>14).slice(0,2)}
}

const q0=process.argv[2]?Number(process.argv[2]):12;
if(q0!==12&&q0!==16)throw Error('choose q=12 or q=16');
const r=boxstates(boxes(q0));
const want={2:512,4:9728,5:12288,6:70656,7:35328,8:99240,9:24048,10:8664,11:800,12:600,13:240,14:40};
if(r.best!==14||r.winner.length||JSON.stringify(r.h)!==JSON.stringify(want))throw Error(JSON.stringify(r));
console.log(r);

function abs(x){return x<0n?-x:x}function gcd(a,b){for(a=abs(a),b=abs(b);b;[a,b]=[b,a%b]);return a||1n;}
function Q(a,b=1n){a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero");if(b<0n){a=-a;b=-b;}let d=gcd(a,b);return [a/d,b/d];}
function add(a,b){return Q(a[0]*b[1]+a[1]*b[0],a[1]*b[1]);}
function sub(a,b){return Q(a[0]*b[1]-a[1]*b[0],a[1]*b[1]);}
function mul(a,b){return Q(a[0]*b[0],a[1]*b[1]);}
function div(a,b){return Q(a[0]*b[1],a[1]*b[0]);}
function neg(a){return [-a[0],a[1]];}
function cmp(a,b){return a[0]*b[1]-a[1]*b[0];}
function ss(a){return a[1]===1n?a[0].toString():a.join("/");}
function ls(M){let rows=M.length,n=M[0].length,k=0,piv=[];
 for(let j=0;j<n;j++){let i=k;while(i<rows&&!M[i][j][0])i++;if(i===rows)continue;[M[i],M[k]]=[M[k],M[i]];let d=M[k][j];for(let l=j;l<n;l++)M[k][l]=div(M[k][l],d);for(let h=0;h<rows;h++)if(h!==k&&M[h][j][0]){let a=M[h][j];for(let l=j;l<n;l++)M[h][l]=sub(M[h][l],mul(a,M[k][l]));}piv.push(j);k++;if(k===rows)break;}
 return {M,piv};}
const Z=Q(0),O=Q(1),L=8;
function arr(f){return Array.from({length:L},(_,i)=>f(i));}
function sig(n,k){let z=0n;for(let i=1;i<=n;i++)if(n%i===0)z+=BigInt(i)**BigInt(k);return z;}
const E2=arr(i=>Q(i===0?1n:-24n*sig(i,1))),E4=arr(i=>Q(i===0?1n:240n*sig(i,3)));
const wide=(a)=>arr(i=>i%2?Z:a[i/2]),su=(a,b)=>a.map((x,i)=>sub(x,b[i])),sc=(a,c)=>a.map(x=>mul(x,c));
const conv=(a,b)=>arr(i=>a.slice(0,i+1).reduce((r,x,j)=>add(r,mul(x,b[i-j])),Z));
function pow(a,e){let r=arr(i=>i?Z:O);for(let j=0;j<e;j++)r=conv(r,a);return r;}
const A=su(sc(wide(E2),Q(2)),E2),AS=su(sc(E2,Q(1,2)),wide(E2)),B=E4,BS=wide(E4);
const E=pow(arr(i=>{if(!i)return O;let j=Math.sqrt(i);return Number.isInteger(j)?Q(j%2?-2:2):Z;}),4);
const ES=sc(pow(arr(i=>{let j=(Math.sqrt(8*i+1)-1)/2;return Number.isInteger(j)?Q(2):Z;}),4),Q(-1,4)); // factor Q^{1/2} out
function basis(odd){
 let wo=odd?11:12,cols=(wo>>1)+1;return [Array.from({length:cols},(_,j)=>conv(odd?E:pow(E,0),conv(pow(A,wo-2*j),pow(B,j)))),Array.from({length:cols},(_,j)=>conv(odd?ES:pow(ES,0),conv(pow(AS,wo-2*j),pow(BS,j))))];
}
function fix(odd,cusp,known=[]){
 let [U,V]=basis(odd),M=[];for(let i=0;i<3;i++)M.push([...U.map(u=>u[i]),i===0?O:Z]);for(let i=0;i<cusp;i++)M.push([...V.map(v=>v[i]),Z]);for(const [i,val] of known)M.push([...V.map(v=>v[i]),Q(val)]);
 let {piv}=ls(M);if(piv.at(-1)===U.length)return {no:true,piv};
 if(piv.length!==U.length)return {under:true,m:M.map(row=>row.map(ss)),piv};
 let v=M.map(row=>row.at(-1)).slice(0,U.length);
 return {v:v.map(ss),co0:arr(i=>ss(v.reduce((r,x,j)=>add(r,mul(x,V[j][i])),Z))),coInf:arr(i=>ss(v.reduce((r,x,j)=>add(r,mul(x,U[j][i])),Z))),coef3:ss(v.reduce((r,x,j)=>add(r,mul(x,U[j][3])),Z))};
}

function req(x,msg){if(!x)throw Error(msg);}
function eq(x,y,msg){req(cmp(x,Q(y))===0n,msg||ss(x)+' != '+y);}
function ast(x,y,msg){req(cmp(x,y)<0n,msg||ss(x)+' >= '+ss(y));}
function qa(x){return Q(abs(x[0]),x[1]);}
function qsum(A){return A.reduce(add,Q(0));}
function row(odd,cusp,known=[]){
  let f=fix(odd,cusp,known);
  req(!f.no&&!f.under,'The cusp linear system must have a unique solution');
  return f;
}
let data=[row(true,1,[[1,2],[2,0]]),
  row(false,2,[[2,2],[3,0]]),row(true,2,[[2,2]]),
  row(false,3,[[3,2]]),row(true,3),row(false,4)];

let s=[6,8,10,12,14,16];
const expected=[
 ["1833472","283100832","8069621760","58859360256"],
 ["380416","24422048","86341632","-1925148672"],
 ["64000","509600","-17782784","80384000"],
 ["-4608","108192","-1443840","30445568"],
 ["-1536","-14688","829440","-28667904"],
 ["-12800","304800","-3311616","37785600"]];
for(let i=0;i<data.length;i++)req(JSON.stringify(data[i].coInf.slice(3,7))===JSON.stringify(expected[i]),'type '+s[i]);
req(data[4].co0[3]==='108' && data[5].co0[4]==='2400','nonzero leading terms of deep cosets');
const zero12=row(false,3,[[3,0]]),unit12=row(false,3,[[3,1]]);
req(JSON.stringify(zero12.coInf)===JSON.stringify(data[5].coInf),'constant term');
let affine12=unit12.coInf.slice(3,7).map((v,i)=>BigInt(v)-BigInt(zero12.coInf[3+i]));
req(affine12.join(',')==='4096,-98304,933888,-3670016','affine norm-12 type');
function delta(){
  let r=arr(i=>i===1?O:Z);
  for(let k=1;k<L;k++){const p=arr(i=>i===0?O:i===k?Q(-1):Z);r=conv(r,pow(p,24));}
  return r;
}
let Del=delta(),Theta=arr(i=>add(sub(pow(E4,6)[i],mul(Q(1440),conv(pow(E4,3),Del)[i])),mul(Q(125280),pow(Del,2)[i])));
req(Theta.slice(0,7).map(ss).join(',')==='1,0,0,52416000,39007332000,6609020221440,437824977408000','ordinary theta');
let [N6,N8,N10,N12]=Theta.slice(3,7),TARGET=Q(52416000),gamma8=Q(12800),gamma10=Q(26700,13),gamma12=Q(4352,5);
const H10=(x,y)=>add(x,mul(Q(11,312),y));
const H12=(x,y,z)=>add(add(x,mul(Q(8,135),y)),mul(Q(1,540),z));
eq(H10(N6,N8),1427674500n);eq(H12(N6,N8,N10),14602887936n);
let M8=div(add(N6,gamma8),gamma8),M10=div(add(H10(N6,N8),gamma10),gamma10),M12=div(add(H12(N6,N8,N10),gamma12),gamma12);
eq(M8,4096n);eq(M12,2n**24n);req(ss(M10)==='61865984/89','minimum index 10');
for(let i=0;i<s.length;i++){
 const [h,g,j,k]=expected[i].map(x=>Q(x)),eps=s[i]%4===0?1:-1;
 const e8=add(h,gamma8),e10=add(H10(h,g),gamma10),e12=add(H12(h,g,j),gamma12);
 req(e8[0]>=0 && e10[0]>=0 && e12[0]>=0,'spectral lower bounds');
 req(cmp(qa(sub(g,Q(304800))),mul(Q(154),e8))<=0n,'norm-8 weight');
 req(cmp(qa(sub(j,Q(-1241088-2070528*eps))),mul(Q(1482),e10))<=0n,'norm-10 weight');
 const a12=s[i]===12?2:0;
 const err=add(sub(k,Q(4558848+33226752*eps)),Q(3670016*a12));
 req(cmp(qa(err),mul(Q(1755),e12))<=0n,'norm-12 weight');
 req(ss(err)===(s[i]===8?'-1962934272':s[i]===6?'58888028160':s[i]===10?'109051904':'0'),'norm-12 signed error');
}
let [ah,ag,aj,ak]=affine12.map(x=>Q(x));
req(ah[0]>0n && H10(ah,ag)[0]>0n,'nonnegative affine spectral slopes');
eq(H12(ah,ag,aj),0);
req(cmp(qa(ag),mul(Q(154),ah))<=0n,'linear type norm8');
eq(div(aj,H10(ah,ag)),1482);
eq(add(ak,Q(3670016)),0);
let safe8=mul(Q(1971200),Q(4096));
ast(add(safe8,Q(304800)),N8);
eq(div(add(N8,mul(Q(4095),Q(304800))),M8),9828000);
eq(div(sub(N8,Q(304800)),M8),9523200);
const t8=BigInt(2*(52416000-304800)-9523200);
req(t8*t8>9523200n**2n+1971200n**2n*2047n,'rank23 periodic norm8');
ast(add(add(div(N10,Q(2n**20n)),Q(3311616)),mul(Q(1482),gamma10)),Q(13000000));
ast(div(Q(2n**26n),M10),Q(100));
ast(qsum([div(add(N10,Q(3311616)),M10),Q(1241088),Q(20705280),Q(15219000)]),Q(47000000));
const C12=Q(3670016),beta=Q(4558848),dd=Q(33226752),ee=Q(1527552),pre=sub(N12,add(beta,dd)),u=div(C12,M12);
eq(mul(Q(1755),gamma12),1527552);
req(ss(u)==='7/32','reciprocal factor');
const optOrigin=div(add(div(pre,M12),add(beta,dd)),add(Q(1),u));
eq(optOrigin,52416000);
const sharpNonlin=div(qsum([div(pre,M12),beta,div(dd,Q(2))]),add(Q(1),u));
eq(sharpNonlin,38784512);
eq(qsum([div(pre,M12),div(dd,Q(2)),mul(u,sharpNonlin)]),51193888);
const double=mul(Q(2),M12);
ast(mul(ee,M12),N12);
const capBig0=qsum([div(N12,double),Q(37785600),div(ee,Q(2))]);
ast(capBig0,TARGET);
const capBigA=div(qsum([div(N12,double),mul(div(dd,Q(2)),sub(Q(1),u)),ee,mul(div(C12,mul(M12,M12)),Q(1962934272n*21521n))]),sub(Q(1),mul(u,u)));
ast(capBigA,Q(30000000));
function Gegen(n,d){
  const a=()=>Array(d+1).fill(Q(0)), One=a(), T=a();
  One[0]=Q(1);T[1]=Q(1);const G=[One,T];
  for(let k=2;k<=d;k++)G.push(a().map((_,i)=>div(sub(mul(Q(2*k+n-4),i?G[k-1][i-1]:Q(0)),mul(Q(k-1),G[k-2][i])),Q(k+n-3))));
  return G;
}
function ev(P,x){return P.reduceRight((a,v)=>add(mul(x,a),v),Q(0));}
let G=Gegen(48,6),F=Array(7).fill(Q(0)),co=[];
F[2]=Q(1,1024);F[4]=Q(-5,64);F[6]=Q(1);
for(let i=6;i>=0;i--){co[i]=div(F[i],G[i][i]);F=F.map((v,j)=>sub(v,mul(co[i],G[i][j])));}
req(co.every(x=>x[0]>=0n),'norm8-in-min12 projective LP');
req(ss(div(qsum(co),co[0]))==='2948400/137','projective LP bound');
let d=2256n,limit=1n+d*(d-1n)/486n;
let G46=Gegen(46,3);
for(let t of [Q(-5,17),Q(1,34)]){
 let val=qsum([Q(1),mul(Q(1173,58),ev(G46[1],t)),mul(Q(4335,58),ev(G46[2],t))]);eq(val,0);
}
for(let t of [Q(-1,2),Q(-1,5),Q(1,10)]){
 let val=qsum([Q(1),mul(Q(65,4),ev(G46[1],t)),mul(Q(1875,4),ev(G46[3],t))]);req(val[0]<=0,'trade LP');
}
req(limit===10468n,'trade cutoff');
console.log('Exact full coefficient table, rooted quotient norm 8,10,12, and bounded trade certificates: PASS');
console.log('minimum binary indices:',ss(M8),ss(M10),ss(M12),'(minimum next power of two for second)');
console.log('norm-12 affine upper bounds:',ss(capBig0),ss(capBigA),'for larger index; index 2^24 at most 52,416,000');


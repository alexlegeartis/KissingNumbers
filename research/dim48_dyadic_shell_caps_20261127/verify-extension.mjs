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
const Z=Q(0),O=Q(1),L=10;
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
function eqq(x,y,msg){req(ss(x)===String(y),msg||ss(x)+' != '+y);}
function qa(x){return x[0]<0n?neg(x):x;}
function power(q,i){return Q(q[0]**BigInt(i),q[1]**BigInt(i));}
function qsum(A){return A.reduce(add,Q(0));}
function ast(x,y,msg){req(cmp(x,y)<0n,msg||ss(x)+' not below '+ss(y));}
function dotq(u,v){return qsum(u.map((x,i)=>mul(x,v[i])));}
const BND=Q(52416000);
const labels=[6,8,10,12,14,16];
const data=[fix(true,1,[[1,2],[2,0]]),fix(false,2,[[2,2],[3,0]]),
 fix(true,2,[[2,2]]),fix(false,3,[[3,2]]),fix(true,3),fix(false,4)];
req(data.every(r=>!r.no&&!r.under),'cusp ranks');
let one12=fix(false,3,[[3,1]]);
const T=new Map(labels.map((s,i)=>[s,data[i].coInf.slice(3,9).map(x=>Q(x))]));
let slope12=one12.coInf.slice(3,9).map((x,i)=>sub(Q(x),T.get(16)[i]));
const type12=a=>T.get(16).map((x,i)=>add(x,mul(Q(a),slope12[i])));
eqq(T.get(14)[4],902734848n);
eqq(T.get(16)[5],19203856800n);
let Del=arr(i=>i===1?O:Z);
for(let k=1;k<L;k++)Del=conv(Del,pow(arr(i=>i===0?O:i===k?Q(-1):Z),24));
let Theta=arr(i=>add(sub(pow(E4,6)[i],mul(Q(1440),conv(pow(E4,3),Del)[i])),mul(Q(125280),pow(Del,2)[i])));
const Ndata=Theta.slice(3,9);
eqq(Ndata[4],15173208925056000n);
eqq(Ndata[5],327259384199748000n);
const candidates=[...labels.map(s=>[s,T.get(s)]),['12:96',type12(96)]];
const w14=[Q(1),Q(185,1758),Q(25,3516),Q(1,9376)],g14=Q(70800,293);
const w16=[Q(1),Q(128,601),Q(1750,48681),Q(8,1803),Q(11,97362)],g16=Q(102400,1803);
for(const [W,g,ind] of [[w14,g14,'23890755584/59'],[w16,g16,'68719476736']]){
  for(const [s,t] of candidates)req(add(dotq(W,t),g)[0]>=0n,'negative margin in '+s);
  eqq(div(add(dotq(W,Ndata),g),g),ind);
}
ast(Q(2n**28n),Q(23890755584n,59));
ast(Q(23890755584n,59),Q(2n**29n));
const R14=[Q(-16187129856n),Q(-189923328),Q(35454976),Q(-3956736)],b14=Q(5351424),d14=Q(-897383424);
const R16=[Q(-716848496640n),Q(22686990336n),Q(-1327497216),Q(95158272)],b16=Q(-119858784),d16=Q(19323715584n);
const masses=s=>[6,8,10].map(j=>Q(s===j?2:0)).concat([Q(s===12?2:0)]);
const gp=2n**48n,mid=2n**24n;
for(const [R,b,d,j,maxJ] of [[R14,b14,d14,4,Q(303,500)],[R16,b16,d16,5,Q(179428077,2097152)]]){
  for(const s of labels){
    const x=T.get(s),m=masses(s),eps=Q(s%4===0?1:-1);
    eqq(qsum([b,mul(eps,d),dotq(R,m)]),ss(x[j]));
    req(cmp(div(dotq(R,x),Q(gp)),maxJ)<=0n,'J bound '+s);
  }
  eqq(slope12[j],ss(R[3]));
  req(dotq(R,slope12)[0]<=0n,'no higher bound on affine branch');
  req(add(b,div(dotq(R,Ndata),Q(gp)))[0]<0n,'zero-class correction');
}
let k14=2n**29n,k16=2n**36n;
const A14=qsum([div(sub(sub(Ndata[4],b14),d14),Q(k14)),div(qa(d14),Q(64)),mul(Q(303,500),Q(gp/k14))]);
ast(A14,Q(43000000));
const A16=qsum([div(sub(sub(Ndata[5],b16),d16),Q(k16)),div(qa(d16),Q(4096)),mul(Q(179428077,2097152),Q(gp/k16))]);
eqq(A16,9830400);
const origin16=add(A16,sub(add(b16,div(dotq(R16,Ndata),Q(gp))),Q(179428077,2097152)));
eqq(origin16,9828000);
ast(A16,BND);
eqq(Q(data[0].co0[4]),36848);eqq(Q(data[2].co0[4]),34300);eqq(Q(data[4].co0[4]),34928);
eqq(Q(data[1].co0[5]),385024);eqq(Q(data[5].co0[5]),393216);eqq(sub(Q(one12.co0[5]),Q(data[5].co0[5])),228);
eqq(add(Q(393216),Q(228*96)),415104);
eqq(Q(data[0].co0[5]),3357774);eqq(Q(data[2].co0[5]),3537800);eqq(Q(data[4].co0[5]),3527064);
eqq(Q(data[1].co0[6]),26300824);eqq(Q(data[5].co0[6]),26099200);eqq(sub(Q(one12.co0[6]),Q(data[5].co0[6])),-896);
console.log('Exact dyadic all-shell extension: PASS.');
console.log('min14, min16 spectral inverse densities 23890755584/59, 68719476736; norm14 <43m, norm16 <=9830400, origin <=9828000.');
console.log('Single residue counts for norms18,20,22,24: 36848, 415104, 3537800, 52416000.');



const g=(a,b)=>b?g(b,a%b):a<0n?-a:a;
const F=(a,b=1n)=>{a=BigInt(a);b=BigInt(b);if(b<0n){a=-a;b=-b}const d=g(a,b);return [a/d,b/d];};
const add=(a,b)=>F(a[0]*b[1]+a[1]*b[0],a[1]*b[1]);
const neg=a=>[-a[0],a[1]];
const sub=(a,b)=>add(a,neg(b));
const mul=(a,b)=>F(a[0]*b[0],a[1]*b[1]);
const div=(a,b)=>F(a[0]*b[1],a[1]*b[0]);
const cp=(a,b)=>a[0]*b[1]-a[1]*b[0];
const fl=a=>a[0]<0n?(a[0]-a[1]+1n)/a[1]:a[0]/a[1];
const ce=a=>-fl(neg(a));
const p=(a,n)=>BigInt(a)**BigInt(n);
const ord=[-3,-2,-1,0,1,2,3];
const pf=(d)=>{let c=[1n], den=1n;for(const x of ord)if(x!==d){const a=Array(c.length+1).fill(0n);c.forEach((v,k)=>{a[k]-=BigInt(x)*v;a[k+1]+=v});c=a;den*=BigInt(d-x);}return c.map(v=>F(v,den))};
function table(b){
 const Q=new Map;
 function wick(i,j){
  if(i===0&&j===0)return 1n;if(i<0||j<0)return 0n;
  const k=i*16+j;if(Q.has(k))return Q.get(k);
  const x=i?BigInt(i-1)*6n*wick(i-2,j)+BigInt(j)*BigInt(b)*wick(i-1,j-1): BigInt(j-1)*6n*wick(0,j-2);
  Q.set(k,x); return x;
 }
 function M(i,j){const z=i+j;
   if(z&1)return F(0);
   let den=1n; for(let k=0;k<z/2;k++)den*=BigInt(48+2*k);
   const P=[[6,b],[-6,-b],[b,6],[-b,-6]];
   return sub(F(52416000n*p(6,z/2)*wick(i,j),den),F(P.reduce((v,[a,c])=>v+p(a,i)*p(c,j),0n)));
 }
 let R=[];
 for(const u of ord){let row=[];for(const v of ord){let a=F(0);const U=pf(u),V=pf(v);
    for(let i=0;i<=6;i++)for(let j=0;j<=6;j++)if(i+j<=10&&M(i,j)[0])a=add(a,mul(mul(U[i],V[j]),M(i,j)));
    row.push([a,mul(U[6],V[6])]);}R.push(row);}
 const op=R[6][6]; R=R.map(row=>row.map(([a,c])=>[sub(a,mul(div(c,op[1]),op[0])),div(c,op[1])]));
 let lo=0n,hi=52416000n;for(const row of R)for(const [a,c] of row){const val=div(neg(a),c);if(c[0]>0n)lo=lo>ce(val)?lo:ce(val);else hi=hi<fl(val)?hi:fl(val);}
 const get=(u,v,k)=>add(R[u+3][v+3][0],mul(R[u+3][v+3][1],F(k)));
 return {lo,hi,get,R};
}

function profile(m,k){
const PM=[[6,k],[k,m]];
const D=(a,b)=>a*a*6+b*b*m+2*a*b*k;
const combos=[];for(let a=-4;a<=4;a++)for(let b=-4;b<=4;b++)if((a||b)&&D(a,b)<=2*m+12)combos.push([a,b,D(a,b)]);
const I=[];
for(let u=0;u<=3;u++)for(let v=-m/2;v<=m/2;v++)if((u||v>=0)&&combos.every(([a,b,d])=>2*Math.abs(a*u+b*v)<=d))I.push([u,v]);
const Q=new Map;
function W(i,j){if(i===0&&j===0)return 1n;if(i<0||j<0)return 0n;const n=i*32+j;if(Q.has(n))return Q.get(n);const a=i?BigInt((i-1)*6)*W(i-2,j)+BigInt(j*k)*W(i-1,j-1):BigInt((j-1)*m)*W(0,j-2);Q.set(n,a);return a;}
const nth=I.length, rows=new Map;
for(let t=0;t<=10;t+=2)for(let a=0;a<=t;a++){
const b=t-a;let den=1n;for(let j=0;j<t/2;j++)den*=BigInt(48+2*j);
let w=I.map(([u,v])=>F(BigInt(u)**BigInt(a)*BigInt(v)**BigInt(b)));
w.push(sub(F(26208000n*6n**BigInt(t/2)*W(a,b),den),F(6n**BigInt(a)*BigInt(k)**BigInt(b))));
for(let i=0;i<nth;i++)if(w[i][0]){if(!rows.has(i)){const p=w[i];rows.set(i,w.map((v,j)=>j<i?F(0):div(v,p)));break;}const q=rows.get(i),r=w[i];for(let j=i;j<=nth;j++)if(q[j][0])w[j]=sub(w[j],mul(r,q[j]));}
if(w.slice(0,nth).every(v=>!v[0])&&w[nth][0])throw [m,k,a,b,"inconsistent"];
}
const N=I.map((_,i)=>i).filter(i=>!rows.has(i));const R=I.map((_,i)=>[...N.map(j=>F(i===j?1:0)),F(0)]);
for(let q=nth-1;q>=0;q--)if(rows.has(q)){const w=rows.get(q);R[q][N.length]=w[nth];for(let j=q+1;j<nth;j++)if(w[j][0])for(let i=0;i<=N.length;i++)if(R[j][i][0])R[q][i]=sub(R[q][i],mul(w[j],R[j][i]));}
const IDs=I.map(([u,v],i)=>u===0&&v>0&&4*v*v>6*m-k*k?i:-1).filter(i=>i>=0);
const T=Array(N.length+1).fill(F(0));for(const i of IDs)for(let j=0;j<T.length;j++)T[j]=add(T[j],R[i][j]);
return {I,N,R,IDs,T};
}

function cyclic(){
 const {I,N,R}=profile(8,1);const C=R.map(v=>v.map(t=>Number(t[0])/Number(t[1])));
 const res=[];for(let a=0;a<=3;a++)for(let b=0;b<=5;b++)for(let c=0;c<=1;c++){const P=C.map(q=>47*(q[0]*a+q[1]*b+q[2]*c)+q[3]);if(P.every(t=>t>=0&&t%47===0))res.push([a,b,c,P]);}
 const folded=I.flatMap(([s,t],i)=>s||t?[[s,t,i],[-s,-t,i]]:[[s,t,i]]);
 const min=res[0][3].slice();const max=min.slice();res.forEach(p=>p[3].forEach((x,j)=>{min[j]=Math.min(x,min[j]);max[j]=Math.max(x,max[j]);}));
 const n=(s,t)=>8*s*s-2*s*t+6*t*t,j=(s,t)=>((8*s-t)%47+47)%47;
 const P=new Map; for(let a=-12;a<=12;a++)for(let b=-12;b<=12;b++){const h=j(a,b);if(!P.has(h)||n(a,b)<P.get(h).n)P.set(h,{n:n(a,b),cells:[a+","+b]});else if(n(a,b)===P.get(h).n)P.get(h).cells.push(a+","+b);}
 const pool=folded.map(([s,t,i])=>({s,t,i,j:j(s,t),q:n(s,t),lo:min[i]*(i?1:2),hi:max[i]*(i?1:2)}));
 return {res,pool,P,folded,I};
}

function shell8(R){
 const prof=R.res.find(a=>a[2]===1)[3];
 const shell=R.pool.map(x=>[x.s,x.t,prof[x.i]*(x.i?1:2)]).filter(r=>r[2]);shell.push([6,1,1],[-6,-1,1]);
 const N6=52416000n,N8=39007332000n;
 const M6=new Map,M8=new Map, pols=new Map;
 const k=(a,b)=>32*a+b;
 function lap(a,b){
  if(pols.has(k(a,b)))return pols.get(k(a,b));
  let A=[new Map([[k(a,b),1n]])];let mat=[[2,0,6],[1,1,2],[0,2,8]];
  for(let j=1;j<=((a+b)>>1);j++){
   const Q=new Map;for(const [t,c] of A.at(-1)){const i=t>>5,l=t%32;for(const [d,e,m] of mat){if(i<d||l<e)continue;let z=BigInt(m)*c;
      for(let s=0;s<d;s++)z*=BigInt(i-s);for(let s=0;s<e;s++)z*=BigInt(l-s);
      const h=k(i-d,l-e);Q.set(h,(Q.get(h)??0n)+z);}}
   A.push(Q);
  }pols.set(k(a,b),A);return A;
 }
 function m6(a,b){if(!M6.has(k(a,b)))M6.set(k(a,b),F(shell.reduce((t,[i,j,n])=>t+BigInt(n)*p(i,a)*p(j,b),0n)));return M6.get(k(a,b));}
 M8.set(0,F(N8));
 function horn(T, h){let Z=F(0);for(const [s,c] of T)Z=add(Z,mul(F(c),h(s>>5,s%32)));return Z;}
 function m8(a,b){const d=a+b;if(d&1)return F(0);if(M8.has(k(a,b)))return M8.get(k(a,b));
   const q={2:0,4:0,6:0,8:0,10:0,12:-72,14:0,16:168,18:-576,20:408,22:-336}[d];
   if(q===undefined)throw "unknown harmonic";
   let aa=F(1), hh=m6(a,b), ss=F(0);
   const L=lap(a,b);
   for(let j=1;j<L.length;j++){aa=div(aa,F(-4*j*(d-j+23)));if(q)hh=add(hh,mul(mul(aa,F(6n**BigInt(j))),horn(L[j],m6)));ss=sub(ss,mul(mul(aa,F(8n**BigInt(j))),horn(L[j],m8)));}
   const out=add(mul(F(q),hh),ss); if(out[1]!==1n)throw "bad integer moment";M8.set(k(a,b),out);return out;
 }
 const grid=(lo,hi)=>Array.from({length:hi-lo+1},(_,i)=>i+lo);
 function lf(roots,t){let C=[1n],den=1n;for(const x of roots)if(x!==t){
  let a=Array(C.length+1).fill(0n);C.forEach((v,k)=>{a[k]-=BigInt(x)*v;a[k+1]+=v});C=a;den*=BigInt(t-x);
 }return C.map(c=>F(c,den));}
 const X=grid(-4,4),Y=grid(-5,5);
 const Ms=Array.from({length:X.length},(_,a)=>Array.from({length:Y.length},(_,b)=>sub(m8(a,b),F(p(1,a)*p(8,b)+p(-1,a)*p(-8,b)))));
 const out=[];for(const s of X)for(const t of Y){const A=lf(X,s),B=lf(Y,t);let v=F(0);
 for(let a=0;a<A.length;a++)for(let b=0;b<B.length;b++)if(A[a][0]&&B[b][0])v=add(v,mul(mul(A[a],B[b]),Ms[a][b]));
 if(v[1]!==1n||v[0]<0n||v[0]%47n)throw [s,t,"bad norm8 census",v.map(String)];
 out.push([s,t,Number(v[0])]);
 }
 const m=out.reduce((v,x)=>v+BigInt(x[2]),2n);if(m!==N8)throw "wrong norm8 count";
 out.push([1,8,1],[-1,-8,1]);
 return out;
}

function mixed47(R,S8){
 const X=R.res.find(x=>x[2]===1)[3];
 const sn=(s,t)=>8*s*s-2*s*t+6*t*t,j=(s,t)=>((8*s-t)%47+47)%47;
 const M6=R.pool.map(x=>[x.s,x.t,X[x.i]*(x.i?1:2)]).filter(r=>r[2]);M6.push([6,1,1],[-6,-1,1]);
 const M=[M6,S8],Sh=M.map(es=>new Set(es.filter(x=>x[2]>0).map(x=>x[0]+","+x[1])));
 const norm=new Map;
 for(const [k,x] of R.P){norm.set(k,{n:x.n,lo:x.cells.some(h=>!Sh[0].has(h)&&!Sh[1].has(h))?470:x.cells.some(h=>!Sh[0].has(h))?376:282});}
 const uniq=new Map;
 M.forEach((array,l)=>array.forEach(([s,t,w])=>{
   const F=sn(s,t),r=(l?376:282)-F,q=282-r; if(!w||r<=0||q<0)return;
   const key=j(s,t)+":"+q,x={s,t,r,q,j:j(s,t),w,l};
   if(uniq.has(key)&&uniq.get(key).w!==w)throw "alias";
   (uniq.get(key)??(uniq.set(key,x),x)).alts?.push([l,s,t]);
 }));
 const A=[...uniq.values()].flatMap(f=>f.q?[{...f,sg:1},{...f,sg:-1}]:[{...f,sg:0}]);
 A.push({s:1,t:8,r:0,q:282,j:0,w:1,l:2,sg:1},{s:-1,t:-8,r:0,q:282,j:0,w:1,l:2,sg:-1});
 const n=A.length,C=Array.from({length:n},()=>new Uint8Array(n));
 for(let a=0;a<n;a++)for(let b=a+1;b<n;b++){
  const x=A[a],y=A[b],z=norm.get((x.j-y.j+47)%47);
  const U=x.r===0||y.r===0?282:x.q+y.q+(x.j===y.j&&x.q===y.q?-282:z.lo-282-z.n);
  const Z=4*x.q*y.q,M=x.sg*y.sg;
  const ok=M===1?U>=0&&U*U>=Z:M===-1?U>=0||U*U<=Z:U>=0;
  if(ok)C[a][b]=C[b][a]=1;
 }
 return {A,C,norm};
}

function exactClique(A,C,initial=[],steps=200000){
 const n=A.length;
 const ids=Array.from({length:n},(_,i)=>i).sort((i,j)=>A[j].w-A[i].w);
 const R=C.map((a,i)=>ids.reduce((x,j,k)=>a[j]?x|(1n<<BigInt(k)):x,0n)); const P=ids.map(i=>R[i]);
 const W=ids.map(i=>A[i].w);
 const bits=ids.map((_,i)=>1n<<BigInt(i));
 const which=new Map(bits.map((b,i)=>[b,i]));
 const pos=v=>which.get(v&(-v));
 const totalMask=(1n<<BigInt(n))-1n;
 const Q=P.map((p,i)=>totalMask^(p|bits[i]));
 let B=initial.reduce((s,i)=>s+A[i].w,0), best=initial,go=0,complete=true;
 const memo=new Map;
 function upper(m,z){
   let k=m,sum=0; while(k){const i=pos(k);sum+=W[i];k^=bits[i];}
   if(sum+z<=B)return sum;
   const w=W.slice();
   let ret=0,assigned=sum;
   while(m){
    let v=m, group=[],q=m,lb=Infinity,ls=0;
    while(v){const i=pos(v); if(w[i]<lb)lb=w[i];ls+=w[i];group.push(i);v&=Q[i];}
    ret+=lb;
    for(const i of group){w[i]-=lb;assigned-=lb;if(!w[i])m^=bits[i];}
    if(ret+assigned+z<=B)return ret+assigned;
   }
   return ret;
 }
 function rec(a,z,u){
   if(++go>steps){complete=false;return;}
   for(let m=a;m; m&=m-1n){const i=pos(m);if(!(Q[i]&a)){z+=W[i];u=[...u,i];a^=bits[i];}}
   if(!a){if(z>B){B=z;best=u.map(i=>ids[i]);}return;}
   if((memo.get(a)??-1)>=z)return;memo.set(a,z);
   if(z+upper(a,z)<=B)return;
   let v=pos(a);
   for(let m=a&(a-1n); m; m&=m-1n){const i=pos(m);if(W[i]<W[v]*.35)break;if(W[i]>W[v]&&((P[v]&a)&~P[i])===0n&&!(P[i]&bits[v])){v=i;}}
   rec(a&P[v],z+W[v],[...u,v]);if(complete)rec(a^bits[v],z,u);
 }
 rec(totalMask,0,[]);
 return {B,best,go,complete};
}

function intervals47(R,S,L=8,center=true){
 const {A,C,norm}=mixed47(R,S);
 const cent=A.findIndex(v=>!v.q);
 const B=(x,y)=>x.r+y.r-(x.j===y.j&&x.q===y.q?0:(()=>{const z=norm.get((x.j-y.j+47)%47);return z.lo-z.n})());
 const le=(x,i,y,j)=>{
  if(!x.r)return x.sg!==y.sg||4*j*j*y.q<=y.r*L*L+j*j*y.q;
  if(!y.r)return le(y,j,x,i);
  const b=BigInt(B(x,y)*L*L),a=BigInt(2*i*j*x.sg*y.sg),Y=BigInt(x.q*y.q),J=BigInt(x.r*L*L+i*i*x.q)*BigInt(y.r*L*L+j*j*y.q);
  const X=a*a*Y;
  if(b<=0n&&(a<=0n||b*b>=X)||b>=0n&&a<0n&&b*b<=X)return true;
  const D=J-b*b-X,V=2n*b*a;
  if(V===0n)return D>=0n;
  if(V>0n)return D>=0n&&D*D>=V*V*Y;
  return D>=0n||D*D<=V*V*Y;
 };
 const nodes=[];
 for(let i=0;i<A.length;i++){const x=A[i];if(i===cent&&center===false)continue;
   const I=(!x.r||!x.q)?[[L,L]]:Array.from({length:L},(_,i)=>[i,i+1]);
   for(let [lo,hi] of I){
     if(center&&i!==cent&&!le(x,hi,A[cent],L))continue;
     if(center&&i!==cent&&x.r&&x.q){
       const v=A[cent],U=B(x,v);
       if(U*U===282*282){lo=L;hi=L;}
     }
     if(nodes.length&&nodes.at(-1).i===i&&nodes.at(-1).lo===lo&&nodes.at(-1).hi===hi)continue;
     nodes.push({...x,i,lo,hi});
   }
 }
 const n=nodes.length;const D=Array.from({length:n},()=>new Uint8Array(n));
 for(let a=0;a<n;a++)for(let b=a+1;b<n;b++){
  const x=nodes[a],y=nodes[b];if(x.i===y.i)continue;
  if(le(x,x.lo,y,y.lo)||le(x,x.hi,y,y.lo)||le(x,x.lo,y,y.hi)||le(x,x.hi,y,y.hi))D[a][b]=D[b][a]=1;
 }
 return {nodes,D,le,cent};
}

function demand(q,msg){if(!q)throw Error(msg);}
const CY=cyclic(),P=CY.res;
demand(P.length===2&&P.map(x=>x.slice(0,3).join(",")).join("|")==="0,0,1|3,5,0","rank-2 minimal shell profiles");
const actual=P[0][3];
for(const [s,t,coef] of [[1,-4,[-4,0,-15,705]],[3,-2,[0,-6,-35,1645]]]){
 const {I,R}=profile(8,1),i=I.findIndex(x=>x[0]===s&&x[1]===t);
 demand(R[i].every(([a,b],k)=>b===1n&&a===BigInt(coef[k])),"finite parameter bounds");
}
demand(Math.max(...Array.from(CY.P.values()).map(h=>h.n))===144,"shortest normal cutoff");
demand(CY.P.size===47,"all discriminant classes");
const S8=shell8(CY);
demand(S8.some(([s,t,n])=>s===0&&t===4&&n===7945350),"large 8-shell");
const init=(A)=>A.map((r,i)=>r.l===0&&r.s===0&&r.sg===Math.sign(r.t)&&(!r.hi||r.hi===4)?i:-1).filter(i=>i>=0);
const eq=mixed47(CY,S8),eq0=init(eq.A);
demand(eq.A.length===107&&eq0.every(i=>eq0.every(j=>i===j||eq.C[i][j])),"mixed fixed norm witness");
const q=exactClique(eq.A,eq.C,eq0);
demand(q.complete&&q.B===23766960,"complete fixed-norm upper");
const e=[];
for(const co of [true,false]){
 const {nodes,D,le}=intervals47(CY,S8,4,co),Z=co?init(nodes):[];
 demand(Z.every(i=>Z.every(j=>i===j||le(nodes[i],4,nodes[j],4))),"literal center-section witness");
 const C=exactClique(nodes,D,Z);
 demand(C.complete&&C.B===(co?23766960:20919607),"continuous interval upper "+co);
 e.push({center:co,nodes:nodes.length,cap:C.B,branch:C.go});
}
const rows=[];for(let b=0;b<=3;b++){
 const {lo,hi,get}=table(b);
 const rng=(u)=>[get(u,3,lo),get(u,3,hi)].map(x=>{demand(x[1]===1n,"exact edge count");return Number(x[0]);}).sort((x,y)=>x-y);
 rows.push({b,lo:Number(lo),hi:Number(hi),range:ord.map(rng)});
}
demand(JSON.stringify(rows.map(r=>r.range[3]))===JSON.stringify([[17216,19136],[15520,15680],[7936,7936],[1656,1656]]),"minimum blocker degrees");
demand([1,2,3].map(j=>rows[0].range[j+3][1]).join(",")==="9720,576,96","maximum one-sided degrees");
const L=add(add(F(9720,15520),F(576,7936)),F(96,1656));demand(cp(L,F(1))<0n,"one side charge");
console.log({status:"PASS",members:shell8(CY).filter(x=>x[2]).length,signedOptions:eq.A.length,fixedNorm:q.B,continuous:e,oneSideLoad:L.map(String)});


const Z={"u":"000120000220020000000120000000000011102021111220","v":"122210000100020000000100000000000010201021001100","z":"000000000000000000000000000000000000000000000000","p":"000020000220020000000020000000000000002020000220","checks":[0,1],"s":[3,4,0,0]};
const R=["222002220210100201201001","022200222021010020120101","002220022202101002012011","111220112010110202000200","011122011201011020200020","220111121000201000221001","022011112100020100022101","002201111210002010002211","111221221211200000102220","200121012001220201211221","212011021110222221022122","210200022021122120000210","102021112022012011102022","202201001112001102011200","020220100111200110201120","221021200221220212221111","100100200112022120021110","121011100101102011101112","120102220100010000212111","120011002100201102120210","120002210000220212011022","201002111210122222102100","020100211121012222210210","110011101202001021020022"].map(z=>[...z].map(Number));
const G=R.map((r,i)=>Array.from({length:24},(_,j)=>+(i===j)).concat(r));
const m=x=>(x%3+3)%3;
const dot=(a,b)=>m(a.reduce((s,x,i)=>s+x*b[i],0));
const word=s=>[...s].map(Number);
const wt=x=>x.reduce((s,k)=>s+(k!==0),0);
function chk(q,label){if(!q)throw Error(label);}

function pg(pairs){const a=Array(24).fill(0);for(const [i,x] of pairs)a[i]=m(x);return a;}
const gg=pg([[23,1],[20,1],[18,-1],[17,1],[15,-1],[12,1],[10,1],[9,-1],[7,-1],[6,-1],[5,-1],[2,-1],[1,-1],[0,-1]]);
const hh=pg([[23,1],[22,1],[21,1],[18,1],[17,1],[16,1],[14,1],[13,-1],[11,-1],[8,1],[6,-1],[5,1],[3,-1],[0,-1]]);
function mul(a,b){const c=Array(a.length+b.length-1).fill(0);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]=m(c[i+j]+a[i]*b[j]);return c;}
const f=mul([2,1],gg),prod=mul(f,hh);
chk(prod.length===48&&prod.every((x,i)=>x===(i===0?2:i===47?1:0)),"x^47-1 factorization");
let ord=1,y=3;while(y!==1){ord++;y=y*3%47;}chk(ord===23,"nontrivial cyclotomic cosets have degree 23");
chk(m(gg.reduce((a,b)=>a+b))!==0&&m(hh.reduce((a,b)=>a+b))!==0,"quadratic-residue factors exclude x-1");
const canon=Array.from({length:23},(_,i)=>Array.from({length:48},(_,j)=>j<47?(f[j-i]??0):0)).concat([Array(48).fill(1)]);
for(let j=0;j<24;j++){
 const k=canon.findIndex((r,i)=>i>=j&&r[j]);chk(k>=0,"systematic QR pivot");
 [canon[k],canon[j]]=[canon[j],canon[k]];if(canon[j][j]===2)canon[j]=canon[j].map(x=>m(2*x));
 for(let i=0;i<24;i++)if(i!==j&&canon[i][j]){
  const h=canon[i][j];canon[i]=canon[i].map((x,c)=>m(x-h*canon[j][c]));
 }
}
chk(JSON.stringify(canon)===JSON.stringify(G),"canonical extended QR systematic conversion, no coordinate relabel");

const u=word(Z.u),v=word(Z.v),phi=word(Z.z);
function inQ(x){return x.length===48&&x.every((t,j)=>j<24||t===m(x.slice(0,24).reduce((s,a,i)=>s+a*R[i][j-24],0)));}
chk(G.every(x=>G.every(y=>dot(x,y)===0)),"QR self-duality [I|R]");
chk(inQ(u)&&inQ(v),"new projective pencil members");
chk(wt(u)===18&&wt(v)===15&&v.slice(0,3).every(Boolean)&&u.slice(0,3).every(x=>x===0),"marked weights");
for(let alpha of [1,2])chk(wt(v.map((x,i)=>m(x+alpha*u[i])))===15,"other pencil weights");
const S=u.map((x,i)=>x?i:-1).filter(i=>i>=0);
const D=G.slice(3).map(z=>S.map(i=>z[i]));
function rank(A){
 A=A.map(r=>r.slice());let r=0;
 for(let j=0;j<A[0].length&&r<A.length;j++){
  const i=A.findIndex((a,k)=>k>=r&&a[j]);if(i<0)continue;
  [A[r],A[i]]=[A[i],A[r]];
  let b=A[r][j];if(b===2)A[r]=A[r].map(x=>m(x*b));
  for(let k=r+1;k<A.length;k++)if(A[k][j]){let a=A[k][j];A[k]=A[k].map((x,i)=>m(x-a*A[r][i]));}
  r++;
 }return r;
}
chk(rank(D)===16,"dense QR projection rank");
const p=phi.map((x,i)=>m(x-(u[i]===2)));
const q=phi.map((x,i)=>m(x-(u[i]===1)));
for(const z of [p,q])chk(dot(z,u)===0&&dot(z,v)===1,"defective syndrome is (0,1)");
chk(S.every(i=>u[i]!==0),"no check column reaches (0,1)");
chk(4*24**2<45*54,"dense strict cross angle");
const tri=[[2,0],[-1,1],[-1,-1]],roots=tri.map((x,i)=>x.map((z,k)=>z-tri[(i+2)%3][k]));
const gram=(a,b)=>a[0]*b[0]+3*a[1]*b[1]; // four times raw
for(let i=0;i<3;i++){
 chk(gram(roots[i],roots[i])===12,"A2 root length");
 chk(tri.map(x=>gram(x,roots[i])).join()===Array.from({length:3},(_,j)=>j===i?6:j===(i+2)%3?-6:0).join(),"root/trit angular table");
 for(let j=0;j<3;j++){
  chk(gram(roots[i],roots[j])===(i===j?12:-6),"same outer roots tangent");
  chk(24+2*gram(roots[i],roots[j])===(i===j?48:12),"opposite outer 4*d^2");
 }
}
const K=7906;
for(let b=3;b<=6;b++){
 const same=b+3,opp=12-2*b;
 chk(same+opp===15-b,"overlap controls");
 chk(K*2**opp<3**(same+opp),"translation expectation <1");
}
chk(4*(5*1.5)**2<54*5&&4*(4**2*3)<54*4&&4*(9*3)<54*3&&4*12<54*2,"old small-shell walls");
const oldDense=2*3**21,oldSparse=130944*K;
const prior=21956420122,oldSmall=prior-oldDense-oldSparse;
chk(oldSmall===443520+16020-3960+864*16+131*8,"old small-shell stock arithmetic");
const total=prior+1;
console.log(JSON.stringify({u:Z.u,v:Z.v,phase:Z.z,qrProjectionRank:16,dense:oldDense,sparse:oldSparse,small:oldSmall,added:1,prior,total}));


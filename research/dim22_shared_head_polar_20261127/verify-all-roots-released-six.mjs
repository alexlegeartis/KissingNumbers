
const N=8,unit=i=>Array.from({length:N},(_,j)=>i===j?1:0),vsum=(...rows)=>rows.reduce((a,b)=>a.map((x,i)=>x+b[i]),Array(N).fill(0)),times=(a,n)=>a.map(x=>x*n);
const P=Array.from({length:N-2},(_,i)=>unit(i)),[a,b,c]=P,z=unit(N-2),Y=unit(N-1);
const row=(a,b)=>[...a,b];
function* ix(n,k,j=0,a=[]){if(k===0){yield a;return;}for(let i=j;i<=n-k;i++)yield*ix(n,k-1,i+1,[...a,i]);}

function mat(zslot,exc,par){
 const S=vsum(...P,z),order=[...P.slice(0,zslot),z,...P.slice(zslot)];
 let A=Array.from({length:N-2},(_,i)=>row(vsum(order[i+1],times(order[i],-1)),0));
 A.push(row(times(order.at(-1),-1),0),row(times(Y,-1),0),
        row(vsum(z,a,b,c),8),row(vsum(a,b,c,P[3]),8),row(vsum(a,Y),8));
 if(exc===null){A.push(row(vsum(a,P.at(-1),Y),8),row(vsum(S,times(order[0],2),Y),16));return A;}
 let E=order.filter(x=>exc.includes(x.indexOf(1))),H=order.filter(x=>!E.includes(x)), HH=vsum(...H),M=H[0],m=H.at(-1),M0=H.find(x=>x!==z),E0=E.find(x=>x!==z);
 if(H.length===5&&H.includes(z)) A.push(row(vsum(M0,E0,Y),8));
 if(H.length!==5)A.push(row(vsum(a,P.at(-1),Y),8));
 const even=H.length%2===0;
 if(even){
  A.push(row(vsum(S,times(M,2),times(m,par?0:-2),Y),16));
  if(E.length)A.push(row(vsum(S,times(E[0],2),times(m,par?-2:0),Y),16));
 }else if(par){
   A.push(row(vsum(S,times(M,2),Y),16));
   if(E.length)A.push(row(vsum(S,times(E[0],2),times(Y,-1)),16),row(vsum(S,times(E[0],2),times(m,-2),Y),16));
 }else{
   A.push(row(vsum(S,times(M,2),times(Y,-1)),16),row(vsum(S,times(M,2),times(m,-2),Y),16));
   if(E.length)A.push(row(vsum(S,times(E[0],2),Y),16));
 }
 if(H.length===5&&E.includes(z))if(par)A.push(row(vsum(HH,times(m,-2),Y),8),row(vsum(HH,times(Y,-1)),8));
 else A.push(row(vsum(HH,Y),8));
 return A;
}

function req(ok,why){if(!ok)throw Error(why);}
function integerSolve(A) {
 const B=A.map(row=>row.map(BigInt)); let prev=1n;
 for(let j=0;j<N-1;j++) {
   let i=j;while(i<N&&!B[i][j])i++;if(i===N)return null;
   if(i!==j)[B[i],B[j]]=[B[j],B[i]];
   const d=B[j][j];
   for(let i=j+1;i<N;i++){
     for(let k=j+1;k<=N;k++){const q=B[i][k]*d-B[i][j]*B[j][k];
      req(q%prev===0n,'Bareiss exact'); B[i][k]=q/prev; }
     B[i][j]=0n;
   }
   prev=d;
 }
 let den=B[N-1][N-1];if(!den)return null;
 let C=Array(N).fill(0n); 
 for(let i=N-1;i>=0;i--){let h=den*B[i][N];for(let j=i+1;j<N;j++)h-=B[i][j]*C[j];
   req(h%B[i][i]===0n,'scaled triangular solve exact');C[i]=h/B[i][i];}
 if(den<0n){den=-den;C=C.map(w=>-w);}
 return [C,den];
}


function chamber(z,q,p){
 let A=mat(z,q,p),n=0,cnt=0,eq=0;
 for(let I of ix(A.length,N)){
  n++;
  let u=integerSolve(I.map(i=>A[i]));if(!u)continue;
  let [C,d]=u;
  if(A.some(row=>C.reduce((u,v,i)=>u+v*BigInt(row[i]),0n)>BigInt(row[N])*d))continue;
  cnt++;
  let z2=C.reduce((s,c,i)=>s+(i===N-1?1n:3n)*c*c,0n);
  req(z2<=96n*d*d,'norm above 32: '+[z,q,p,C,d].join(' '));
  if(z2!==96n*d*d)continue;
  eq++;
  const is=w=>C.every((x,i)=>x===BigInt(w[i])*d);
  req(is([4,4,0,0,0,0,0,0]) || is([4,0,0,0,0,0,4,0]) || (is([2,2,2,2,0,0,2,6]) && q!==null && q.join()==='4,5' && p===0),'bad norm-32 ray');
 }
 return [n,cnt,eq];
}
const verify=[];
for(let z=0;z<=6;z++){
 for(let q of [null,...[0,1,2].flatMap(k=>[...ix(7,k)])]){
  for(let p of q===null?[0]:[0,1])verify.push(chamber(z,q,p));
 }
}
req(verify.length===413,'413 chambers');
console.log(JSON.stringify({chambers:verify.length,bases:verify.reduce((s,r)=>s+r[0],0),
 feasibleBases:verify.reduce((s,r)=>s+r[1],0),norm32Bases:verify.reduce((s,r)=>s+r[2],0),
 verdict:'PASS: only original roots or line-cap rays, even after releasing all roots through support 6'},null,2));


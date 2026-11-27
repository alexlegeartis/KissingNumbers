const assert=require('node:assert/strict');
const out=(function(){

const D=16384,M=262144,D2=D*D,M2=M*M;
let visited=0,queries=0,count=0,maxd=0,leaf=[];
function Box(x,y,s){
 const X=x+s,Y=y+s, L=19*D2+x*x+y*y,H=19*D2+X*X+Y*Y;
 const dis=(p,a,b)=>p<a?a-p:p>b?p-b:0;
 const g=Math.max(6*D+dis(D,x,X),6*D+dis(D,y,Y),4*D+x+y,6*D+Math.max(0,x-Y,y-X),2*(x+y),2*(D+Math.max(x,y)));
 if(g*g>2*H)return null;
 const lo=Math.floor(Math.sqrt(L)),hi=Math.floor(Math.sqrt(H))+1;
 if(lo*lo>L||hi*hi<=H)throw Error('sqrt bounds');
 const vf=(a,b)=>{const q=Math.floor(a/b);if(q*b>a||(q+1)*b<=a)throw Error('quotient');return q;};
 const v=(a,b)=>vf(a*M,b),ceil=(a,b)=>vf(a*M+b-1,b);
 return {x,y,s,loF:v(D,hi),loX:v(x,hi),loY:v(y,hi),hiF:ceil(D,lo),hiX:ceil(X,lo),hiY:ceil(Y,lo)};
}
function split(A){
 if(!A.s||A.s%2)throw Error('split');
 if(A.child===undefined){const t=A.s/2;A.child=[Box(A.x,A.y,t),Box(A.x+t,A.y,t),Box(A.x,A.y+t,t),Box(A.x+t,A.y+t,t)].filter(Boolean);}
 return A.child;
}
function query(A,B,t,lim){
 queries++;
 const u=t===0?9*A.loF*B.loF+A.loX*B.loX+A.loY*B.loY:13*A.loF*B.loF-A.hiX*B.hiX-A.hiY*B.hiY;
 if(2*u>M2)return false;
 if(B.s<=lim)return true;
 const C=split(B);
 if(t===0){for(const b of C)if(query(A,b,t,lim))return true;}
 else {for(let i=C.length-1;i>=0;i--)if(query(A,C[i],t,lim))return true;}
 return false;
}
const root=Box(0,0,2*D);
function go(A,d){
 visited++;maxd=Math.max(maxd,d);
 if(!query(A,root,0,Math.max(1,A.s/2)))return;
 if(!query(A,root,1,Math.max(1,A.s/2)))return;
 if(A.s===1){count++;if(leaf.length<8)leaf.push([A.x,A.y]);return;}
 for(const B of split(A))go(B,d+1);
}
go(root,0);return {visited,queries,count,maxd,leaf};
})();
assert.deepEqual(out,{visited:2673,queries:462334,count:0,maxd:13,leaf:[]});
console.log('Exact 21D two-mark continuous non-intersection:',out);
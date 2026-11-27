const assert=require('node:assert/strict');
const b=[1973n,1664n,1541n,1473n,1431n,1403n,1384n,1371n,1364n,1360n],d=1000n;
for(let k=1n;k<=10n;k++){
 const n=b[Number(k-1n)],X=n*n,D=d*d;
 assert(k*(21n-k)*(X-D)**2n<43n*X*D);
 const cut=(j,m)=>{assert(2n*((21n-j)*D+j*X)<(6n*d+m*(n-d))**2n);assert(2n*j<m*m);};
 if(k===1n){assert(9n*X-14n*n*d-D>0n);assert(18n*n>14n*d);continue;}
 if(k===2n)assert(8n*n>11n*d);else cut(k,k<=7n?(k<4n?k:4n):5n);
 assert((18n-k)*X>(21n-k)*D);
 if(k>=8n){
  cut(21n-k,k===10n?5n:6n);
  assert((k-3n)*X>k*D);
 }
}
const moment=k=>28*k*(k-1)-6*80*k+12*210;
assert.deepEqual([moment(8),moment(9)],[248,216]);assert.notEqual(80*10,4*210);
console.log('Exact reciprocal 21D no-go checked.');

for(const a of [24n,75n])for(const b of [24n,75n]){
 const d=40n,D=d*d;
 assert(4n*(10n*D+a*b)**2n>(20n*D+a*a)*(20n*D+b*b));
}
assert(58n*25n>38n*38n); // sqrt58-7>3/5
assert(10n*64n>25n*25n); // 5-sqrt10<15/8
console.log('Exact arbitrary pointwise single-coordinate no-go checked.');

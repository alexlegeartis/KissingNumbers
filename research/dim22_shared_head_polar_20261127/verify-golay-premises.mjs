const pc=x=>{let y=0;for(;x;x&=x-1)y++;return y;};const g=[0,1,5,6,7,9,11].reduce((a,k)=>a|(1<<k),0);let G=[0];for(let k=0;k<12;k++){let p=g<<k;p|=(pc(p)%2)<<23;G.push(...G.map(x=>x^p));}
const O=G.filter(w=>!(w>>>21)&&pc(w)==8), C=G.filter(w=>!(w>>>21)), P=G.filter(w=>(w>>>21)==7&&pc(w)==8).map(w=>w&0x1fffff),F=0x1fffff;
const seq=n=>Array.from({length:n},(_,i)=>i);function* comb(n,k,j=0,b=0){if(!k){yield b;return;}for(let i=j;i<=n-k;i++)yield*comb(n,k-1,i+1,b|(1<<i));}
let S=G.map(x=>x&F),L=S.filter(x=>pc(x)===5),wt={};for(let x of G)wt[pc(x)]=(wt[pc(x)]||0)+1;
const ok=(b,msg)=>{if(!b)throw Error(msg);};
ok(G.length===4096&&new Set(G).size===4096&&new Set(S).size===4096,'Golay counts');
ok(JSON.stringify(wt)==='{"0":1,"8":759,"12":2576,"16":759,"24":1}','Golay weights');
for(let i=0;i<12;i++)for(let j=0;j<12;j++)ok(pc(G[1<<i]&G[1<<j])%2===0,'Golay self dual');
ok(C.length===512&&O.length===210&&L.length===21,'shortened/punctured counts');
for(let i=0;i<L.length;i++)for(let j=i+1;j<L.length;j++)ok(pc(L[i]&L[j])===1,'projective line intersections');
for(let a of comb(21,4))ok(O.filter(x=>(x&a)===a).length>=2,'four-set lower incidence');
ok(4*(21*20/2)+128*O.length+2*21*C.length+2*16*L.length===49896,'incumbent size');
console.log('Finite Golay, line and shortened-octad premises PASS.');


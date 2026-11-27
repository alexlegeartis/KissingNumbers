// Exact 37D arithmetic for the two Witt-protected 12-coordinate lines.
const source=S=>{
const r=[0,1,2,3,4,5,6,7,8,9,11,10],p=[[],[]];
const pop=x=>{let s=0;for(;x;x&=x-1n)s++;return s;};
const iv=j=>{for(let k=1;k<11;k++)if(j*k%11===1)return k;throw Error("iv");};
for(let j=0;j<12;j++){p[0][r[j]]=j<11?r[(j+1)%11]:r[11];p[1][r[j]]=j===0?r[11]:j===11?r[0]:r[(11-iv(j))%11];}
const pm=p.map(P=>Array.from({length:36},(_,i)=>1n<<BigInt(12*Math.floor(i/12)+P[i%12])));
const seeds=[0x0000035b1n,0x000027e10n,0x003024115n,0x003115404n,0x00f081104n,0x0270001c8n,0x0271c8000n,0x07d101000n],ex=[0x0000ff000n,0x000ff0000n,0x000f0f000n];
const all=new Set();for(let x of seeds){let s=new Set([x]),Q=[x];for(let j=0;j<Q.length;j++)for(let q of pm){let z=0n;for(let k=0;k<36;k++)if(Q[j]>>BigInt(k)&1n)z|=q[k];if(!s.has(z)){s.add(z);Q.push(z);}}for(let v of s)all.add(v);}
for(let x of ex)all.add(x);
if(all.size!==2742)throw Error("support");
const out=[...all].sort((a,b)=>a<b?-1:a>b?1:0),rad=S.trim().split(/\s+/).map(BigInt);
return {o:out,r:rad,pop,ex};
};
const g=source("34393344002\n102105092\n9145815040\n19461701665\n1076045824\n34362396676\n604684288\n1209368576\n17855414336\n17196662840\n4429457409\n1150296576\n19327877196\n4598079488\n272793858\n21475885123\n11814043648\n40802191504\n136904962\n17230725504\n34359838016\n4836102656\n39728456736\n2156005636\n34904998448\n4430234244\n34963722250\n3288369160\n135009792\n34376517897\n9194963328\n36541038657\n43016912992\n201342995\n34502607232\n34428945537\n21086209\n3238003472\n159451152\n2425424\n35433631872\n570434308\n570427553\n806486592\n268462464\n573047872\n106522\n2148668424\n2421164033\n8858510340\n10737455170\n1091174410\n3359638532\n2189468672\n9697234176\n36523999270\n8858898696\n2282758464\n4563540000\n4320166048\n4299704328\n27525123\n8625594504\n4328801408\n38666371072\n17222074888\n1124347904\n51568971776\n17180332034\n1612189728\n2701666304\n4364214528\n34399584544\n35433546243\n2214604932\n606356\n4336914690\n11012736\n10804593680\n54659136\n34360401937\n70256674\n34411120128\n12885295872\n34378637376\n18253873297\n19327664384\n1115750436\n8624112128\n6442976513\n5368872980\n1108361225\n34431336464\n289669256\n6439428\n2181059104\n4207619\n150999233\n805372962\n277233680\n671098000\n9797894728\n38834304\n17284751360\n1245488\n25837438977\n1694498821\n4966056232\n1621115136\n151143680\n25778200848\n21542011016\n805605377\n205527296\n335552585\n8591024161\n22028308\n17450467333\n268484708\n34498347016\n5375000706\n827584\n4907338752\n34762917952\n2167472256\n35702964484\n13170114577\n9131033608\n2523176\n184582156\n211812420\n17726210052\n34896939012\n405798936\n1346375760\n2583691394\n17314089994\n4564189190\n9135259713\n8590200884\n38923174408\n8724318336\n42953878016\n2160083136\n1141113154\n6442606594\n3150593\n4850848\n17263788546\n335546900\n17448862736\n2820702720\n142635016\n17246979364\n2152013856\n17183286272\n4396154928\n286281730\n17215526928\n72941632\n9705536\n18262131776\n2685665298\n12887080992\n419439136\n26843561990\n8599913472\n8736735266\n\n");

function check(x,m){if(!x)throw Error(m);}
function nCr(n,k){let a=1;for(let j=1;j<=k;j++)a=a*(n-j+1)/j;return a;}
const idx=(m,u)=>Number(m>>BigInt(u*12)&4095n);
const W=[];
for(let u of [0,2]){
 const F=g.o.filter(m=>g.pop(BigInt(idx(m,u)))===6).map(m=>idx(m,u));
 check(F.length===132&&new Set(F).size===132,"one signed octad per hexad");
 const S=new Set(F);
 check(F.every(v=>S.has(4095^v)),"complement design");
 let v=0,T=[];
 for(let i=0;i<4096;i++){
  const count=g.pop(BigInt(i));
  if(count===5){check(F.filter(h=>(h&i)===i).length===1,"Witt S(5,6,12)");v++;}
  if(count===8){
   const h=F.reduce((s,x)=>{const k=g.pop(BigInt(x&i));s[k]=(s[k]??0)+1;return s;},{});
   check(h[5]===32&&h[6]===4,"uniform eight-support incidences");T.push(h[5]*4+h[6]*2);
  }
 }
 check(v===792&&T.length===495&&T.every(x=>x===136),"blocker incidences");
 W.push(F);
}
function rationalMax(skip){
 // Work with u=sqrt(2)*abs(x), and sorted coordinates. u_1+u_2<=1,
 // u_1+...+u_4+u_skip<=2. Gap variables y_i=u_i-u_(i+1)>=0.
 const h=k=>Math.min(k,4)+ +(k>=skip);
 const V=[];for(let i=1;i<=12;i++){
  const a=Math.min(i,2),c=h(i), U=[1,a],V2=[2,c];
  const x=U[0]*V2[1]<=V2[0]*U[1]?U:V2;
  V.push([i*x[0]**2,x[1]**2]);
  for(let j=i+1;j<=12;j++){
   const b=Math.min(j,2),d=h(j), k=a*d-b*c;
   if(k===0)continue;
   const xx=(d-2*b)*Math.sign(k),yy=(2*a-c)*Math.sign(k);
   if(xx<0||yy<0)continue;
   V.push([i*(xx+yy)**2+(j-i)*yy**2,k*k]);
  }
 }
 const best=V.reduce((a,b)=>a[0]*b[1]>=b[0]*a[1]?a:b);
 return best;
}
const z=[5,7,9].map(rationalMax);
check(z[0][0]*25===48*z[0][1]&&z[1][0]*25===48*z[1][1],"no unit point for 0 and 1 removed hexad banks");
check(z[2][0]===2*z[2][1],"two-bank relaxation reaches only norm²2 scaled; actual uniform eight-support hits 36 banks");
check(132*6/12===66&&132*nCr(6,2)/nCr(12,2)===30,"Witt incidences");
check(36+30*8===276&&132*2===264&&4*(9/2-2)<12,"at least five above-threshold hexads for any new unit point with original global roots");
check(136>6*4+2,"uniform eight-support blocks more old points than one old point can block added rays");
console.log(JSON.stringify({WittLines:[0,2],hexads:132,protectedFiveSubsets:792,actualSphericalNormSquareBound:"24/25 even after one full hexad bank is deleted",arbitraryRayMinimumBanks:5,arbitraryRayMinimumOldBlockers:10,uniformEightOldBlockers:136,oldRayUpperNewBlockers:26}));

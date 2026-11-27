// Node 20+, inputs (in order): original CT skeleton, triangle-refill certificate, double-trine witness, live neutral-reheat witness.
import {readFileSync} from "node:fs";
import {createHash} from "node:crypto";
const pins=["c4598a0dfc49365ae2b082745454c38b53b76703e5eef5dafbd71e2d06c59180","8ef795fd0b034e913e28c05b36f43c80cf274f78dd2af2bb8a07b1488bc8c306","ad81a8abdfb10889b7f7826b5f4823fa355f415e7cca4083d092167faed13b08","5327877dfd99b826787391c9ecdb76db9f58747ae693f8540ef03f0b3b678c93"];
const [C,T,D,W]=pins.map((pin,i)=>{const v=readFileSync(process.argv[i+2]);if(createHash("sha256").update(v).digest("hex")!==pin)throw Error("source "+i+" hash");return JSON.parse(v);});
function loadCT(C) {
 const G=C.gram_12, M=C.sublattice_matrix, n=12, range=n=>Array.from({length:n},(_,i)=>i);
 const need=(x,label)=>{if(!x)throw Error(label)};
 const mod=x=>(x%3+3)%3, pw=range(6).map(i=>3**i);
 const gc=(a,b)=>b?gc(b,a%b):Math.abs(a),lcm=(a,b)=>a/gc(a,b)*b;
 function F(a,b=1){if(b<0){a=-a;b=-b;}const g=gc(a,b);need(Number.isSafeInteger(a)&&Number.isSafeInteger(b)&&b>0,"rational exact");return [a/g,b/g];}
 const add=([a,b],[c,d])=>F(a*d+b*c,b*d),sub=([a,b],[c,d])=>F(a*d-b*c,b*d),mul=([a,b],[c,d])=>F(a*c,b*d),div=([a,b],[c,d])=>F(a*d,b*c);
 let L=range(n).map(i=>range(n).map(j=>F(+(i===j)))),D=[];
 for(let j=0;j<n;j++){
  let z=F(G[j][j]);for(let k=0;k<j;k++)z=sub(z,mul(D[k],mul(L[j][k],L[j][k])));
  need(z[0]>0,"PD");D.push(z);
  for(let i=j+1;i<n;i++){let t=F(G[i][j]);for(let k=0;k<j;k++)t=sub(t,mul(D[k],mul(L[j][k],L[i][k])));L[i][j]=div(t,z);}
 }
 const denL=range(n).map(j=>L.reduce((p,x)=>lcm(p,x[j][1]),1));
 const pp=range(n).map(j=>D[j][1]*denL[j]**2), J=pp.reduce(lcm,1);
 const Q=range(n).map(j=>J/pp[j]*D[j][0]), K=range(n).map(j=>range(n).map(i=>L[i][j][0]*(denL[j]/L[i][j][1])));
 need(Number.isSafeInteger(8*J),"global integer bound");
 for(let a=0;a<n;a++)for(let b=0;b<n;b++)need(K.reduce((s,k,i)=>s+Q[i]*k[a]*k[b],0)===J*G[a][b],"LDL equals G");
 const GP=v=>G.map(r=>r.reduce((s,g,j)=>s+g*v[j],0));
 const ip=(x,gy)=>x.reduce((s,v,j)=>s+v*gy[j],0);
 const label=v=>pw.reduce((s,p,i)=>s+p*mod(v[C.label_first[i]]+v[C.label_second[i]]),0);
 const R=v=>v.map((t,j)=>t-v.reduce((s,x,i)=>s+x*M[i][j],0));
 const neg=k=>pw.reduce((s,p)=>s+p*mod(-(Math.floor(k/p)%3)),0);
 const diff=(a,b)=>pw.reduce((s,p)=>s+p*mod(Math.floor(a/p)%3-Math.floor(b/p)%3),0);
 const v=range(n).fill(0), S=range(729).map(()=>[]),sh=range(9).map(()=>[]);
 const chosen=[];
 function isqrt(x){let y=Math.floor(Math.sqrt(x));while((y+1)**2<=x)y++;while(y*y>x)y--;return y;}
 function walk(j,rem){
  if(j<0){const u=v.slice(), N=ip(u,GP(u)), ell=label(u);need(J*N===8*J-rem,"valid leaf");sh[N].push(u);S[ell].push([N,u]);return;}
  const s=K[j].reduce((a,k,i)=>a+(i>j?k*v[i]:0),0), p=isqrt(Math.floor(rem/Q[j])), a=K[j][j];
  for(let t=Math.ceil((-p-s)/a);t<=Math.floor((p-s)/a);t++) {v[j]=t;walk(j-1,rem-Q[j]*(a*t+s)**2);}
 }
 walk(11,8*J);
 need(JSON.stringify(sh.map(a=>a.length))===JSON.stringify([1,0,0,0,756,0,4032,0,20412]),"K12 short shells");
 const mn=S.map(x=>Math.min(...x.map(u=>u[0]))), banks=S.map((x,i)=>x.filter(u=>u[0]===mn[i]).map(u=>u[1]));
 const reps=banks.map((b,i)=>{if(mn[i]!==4)return null;let a=b[0];return [a,R(a),R(R(a))];});
 for(let i=0;i<729;i++)if(mn[i]===4)need(new Set(reps[i].map(v=>v.join(","))).size===3&&reps[i].every(v=>label(v)===i&&ip(v,GP(v))===4),"well-formed phase numbering");
 return {G,M,ip,GP,label,R,neg,diff,sh,bank:banks,mn,reps,J,D,L,Q,K};
}
function currentTriples(C,X,T,D,W) {
 const {mn,reps}=X, need=(x,m)=>{if(!x)throw Error(m);};
 const ql=mn.flatMap((m,i)=>m===4?[i]:[]), ix=new Int16Array(729).fill(-1);ql.forEach((c,i)=>ix[c]=i);
 const root=ql.flatMap(c=>reps[c]), enc=(a,b,c)=>a+756*b+756*756*c;
 function code(r){let u=ix[r[0]],v=ix[r[1]],w=ix[r[2]];
  need(u>=0&&v>=0&&w>=0,"triple labels");return enc(u*3+r[3],v*3+r[4],w*3+r[5]);
 }
 const rootRows=new Set();
 function full(c){for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let k=0;k<3;k++)add([...c,i,j,k]);}
 function add(r){const k=code(r);need(!rootRows.has(k),"new triple duplicate");rootRows.add(k);}
 function del(r){need(rootRows.delete(code(r)),"missing old triple before removal");}
 for(const c of C.outer_cosets)if(c.every(k=>mn[k]===4))full(c);
 for(const r of T.removed_points)del(r);for(const r of T.added_points)add(r);
 for(const r of D.remove_rows)del(r);for(const c of D.insert_full_triples)full(c);
 need(rootRows.size===99225,"previous actual triple count");
 for(const r of W.remove_root)del(r);
 for(const [a,b,c,m] of W.add)for(let k=0;k<27;k++)if((m>>k)&1)add([a,b,c,k%3,k/3%3|0,k/9|0]);
 need(rootRows.size===106326,"actual current triple count");
 const V3=[...rootRows].map(v=>[v%756,v/756%756|0,v/(756*756)|0]);
 return {ql,ix,root,V3};
}
function orphanContactFlags(C,X,T,D,W) {
 const {ix,root,V3}=currentTriples(C,X,T,D,W),{mn,bank,ip,GP,neg}=X;
 const ng=root.map(GP), RR=root.map(q=>Int8Array.from(ng,w=>ip(q,w)));
 const buckets=Array.from({length:3},(_,h)=>Array.from({length:756},(_,c)=>[]));
 for(const x of V3)for(let h=0;h<3;h++)buckets[h][x[h]].push(x);
 const used=new Set(W.restore_base.map(([i,j])=>i+":"+j)), orph=W.remove_base.filter(([i,j,m])=>m===7&&!used.has(i+":"+j));
 if(orph.length!==1134)throw Error("orphan heavy directions");
 const H=new Uint32Array(16);
 const pair=Array.from({length:3},()=>Array.from({length:3},()=>new Int16Array(729).fill(-1)));
 const q6=mn.flatMap((m,i)=>m===6?[i]:[]), gp6=q6.map(i=>bank[i].map(GP));
 for(const c of C.outer_cosets)if(c.filter(i=>mn[i]===6).length===2&&c.includes(0)){
   const a=c.findIndex(i=>i!==0),b=c.findLastIndex(i=>i!==0);
   pair[a][b][c[a]]=c[b];pair[b][a][c[b]]=c[a];
 }
 for(const [owner,j] of orph){
  const c=C.outer_cosets[owner],h=c.findIndex(k=>mn[k]===8),p=c.findIndex(k=>mn[k]===4),e=c.indexOf(0),v=bank[c[h]][j],vv=GP(v),rx=3*ix[c[p]],ry=3*ix[neg(c[p])];
  if(owner<=0||rx<0||ry<0)throw Error("legal old labels");
  const near=q6.filter((c,i)=>gp6[i].some(w=>ip(v,w)===5));
  if(near.length!==32)throw Error("correct fixed q6 contact census");
  for(const [side,r] of [[p,root[rx]],[e,root[ry]]]){
    const g=GP(r),P=near.flatMap(c=>bank[pair[h][side][c]].map(w=>ip(w,g)));
    if(P.length!==576||Math.min(...P)!==-2||Math.max(...P)!==1)throw Error("expected positive and opposite root pressures");
  }
  const masks=new Uint8Array(9);
  for(let k=0;k<756;k++)if(ip(root[k],vv)===4)for(const x of buckets[h][k])for(let a=0;a<3;a++)for(let b=0;b<3;b++){
    const u=RR[rx+a][x[p]],w=RR[ry+b][x[e]];
    if(u>=2&&w>=1)masks[a+3*b]|=1;
    if(u>=1&&w>=2)masks[a+3*b]|=2;
    if(u>=2&&w<=-1)masks[a+3*b]|=4;
    if(u<=-1&&w>=2)masks[a+3*b]|=8;
  }
  for(let i=0;i<9;i++)H[masks[i]]++;
 }
 if(H[15]!==10206)throw Error("needed old 444 contact flags");
 return {baseHeavyOrphans:orph.length,tripleStock:V3.length,bothSignedAxes:10206,flagCounts:[...H]};
}

console.log(JSON.stringify(orphanContactFlags(C,loadCT(C),T,D,W)));

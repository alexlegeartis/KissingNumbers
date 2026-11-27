'use strict';
const assert = require('node:assert/strict');
const wt=x=>{x>>>=0; let n=0;while(x){x=(x&(x-1))>>>0;n++}return n;};
const ar=(n,f)=>Array.from({length:n},(_,i)=>f(i));
const xor=(a,k)=>a.reduce((x,v,j)=>x^((k>>j)&1?v:0),0);
const mon=I=>ar(16,i=>I.every(j=>(i>>j)&1)?1<<i:0).reduce((x,y)=>x|y,0);
const M=[[],[0],[1],[2],[3],[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]].map(mon);
const A=ar(32,i=>xor(M.slice(0,5),i)),Q=ar(64,i=>xor(M.slice(5),i));
const S=Q.map(q=>A.map(a=>a^q)),min=S.map(h=>Math.min(...h.map(wt)));
const bent=min.map(x=>x===6),B=Q.map(q=>Q.map(r=>wt(q&r)&1));
const V=ar(64,i=>i).filter(i=>bent[i]);
assert.equal(new Set(S.flat()).size,2048);assert.equal(V.length,28);
for(let i=0;i<64;i++)for(let j=0;j<64;j++){
  assert.equal(B[i][j],(Number(bent[i])+Number(bent[j])+Number(bent[i^j]))&1);
  assert(A.every(a=>!(wt(a&Q[i])&1)));
}
for(let i=1;i<64;i++)assert.equal(min[i],bent[i]?6:4);
// Exact finite shoulder-overlap test: adjacent banks must have bent difference.
const H=V.map(i=>S[i].filter(s=>wt(s)===6));assert(H.every(h=>h.length===16));
for(let i=0;i<28;i++)for(let j=i+1;j<28;j++){
  let mx=Math.max(...H[i].flatMap(s=>H[j].map(t=>wt(s&t))));
  assert.equal(mx,bent[V[i]^V[j]]?3:4);
}
// For each anisotropic plane, partition any affine slice into 32 Clebsch sections.
let nClebsch=0;
for(let ai=0;ai<28;ai++)for(let bi=ai+1;bi<28;bi++){
  const a=V[ai],b=V[bi];if(!B[a][b])continue;
  const sing=ar(64,i=>i).filter(i=>i&&!B[a][i]&&!B[b][i]&&!bent[i]);
  assert.equal(sing.length,5);assert.equal(new Set(ar(16,i=>xor(sing.slice(0,4),i))).size,16);
  assert.equal(xor(sing,31),0);
  const W=sing.map(i=>S[i].filter(w=>wt(w)===4));assert(W.every(x=>x.length===4));
  function pick(i,z){if(i===5)return !z;return W[i].some(v=>pick(i+1,z^v));}
  assert(pick(0,0));nClebsch++;
}
assert.equal(nClebsch,168);
const E=[];for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)if((i>>1)!==(j>>1))E.push([i,j]);
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),face=f=>ar(3,i=>i===(f>>1)?1-2*(f&1):0);
const normals=E.map(([a,b])=>face(a).map((x,i)=>x+face(b)[i]));
const T=[...ar(8,k=>ar(3,i=>1-2*((k>>i)&1))),...ar(6,face)];
const dirs=T.map(t=>normals.reduce((m,v,i)=>dot(t,v)>0?m|(1<<i):m,0));
assert(dirs.slice(0,8).every(m=>wt(m)===3));assert(dirs.slice(8).every(m=>wt(m)===4));
const opp=E.map(([a,b])=>E.findIndex(([c,d])=>a===(c^1)&&b===(d^1)));
const tags=ar(32,t=>E.reduce((m,[a,b],i)=>((t<<1>>a)^(t<<1>>b))&1?m|(1<<i):m,0));
const rays=ar(4096,h=>dirs.reduce((m,d,i)=>h&d?m:m|(1<<i),0));
const eligible=ar(4096,p=>tags.reduce((m,d,i)=>rays[p^d]?m|(1<<i):m,0));
assert.equal(Math.max(...eligible.map(wt)),24);
const near=T.map(t=>T.reduce((m,u,i)=>{const d=dot(t,u);return d>0&&4*d*d>dot(t,t)*dot(u,u)?m|(1<<i):m;},0));
const nearAll=ar(16384,m=>near.reduce((a,l,i)=>m&(1<<i)?a&l:a,16383));
let nFaces=0, exceptional={}, matched=0, diag=0, square=0, radialMatched=0, rank2lines=0, radialParallel=0, radialDisjoint=0, dm={};
for(const y of V) for(const z of V)if(bent[12^z]&&bent[y^z])
for(const w of V)if(bent[12^w]&&bent[y^w])
for(let x=0;x<64;x++)if(bent[x^12]&&bent[x^y]&&bent[x^z]&&bent[x^w]){
  const f=[0,x,12,y,z,w],labels=E.map(([a,b])=>f[a]^f[b]);nFaces++;
  const Eq=[];
  for(let i=0;i<12;i++)for(let j=i+1;j<12;j++)if(labels[i]===labels[j]&&opp[i]!==j)Eq.push((1<<i)|(1<<j));
  const tt=B.map(row=>f.slice(1).reduce((r,a,j)=>r|(row[a]<<j),0));
  const P=tt.reduce((r,t)=>r|(1<<t),0), rank=Math.log2(wt(P)),ev=tt.map(t=>tags[t]),gauges=[...new Set(ev)];
  assert(Number.isInteger(rank)&&rank>=2&&rank<=5); if(rank===5)continue;
  for(let p=0;p<4096;p++){
    const N=wt(eligible[p]&P)*(64/wt(P));
    if(N<=48||Eq.some(e=>wt(p&e)!==1)||gauges.some(g=>(g^p)<p))continue;
    const D=ev.map(e=>rays[p^e]);assert.equal(D.filter(Boolean).length,N);
    let mt=0, allSingletonCharges=true;const used=new Set();
    for(let i=0;i<64;i++)if(D[i]&&!used.has(i)){
      for(let j=i+1;j<64;j++)if(D[j]&&!used.has(j)&&!bent[i^j]&&!(D[j]&~nearAll[D[i]])){
        used.add(i);used.add(j);mt++;allSingletonCharges&&=D[i]===D[j]&&D[i]<256&&wt(D[i])===1;break;
      }
    }
    const key=rank+'/'+N+'/'+mt;exceptional[key]=(exceptional[key]??0)+1;
    if(N-mt<=48){matched++;if(allSingletonCharges)radialMatched++;continue;}
    if(rank<=3){
      const edges=new Set();
      for(const d of D){
        assert(d<256);if(wt(d)<=1)continue; assert.equal(wt(d),2);
        const [i,j]=ar(8,i=>i).filter(i=>d&(1<<i));
        if((i^j)===7)continue;assert.equal(wt(i^j),1);edges.add(d);
      }
      const list=[...edges];function owner(i,a){return i===list.length||ar(8,j=>j).some(j=>(list[i]&(1<<j))&&!(a&(1<<j))&&owner(i+1,a|(1<<j)));}
      assert(owner(0,0));diag++;
      
      let activeAxis=-1;const axes=new Set();
      for(let q=0;q<64;q++){
        const d=D[q];assert.equal(wt(d),2);
        const [i,j]=ar(8,i=>i).filter(i=>d&(1<<i));
        const pbad=p^ev[q];assert.equal(wt(pbad),6);
        if((i^j)===7){
          assert(normals.every((v,k)=>!!(pbad&(1<<k))===(dot(v,T[i])===0)));
        }else{
          const ax=Math.log2(i^j);assert(Number.isInteger(ax));
          axes.add(ax);activeAxis=ax;
          const normalSet=normals.map(v=>v.map((x,k)=>k===ax?x:x*T[i][k]));
          // Two transverse difference normals and both x signs toward each negative transverse coordinate.
          assert(normalSet.every((v,k)=>!!(pbad&(1<<k))===(v[ax]===0?v.reduce((a,b)=>a+b,0)===0:v.some((z,h)=>h!==ax&&z===-1))));
        }
      }
      if(rank===2){assert(activeAxis<0);rank2lines++;}
      else {assert(activeAxis>=0);if(axes.size===1)radialParallel++;const es=[...edges]; 
let deg=0;for(let b=0;b<8;b++)deg=Math.max(deg,es.filter(x=>(x>>b)&1).length);
dm[es.length+'/'+deg]=(dm[es.length+'/'+deg]??0)+1;if(deg<=1)radialDisjoint++;}

      continue;
    }
    assert.equal(rank,4);
    const U=D.reduce((a,b)=>a|b,0)>>8,ax=[0,1,2].find(i=>U===(3<<(2*i)));assert(ax!==undefined);
    const caps=[0,1].map(s=>ar(8,i=>i).reduce((m,i)=>m|(((i>>ax)&1)===s?1<<i:0),1<<(8+2*ax+s)));
    for(const d of D){
      assert(!(d&~(caps[0]|caps[1])));
      const a=d&caps[0],b=d&caps[1]; if(a&&b){
        assert(wt(a)===1&&wt(b)===1&&a<256&&b<256&&(Math.log2(a)^Math.log2(b))===7);
      } else if(d>=256) assert(d&255);
    }
    square++;
  }
}
assert.equal(nFaces,3168);assert.equal(diag,552);assert.equal(matched,3360);assert.equal(square,3600);
assert.deepEqual(exceptional,{"2/64/0":12,"3/64/0":540,"3/64/16":2160,"3/64/32":1200,"4/56/0":3600});
assert.equal(radialMatched,3360);assert.equal(rank2lines,12);assert.equal(radialDisjoint,360);assert.deepEqual(dm,{'4/1':360,'8/2':180});
console.log(JSON.stringify({nFaces,nClebsch,exceptional,diag,matched,square,radialMatched,rank2lines,radialParallel,radialDisjoint,dm,allRealClosedFlowCap:1536,totalWithSparse:12204}));

{
// In the notation of the closed-flow verifier: enumerate all two-adjacent-defect labels.
const adj=E.map((_,i)=>ar(12,j=>j).filter(j=>j!==i&&dirs.slice(0,8).some(d=>(d>>i&1)&&(d>>j&1))));
const seen=new Set(),cache=new Map(),sources={},st={};
let Nprofiles=0,Ncandidate=0,Nphases=0;
for(const y of V)for(const z of V)if(bent[12^z]&&bent[y^z])
for(const w of V)if(bent[12^w]&&bent[y^w])
for(let x=0;x<64;x++)if(bent[x^12]&&bent[x^y]&&bent[x^z]&&bent[x^w]){
  const f=[0,x,12,y,z,w],L=E.map(([a,b])=>f[a]^f[b]);
  for(let e=0;e<12;e++)for(const z2 of V)if(z2!==L[e]&&adj[e].every(j=>bent[z2^L[j]])){
    const labels=L.slice();labels[e]=z2;
    const key=labels.join('.');if(seen.has(key))continue;seen.add(key);Nprofiles++;
    const Eq=[];let consistent=true;
    for(let i=0;i<12;i++){
      const same=ar(12,j=>j).filter(j=>j!==i&&labels[i]===labels[j]&&opp[i]!==j);
      for(const j of same)if(j>i)Eq.push((1<<i)|(1<<j));
      if(same.length>2)consistent=false;
    }
    if(!consistent)continue; // any three same-label non-antipodal neighbours give an odd triangle.
    const ev=B.map(row=>labels.reduce((r,l,i)=>r|(row[l]<<i),0));
    const im=[...new Set(ev)].sort((a,b)=>a-b),ik=im.join('.');let phases=cache.get(ik);
    if(phases===undefined){
      phases=[];const vis=new Uint8Array(4096);
      for(let p=0;p<4096;p++)if(!vis[p]){
        let n=0;for(const g of im){const h=g^p;vis[h]=1;if(rays[h])n++;}
        const N=n*64/im.length;if(N>56)phases.push([N,p]);
      }
      cache.set(ik,phases);
    }
    if(!phases.length)continue;Ncandidate++;
    for(const [N,p] of phases){
      if(Eq.some(m=>wt(p&m)!==1))continue;Nphases++;
      st[N]=(st[N]??0)+1;
      const D=ev.map(g=>rays[p^g]);assert.equal(D.filter(Boolean).length,N);
      let mt=0;const used=new Set();
      for(let i=0;i<64;i++)if(D[i]&&!used.has(i))
        for(let j=i+1;j<64;j++)if(D[j]&&!used.has(j)&&!bent[i^j]&&!(D[j]&~nearAll[D[i]])){
          used.add(i);used.add(j);mt++;break;
        }
      if(mt>=N-56){sources.match=(sources.match??0)+1;continue;}
      const edge=new Set();
      for(const d of D){
        assert(d<256);if(wt(d)<=1)continue;assert.equal(wt(d),2);
        const [i,j]=ar(8,i=>i).filter(i=>d&(1<<i));
        if((i^j)===7)continue;assert.equal(wt(i^j),1);edge.add(d);
      }
      const ls=[...edge];function own(i,m){return i===ls.length||ar(8,b=>b).some(b=>(ls[i]&(1<<b))&&!(m&(1<<b))&&own(i+1,m|(1<<b)));}
      assert(own(0,0));sources.diag=(sources.diag??0)+1;
    }
  }
}
assert.equal(Nprofiles,86400);assert.equal(cache.size,516);
assert.equal(Ncandidate,56160);assert.equal(Nphases,52560);
assert.deepEqual(st,{'60':11520,'64':41040});assert.deepEqual(sources,{diag:14400,match:38160});
console.log(JSON.stringify({Nprofiles,Ncandidate,Nphases,st,sources,twoAdjacentDefectsRealTailCap:1792}));
}
{
const sub=(a,b)=>a.map((x,i)=>x-b[i]);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
let strict=0,equal=0;
for(let i=0;i<12;i++)for(let j=i+1;j<12;j++)for(let k=j+1;k<12;k++){
  const a=normals[i],b=normals[j],c=normals[k],N=cross(sub(b,a),sub(c,a)),h=dot(N,a);
  if(!h)continue;assert(5*h*h>=dot(N,N));strict++;if(5*h*h===dot(N,N))equal++;
}
assert.equal(strict,128);assert.equal(equal,24);
console.log(JSON.stringify({facetsThroughThreeNonzero:strict,sharp:equal,minimumPositiveSquaredPlaneDistance:"1/5"}));
}
{
  const f=[0,1,12,12,22,23],p=798,L=E.map(([a,b])=>f[a]^f[b]),g=B.map(row=>L.reduce((h,v,i)=>h|(row[v]<<i),0));
  const D=g.map(h=>rays[h^p]),Kd=ar(64,q=>q).filter(q=>!g[q]);
  assert.deepEqual(Kd,[0,1,14,15,18,19,28,29]);assert.deepEqual(Kd.filter(q=>q&&!bent[q]),[1]);
  assert.deepEqual([...new Set(D)].sort((a,b)=>a-b),[3,12,17,34,48,68,136,192]);
  for(let q=0;q<64;q++) for(const d of new Set(D)) if(d!==D[q]&&(d&D[q])){
    assert.equal(ar(64,t=>t).filter(t=>D[t]===d&&bent[q^t]).length,4);
  }
  console.log(JSON.stringify({residualFaces:f,residualPhase:p,residualKernel:Kd}));
}
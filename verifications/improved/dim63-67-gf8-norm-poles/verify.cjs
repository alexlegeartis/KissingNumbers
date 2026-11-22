"use strict";
const fs = require("node:fs"), path = require("node:path"), crypto = require("node:crypto");
function pin(rel, expected) {
    const b=fs.readFileSync(path.join(__dirname,"..",rel));
    const sha=crypto.createHash("sha1").update("blob "+b.length+"\0").update(b).digest("hex");
    if(sha!==expected) throw Error("Parent blob mismatch: "+rel);
    return JSON.parse(b);
}
const old63=pin("dim63-balanced-nine-refill/claims.json","e85ed9d87813b767a188823fbd144a6a7b355d1f");
const old64=pin("dim64-65-partial-nonary-signs/claims.json","96a743ef12f41a0c144c7ca9c86cf2c52749f573");
const old65=pin("dim65-67-bch-exterior-caps/claims.json","1dafcee005375381c2dad801cf50d942cfde1afc");
const old6667=pin("dim66-67-netto-signed-trades/claims.json","c05a97f686b77e669f784d66e51442fbb4143a44");
const oldSection=pin("dim63-65-affine-reconstruction/section63.json","cdc01424d1789116ac83c5bb6509d88269bfc16b");

const assert=(x,msg)=>{if(!x)throw Error(msg);};
const bit=i=>1n<<BigInt(i);
function wt(x){let k=0;while(x){x&=x-1n;k++;}return k;}
function mul(x,y){let z=0;while(y){if(y&1)z^=x;x<<=1;if(x&64)x^=0x43;y>>=1;}return z;}
function bch(d){
    const powers=[1];for(let i=1;i<63;i++)powers.push(mul(powers.at(-1),2));
    assert(new Set(powers).size===63&&mul(powers.at(-1),2)===1,"field order");
    const roots=new Set();
    for(let a=1;a<d;a++)for(let r=a;!roots.has(r);r=(2*r)%63)roots.add(r);
    let g=[1];
    for(const r of roots){
        const h=Array(g.length+1).fill(0);
        g.forEach((a,i)=>{h[i]^=mul(a,powers[r]);h[i+1]^=a;});
        g=h;
    }
    assert(g.every(x=>x===0||x===1),"binary generator");
    const poly=g.reduce((w,a,i)=>w^(a?bit(i):0n),0n);
    const rows=Array.from({length:64-g.length},(_,i)=>{let w=poly<<BigInt(i);return w|((wt(w)%2)?bit(63):0n);});
    return {rows,poly,powers};
}
function rref(words, width){
    const rows=[...words], piv=[];
    for(let col=0;col<width;col++){
        const n=piv.length;
        const index=rows.findIndex((r,i)=>i>=n&&(r&bit(col)));
        if(index<0)continue;
        [rows[n],rows[index]]=[rows[index],rows[n]];
        for(let k=0;k<rows.length;k++)if(k!==n&&(rows[k]&bit(col)))rows[k]^=rows[n];
        piv.push(col);
    }
    return {rows:rows.slice(0,piv.length),piv};
}
function residual(word,basis){
    for(let i=0;i<basis.piv.length;i++)if(word&bit(basis.piv[i]))word^=basis.rows[i];
    return word;
}
const span=rows=>rows.reduce((W,row)=>[...W,...W.map(x=>x^row)],[0n]);
const P=bch(13), E=bch(21);
const basisP=rref(P.rows,64);
assert(basisP.rows.length===30&&rref(E.rows,64).rows.length===18,"BCH ranks");
assert(E.rows.every(w=>residual(w,basisP)===0n),"nested extended BCH");
const gram=P.rows.map(a=>P.rows.reduce((v,b,i)=>v^((wt(a&b)%2)?bit(i):0n),0n));
const gb=rref(gram,30);
const nulls=Array.from({length:30},(_,j)=>j).filter(j=>!gb.piv.includes(j)).map(j=>
    gb.rows.reduce((w,r,i)=>w^((r&bit(j))?bit(gb.piv[i]):0n),bit(j)));
const D=rref(nulls.map(w=>P.rows.reduce((u,r,i)=>u^((w&bit(i))?r:0n),0n)),64);
assert(D.rows.length===28&&D.rows.every(a=>wt(a)%4===0&&D.rows.every(b=>wt(a&b)%2===0)),"doubly even radical");
const outside=P.rows.find(w=>residual(w,D)!==0n);
assert(typeof outside==="bigint","an exterior BCH translate");

const X=[...P.powers,0];
function gfMul(x,y,modulus){
    let s=0;const top=1<<(31-Math.clz32(modulus));
    while(y){if(y&1)s^=x;y>>=1;x<<=1;if(x&top)x^=modulus;}return s;
}
const pairs = [];
for (let i=0;i<6;i++) for (let j=i+1;j<6;j++) pairs.push([i,j]);
const points=X;
const linear = Array.from({length:6},(_,i) => points.reduce(
    (w,x,j)=>w ^ ((x&(1<<i)) ? bit(j) : 0n),0n));
const affineBasis = [(1n<<64n)-1n,...linear];
const quadraticBasis = pairs.map(([i,j])=>linear[i]&linear[j]);
const rmBasis = [...affineBasis,...quadraticBasis];
assert(rref(rmBasis,64).rows.length===22 && rmBasis.every(w=>residual(w,D)===0n), 'RM(2,6) embedding');
const affine = affineBasis.reduce((list,w)=>[...list,...list.map(x=>x^w)],[0n]);
assert(affine.length===128 && Math.min(...affine.slice(1).map(wt))===32, 'Affine distance');
function trace(x) {
    let t=0;
    for(let i=0;i<5;i++) {
        t ^= x;
        x=gfMul(x,x,0x25);
    }
    assert(t===0 || t===1, 'GF(32) trace');
    return t;
}
for(let x=1;x<32;x++)
    assert(new Set(Array.from({length:32},(_,y)=>gfMul(x,y,0x25))).size===32, 'GF(32)');
const forms = Array.from({length:32},(_,a)=>{
    const B=Array.from({length:5},(_,i)=>Array.from({length:5},(_,j)=>
        trace(gfMul(a,gfMul(1<<i,1<<j,0x25),0x25))));
    return pairs.reduce((m,[i,j],k)=>m | ((j===5 ? B[i][i] : (B[i][j]^(B[i][i]&B[j][j])))<<k),0);
});
assert(new Set(forms).size===32, 'Kerdock forms');
const difference = new Set(forms.flatMap(a=>forms.map(b=>a^b)));
difference.delete(0);
const evaluate = m => quadraticBasis.reduce((w,q,i)=>w ^ ((m&(1<<i)) ? q : 0n),0n);
let distance=32;
for(const m of difference) {
    const w=evaluate(m);
    for(const a of affine) distance=Math.min(distance,wt(w^a));
}
assert(distance===28, 'Kerdock distance');
const offsets=[], used=new Set();
for(let d=0;d<32768 && offsets.length<4;d++) {
    const translated=forms.map(m=>m^d);
    if(translated.some(m=>used.has(m))) continue;
    offsets.push(d);
    translated.forEach(m=>used.add(m));
}
assert(offsets.length===4 && used.size===128 && offsets.slice(0,3).join(",")==="0,1,2", 'Disjoint quadratic translates');

const one=(1n<<64n)-1n, eb=rref(E.rows,64);
const oldCaps=new Set(offsets.slice(0,3).flatMap(t=>
    forms.flatMap(q=>affine.map(a=>outside^evaluate(q^t)^a))));
assert(oldCaps.size===12288&&[...oldCaps].every(w=>residual(w,basisP)===0n),"old cap membership");
function norm(x){
    let y=x;for(let i=0;i<3;i++)y=mul(y,y);
    return mul(y,x);
}
const subfield=X.filter(x=>norm(x)===mul(x,x));
assert(subfield.length===8,"subfield size");
const fieldOctad=X.reduce((h,x,i)=>h|(subfield.includes(x)?bit(i):0n),0n);
assert(fieldOctad===0x8040201008040201n&&D.rows.every(w=>wt(w&fieldOctad)%2===0),"dual octad");
const physicalOctad=BigInt(oldSection.octad);
assert(physicalOctad===9606450284331136n&&oldSection.terms.poles===112&&oldSection.terms.dense===146800640,"pinned section terms");
const indices=mask=>X.map((_,i)=>i).filter(i=>(mask&bit(i))!==0n);
const permutation=new Array(64);
for(const flag of [true,false]){
    const from=indices(flag?fieldOctad:one^fieldOctad), to=indices(flag?physicalOctad:one^physicalOctad);
    from.forEach((i,j)=>permutation[i]=to[j]);
}
function permute(w){return indices(w).reduce((m,i)=>m|bit(permutation[i]),0n);}
assert(permute(fieldOctad)===physicalOctad&&new Set(permutation).size===64,"section coordinate map");
const projection=new Map(), nullRows=[];
for(let w of D.rows){
    for(let k of indices(fieldOctad)){
        if((w&bit(k))===0n)continue;
        if(projection.has(k))w^=projection.get(k);else{projection.set(k,w);w=0n;break;}
    }
    if(w!==0n)nullRows.push(w);
}
assert(projection.size===7&&nullRows.length===21&&rref(nullRows,64).rows.length===21&&
    nullRows.every(w=>(w&fieldOctad)===0n),"dense section projection and kernel");
assert(70*2**nullRows.length===146800640,"section dense count");
function chooseProjection(target){
    let w=0n;
    for(const k of indices(fieldOctad)){
        if(((w^target)&bit(k))===0n)continue;
        assert(projection.has(k),"section projection consistency");
        w^=projection.get(k);
    }
    assert((w&fieldOctad)===target,"section target");
    return w;
}
const phase=span(rref(P.rows.map(w=>residual(w,D)),64).rows)
    .find(w=>w!==0n&&wt(w&fieldOctad)%2===0);
assert(typeof phase==="bigint"&&residual(phase,basisP)===0n&&residual(phase,D)!==0n,
    "exterior section phase");
function simplex(heavy){
    const labels=X.map(x=>norm(x^X[heavy]));
    const levels=new Set(labels);
    assert(levels.size===8&&labels[heavy]===0,"norm label range");
    assert([...levels].every(v=>labels.filter(x=>x===v).length===(v===0?1:9)),"norm multiplicities");
    const gens=rref(Array.from({length:6},(_,i)=>labels.reduce(
        (v,x,j)=>v|((x&(1<<i))?bit(j):0n),0n)),64).rows;
    assert(gens.length===3&&gens.every(w=>residual(w,D)===0n),"norm rows are inside the radical");
    const group=span(gens);
    assert(group.slice(1).every(w=>wt(w)===36)&&group.every(w=>(w&bit(heavy))===0n),
        "one marked sign and distance 36");
    return [...group,...group.map(w=>w^one)];
}
function sphere(section){
    const shifts=span((section?nullRows:D.rows).slice(0,12));
    const taken=new Set(), out=[];
    const active=section?indices(one^fieldOctad):indices(one);
    for(const i of active){
        const group=simplex(i);
        let translated=outside;
        if(section){
            // Outside the subfield, the norm takes four values twice each on the subfield.
            const fibers=new Map();
            for(const j of indices(fieldOctad)){
                const label=norm(X[j]^X[i]);
                fibers.set(label,[...(fibers.get(label)??[]),j]);
            }
            assert(fibers.size===4&&[...fibers.values()].every(a=>a.length===2),"section pairs");
            const target=[...fibers.values()].reduce((m,[j])=>m|bit(j),0n);
            translated=phase^chooseProjection(target^(phase&fieldOctad));
        }
        let accepted=false;
        for(const shift of shifts){
            const batch=group.map(w=>w^translated^shift);
            if(batch.some(w=>taken.has(w)||(!section&&(residual(w,eb)===0n||oldCaps.has(w)))))continue;
            assert(batch.every(w=>residual(w,basisP)===0n&&residual(w,D)!==0n),"exterior code membership");
            assert(!section||batch.every(w=>wt(w&fieldOctad)===4&&wt(permute(w)&physicalOctad)===4),"63D balance in original coordinates");
            for(const w of batch){taken.add(w);out.push({i,w});}
            accepted=true;break;
        }
        assert(accepted,"a disjoint group exists for every marked coordinate");
    }
    assert(out.length===active.length*16&&taken.size===out.length,"spike cardinality");
    return out;
}
const full=sphere(false), section=sphere(true);
function replay(points){
    const norm=18**2+63*2**2, limit=norm/2;
    let pairs=0, maximum=-Infinity;
    for(let j=0;j<points.length;j++)for(let k=0;k<j;k++){
        const {i,w}=points[j], {i:l,w:v}=points[k], diff=w^v;
        let d;
        if(i===l)d=((diff&bit(i))?-1:1)*324+4*(63-2*wt(diff&(one^bit(i))));
        else{
            const special=bit(i)|bit(l), m=wt(diff&special);
            d=36*(2-2*m)+4*(62-2*wt(diff&(one^special)));
        }
        maximum=Math.max(maximum,d);assert(d<=limit,"new pair dot");pairs++;
    }
    assert(norm===576&&maximum===288,"norm and tight pair limit");
    return {points:points.length,pairs,scaledMaxDot:maximum};
}
const resultFull=replay(full), resultSection=replay(section);
// Scale actual spikes by 3 while old sparse and dense points keep norm squared 64.
assert(2*(18+15*2)===96&&4*(18+3*2)===96&&8*(18+8*2)<3*96,"old sparse envelopes");
assert(8*2<96,"ordinary poles on the other coordinates");
assert(18+2*(63-2*14)<96&&-18+2*(63-2*13)<96,"flat bodies from P");
const bases=[old63,old64,old65,old6667,old6667].map((stock,i)=>stock.dimensions.find(x=>x.dimension===63+i).lower_bound);
const results=bases.map((previous,i)=>{
    const replaced=i===0?112:128, inserted=i===0?section.length:full.length;
    return {dimension:63+i,previous,replaced,inserted,gain:inserted-replaced,lower_bound:previous+inserted-replaced};
});
assert(JSON.stringify(results)===JSON.stringify(JSON.parse(fs.readFileSync(path.join(__dirname,"claims.json"))).dimensions),"claims ledger");
console.log(JSON.stringify({full:resultFull,section:resultSection,sectionPhase:phase.toString(16),projectionRank:projection.size,results},null,2));



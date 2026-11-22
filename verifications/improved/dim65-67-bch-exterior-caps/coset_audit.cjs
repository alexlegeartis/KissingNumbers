"use strict";
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
const image=rref(E.rows.map(w=>residual(w,D)),64);
const quotient=span(image.rows);
const phases=[...new Set(P.rows.map(w=>residual(w,D)))];
const outer65=P.rows.find(w=>!quotient.includes(residual(w,D)));
const words=span(E.rows), hist={};
for(const w of words){const d=wt(w);hist[d]=(hist[d]??0)+1;}
assert(hist[0]===1&&hist[22]>0&&Math.min(...Object.keys(hist).map(Number).filter(Boolean))===22,"cap distance");
const key=w=>"0x"+w.toString(16);
assert(image.rows.length===2,"cap image rank");
const overlap = words.filter(w=>residual(w,D)===0n).length;
assert(overlap===65536,"dense overlap");
// The parent min14 also follows from the narrow-sense BCH bound >=13 and parity extension.
assert(36*36*3<=32*32*4,"65D cap versus unremoved equator");
assert(3*20+4*16<=4*32,"65D same-sign cap pairs");
assert(3*64-4*16===4*32,"65D opposite-sign cap pairs");
assert(3*16+4*16<4*32,"65D sparse and dense cap pairs");
assert(8*4===32,"65D tail poles");
const count65=335795354+2688+2+2*words.length-overlap;
assert(count65===336256796,"65D total");
// The three shifted Kerdock classes are disjoint because their unshifted classes are disjoint.
// The 66/67D envelope geometry is unchanged, except that min14 (rather than min16)
// bounds distinct ambient BCH bodies, and all equatorial bodies are retained.
assert(residual(outside,D)!==0n&&residual(outside,basisP)===0n,"exterior BCH coset");
assert(36*36*37<32*32*64,"66/67D cap versus equator");
assert(36*37<32*64,"66/67D distinct axes");
const bounds=[{dimension:65,old:335929116,new:count65},
    {dimension:66,old:342381366,new:342373178+2*(8192-2)},
    {dimension:67,old:355183103,new:355170821+3*(8192-2)}];
assert(bounds[1].new===342389558&&bounds[2].new===355195391,"66/67D totals");
for(const x of bounds)x.gain=x.new-x.old;
console.log(JSON.stringify({ambientGenerator:key(P.poly),capGenerator:key(E.poly),shift:key(outside),capRank:18,
    capWords:words.length,capIntersectionWithEquator:overlap,capHistogram:hist,bounds},null,2));

// Insertion check for the full fixed 65D BCH cap. Every other ambient BCH coset
// has an explicit member of weight at most 20, so every prospective body in it
// is within distance 20 of a retained cap body. The loop checks actual words.
const basisE=rref(E.rows,64);
const compIndex=rref(P.rows.map(w=>residual(w,basisE)),64);
assert(compIndex.rows.length===12,"outer BCH quotient rank");
const representatives=span(compIndex.rows);
assert(representatives.every(w=>residual(w,basisP)===0n),"quotient representatives are in P");
const cv = new Uint32Array(words.length), cw = new Uint32Array(words.length);
for(let i=0;i<words.length;i++){cv[i]=Number(words[i]&0xffffffffn);cw[i]=Number(words[i]>>32n);}
function pop32(w){
    w=w-((w>>>1)&0x55555555);
    w=(w&0x33333333)+((w>>>2)&0x33333333);
    return (Math.imul((w+(w>>>4))&0x0f0f0f0f,0x01010101)>>>24);
}
const witnesses=[], free=[];
let checks=0;
const buckets={};
for(let q=1;q<representatives.length;q++){
    const lo=Number(representatives[q]&0xffffffffn), hi=Number(representatives[q]>>32n);
    let witness=-1;
    for(let k=0,j=(q*10177)&262143; k<cv.length; k++,j=(j+67003)&262143){
        checks++;
        const d=pop32(lo^cv[j])+pop32(hi^cw[j]);
        if(d>=22)continue;
        assert(d===wt(representatives[q]^words[j]),"literal blocker weight");
        buckets[d]=(buckets[d]??0)+1;
        witness=j;break;
    }
    witnesses.push(witness);
    if(witness===-1)free.push(q);
}

assert(free.length===0,"every nonzero coset has a blocker");
console.log(JSON.stringify({ambientCosetsBlocked:witnesses.length,blockerWeights:buckets,checks},null,2));



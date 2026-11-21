// Run node flat67.cjs patterns.json weight10.json patterns_extra.json.
// Obtain the three pinned input files listed in the heartbeat checkpoint.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const [patterns, weight10, extra] = process.argv.slice(2).map(path => JSON.parse(fs.readFileSync(path, "utf8")));

const six = [[3984,2816],[4277,164],[7938,3584],[19978,17416],[27681,25600],[29444,17152],[34891,74],[35372,33312],[37570,33408],[38153,33032],[53262,53248],[66243,65728],[66339,66050],[66792,1120],[67667,67586],[67747,162],[84998,83968],[90258,82048],[98644,336],[105474,70656],[131470,140],[132677,1540],[137568,6208],[141443,10368],[148130,672],[164138,131336],[165010,163856],[188488,49216],[197894,131332],[201240,4624],[231553,165888],[278977,278592],[303111,303104],[311818,33288],[338184,272384],[380960,53248],[393442,393248],[394770,394240],[401520,262240],[406656,13312],[411668,133136],[491554,360448],[525510,525316],[534724,526400],[545872,17472],[557258,32840],[575746,526592],[592899,527360],[602626,598016],[606221,606208],[632320,558080],[655690,524608],[663825,655616],[666632,656384],[721414,66052],[787748,263424],[793089,786944],[795394,270848],[819352,557072],[1051810,1049728],[1054790,1052676],[1056843,1056776],[1057417,648],[1058114,1049664],[1065272,16656],[1070596,1053184],[1091602,34832],[1098884,1065984],[1133120,68096],[1180722,1072],[1191939,139266],[1253412,1245184],[1310890,168],[1311366,1049092],[1314892,1310728],[1328642,279040],[1377297,1311744],[1475328,1441792],[1573398,524308],[1642880,1052928],[1671201,1638400],[1909248,1376256],[2097351,70],[2097510,2097220],[2098826,2097672],[2106396,9224],[2119816,2099328],[2138528,41216],[2146833,49168],[2175041,2101312],[2179170,2162720],[2230283,133122],[2238978,139776],[2294026,2293760],[2311424,148480],[2363666,2359552],[2363914,266248],[2376714,2360320],[2395200,265216],[2424839,2162692],[2425476,2097792],[2490409,2359328],[2622280,2621952],[2622626,524448],[2627590,2099204],[2639906,2637824],[2689072,2162720],[2789380,2654208],[3148037,1050880],[3151392,1050112],[3244042,2129928],[3277008,3145792],[3295234,1067008],[3956736,802816]];
const seven = [[799,772],[1905,1344],[11572,3076],[12760,200],[14365,10256],[17014,16912],[19788,16452],[21097,16960],[38244,4356],[43841,32833],[46736,1552],[53523,32785],[71738,65554],[74318,8768],[75141,65668],[78242,12416],[103473,4144],[115592,896],[145572,4256],[147673,131264],[172932,8576],[180517,49408],[181416,33800],[201924,65604],[207952,3088],[215592,131592],[263227,262168],[282663,4130],[286892,8228],[288544,8960],[292930,24578],[296609,32801],[301444,37120],[329829,327744],[329884,262168],[333632,4672],[344372,262164],[348816,279040],[395976,262784],[492616,262216],[526761,526464],[530060,5128],[542593,1537],[549217,540673],[559189,557072],[561968,525056],[566313,1064],[585857,49153],[606490,65794],[612640,528416],[615072,590336],[655484,88],[661650,655488],[668200,143360],[793640,266272],[815124,794624],[836164,262660],[871552,264320],[888897,98368],[1049572,832],[1054296,1048656],[1071440,4416],[1073237,8197],[1085989,1048609],[1090574,9220],[1128705,2305],[1143816,1057792],[1156368,66816],[1188632,139776],[1204672,8384],[1218696,1050632],[1294402,1081346],[1321505,1049089],[1343875,262401],[1361929,1343488],[1384504,1048624],[1409616,294928],[1445128,263176],[1459296,131168],[1507744,1114144],[1582256,8336],[1589411,1064992],[1706081,524384],[1761284,1572868],[1771544,1572872],[1837330,1048848],[1847488,270400],[1872898,528386],[1966213,655361],[1986568,1589248],[2107992,2624],[2121995,2097162],[2123456,2097344],[2131718,518],[2136272,2101376],[2166281,2057],[2167084,69640],[2205736,8232],[2212016,32928],[2233105,4353],[2256944,2244608],[2263320,163848],[2327136,197120],[2362208,2099456],[2369681,2361344],[2372704,2098208],[2441304,2097224],[2498884,2367488],[2507781,393220],[2630181,2621472],[2637056,2629632],[2638228,540800],[2639945,524353],[2658329,2134016],[2753681,525328],[2826376,196616],[2924610,532544],[2982408,295424],[3063936,49280],[3158292,3149824],[3199072,2129984],[3231873,1065088],[3237904,1114128],[3293708,524],[3318785,2105345],[3440692,2359328],[3610624,137216],[3703426,3145856],[3736644,2162692],[3809856,1581056],[3933441,1310721]];

function weight(mask) {
  let count = 0;
  for (; mask; mask &= mask - 1) count++;
  return count;
}
function masksFor(values, columns) {
  return columns.reduce((mask, column, i) => mask | ((weight(values & column) % 2) << i), 0);
}

const golay = [0];
const coefficients = [7, 1, 4, 4, 2, 7, 3, 3, 4, 7, 2, 1];
for (let i = 0; i < 12; i++) {
  const v = coefficients.reduce((sum, c, j) => sum + ((c & 1) << (i + j)), 1 << 23);
  golay.push(...golay.map(word => word ^ v));
}
assert.equal(new Set(golay).size, 4096);
assert.equal(patterns[0].ids.length, 629);
assert.equal(patterns[0].pivot, 3);
assert.equal(patterns[0].mask, 3);
const old = patterns[0].ids.map(id => (golay[id % 4096] | golay[id >> 12]) & ~8);
for (const [a, b] of [...weight10, ...extra]) old.push(a | b);
assert.equal(old.length, 993);
assert.equal(old.filter(m => weight(m) === 10).length, 11);
assert.equal(old.filter(m => weight(m) === 11).length, 982);

const timesOmega = x => (x << 1) ^ ((x & 2) ? 7 : 0);
const axes6 = [];
for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) for (let c = 0; c < 4; c++) {
  axes6.push([a, b, c, a ^ b ^ c, a ^ timesOmega(b) ^ timesOmega(timesOmega(c)),
    a ^ timesOmega(timesOmega(b)) ^ timesOmega(c)]);
}
const axes7 = [];
const hamming = [1, 2, 3, 4, 5, 6, 7];
const permuted = [7, 3, 4, 6, 2, 1, 5];
const phase = [15, 0, 8, 1];
for (let x = 0; x < 16; x++) {
  const a = masksFor(x >> 1, hamming) ^ ((x & 1) ? 127 : 0);
  const b = phase.reduce((v, g, i) => v ^ ((x >> i & 1) ? g : 0), 0);
  for (let y = 0; y < 8; y++) {
    const c = b ^ masksFor(y, permuted);
    axes7.push(hamming.map((_, i) => (a >> i & 1) + 2 * (c >> i & 1)));
  }
}
for (const [axes, cardinality] of [[axes6, 64], [axes7, 128]]) {
  assert.equal(axes.length, cardinality);
  let min = axes[0].length;
  for (let i = 0; i < axes.length; i++) for (const v of axes.slice(0, i)) {
    min = Math.min(min, v.filter((c, j) => c !== axes[i][j]).length);
  }
  assert.equal(min, 4);
}

const dot = (a, b) => a.reduce((sum, x, i) => sum + x * b[i], 0);
const row = (index, value) => [0, 1, 2, 3].map(i => i === index ? value : 0);
const parityRows = [[], []];
for (let m = 0; m < 16; m++) parityRows[weight(m) % 2].push([0, 1, 2, 3].map(i => (m >> i & 1) ? -1 : 1));
for (const family of parityRows) {
  assert.equal(family.length, 8);
  for (let i = 0; i < family.length; i++) for (const v of family.slice(0, i)) {
    assert.ok([0, -4].includes(dot(family[i], v)));
  }
}
for (const a of parityRows[0]) for (const b of parityRows[1]) assert.equal(Math.abs(dot(a, b)), 2);
for (const p of parityRows.flat()) for (let i = 0; i < 4; i++) for (const x of [-2, 2]) {
  assert.ok(Math.abs(dot(p, row(i, x))) <= 2);
  for (let j = i + 1; j < 4; j++) for (const y of [-2, 2]) {
    assert.ok(Math.abs(dot(p, row(i, x).map((a, k) => a + row(j, y)[k]))) <= 4);
  }
}

const stocks = [
  ...six.map(([m, p]) => [m, p, 6]),
  ...seven.map(([m, p]) => [m, p, 7]),
];
assert.equal(six.length, 114);
assert.equal(seven.length, 131);
assert.equal(new Set(stocks.map(([m]) => m)).size, stocks.length);
for (let i = 0; i < stocks.length; i++) {
  const [m, p, w] = stocks[i];
  assert.ok(m >= 0 && m < 1 << 22);
  assert.equal(weight(m), w);
  assert.equal(weight(p), Math.floor(w / 2));
  assert.equal(m & p, p);
  // Same mode: distance-four axes or same axes with even sign distance at least two.
  assert.ok(4 * (w - 4) <= 2 * w);
  assert.ok(4 * w - 16 <= 2 * w);
  // Complementary modes achieve at most two per occupied row.
  const denseCap = 2 * w + 2 * Math.ceil(w / 2);
  assert.ok(4 * denseCap ** 2 <= 4 * w * 88);
  for (const o of old) {
    const sparseCap = 4 * weight((m << 2) & o);
    assert.ok(4 * sparseCap ** 2 <= 4 * w * 8 * weight(o));
  }
  for (let j = 0; j < i; j++) {
    const [n, q, k] = stocks[j];
    const shared = m & n;
    const overlap = weight(shared);
    const mismatches = weight(shared & (p ^ q));
    const cap = 4 * overlap - 2 * Math.min(mismatches, overlap - mismatches);
    assert.ok(4 * cap ** 2 <= 16 * w * k);
  }
}

const flat = six.length * 2 * axes6.length * 2 ** 5 + seven.length * 2 * axes7.length * 2 ** 6;
const sparse = 145004544;
const gain = flat - 79 * 4096;
for (const [dimension, fibers] of [[91, 24], [92, 32], [93, 44]]) {
  console.log(JSON.stringify({dimension, lowerBound: fibers * 8 ** 10 + sparse + flat, gain, flat}));
}

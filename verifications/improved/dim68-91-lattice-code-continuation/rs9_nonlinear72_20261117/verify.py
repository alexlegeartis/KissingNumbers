import itertools,json
from pathlib import Path
P=Path(__file__).resolve().parent
S=[0,7,25,42,63,102,105,115,124,163,173,182,184,203,204,210,213,224,269,274,305,323,328,340,351,357,378,385,394,407,412,420,443,454,473,495,542,545,556,562,578,589,599,600,646,648,659,669,705,743,746,756,761,772,795,807,824,846,849,864,875,886,893,911,912,930,937,949,958,986,1004,1011]
assert len(S)==len(set(S))==72 and all(0<=x<1024 for x in S)
assert min((a^b).bit_count()for a,b in itertools.combinations(S,2))==3
rows=[]
for n,M in [(87,9),(88,12),(89,16),(90,27),(91,27)]:
 k=min(9,n-81);supports=k*9**5;sparse=72*supports+2*n*(n-1);dense=2*M*3**18
 rows.append(dict(dimension=n,supports=supports,signs_per_support=72,dense=dense,sparse=sparse,total=dense+sparse))
json.dump({'source':'https://users.aalto.fi/~pat/opt/72.txt','source_description':'First code (first72 words after file header562,10,72), authors electronic supplement to Östergård–Baicheva–Kolev1999.','sign_masks':S,'minimum_hamming_distance':3,'antipodal':False,'rows':rows},open(P/'certificate.json','w'),indent=2)
print(json.dumps(rows,indent=2))

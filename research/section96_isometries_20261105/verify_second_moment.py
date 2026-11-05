import numpy as np
G=np.zeros((8,20),np.int32);g=np.array([7,1,4,4,2,7,3,3,4,7,2,1])
full=np.zeros((12,24),np.int32)
for i in range(12):full[i,i:i+12]=g;full[i,23]=3
G[:]=full[4:,4:];hist=np.zeros(21,np.int64);excluded=0
for b in range(256):
 t=np.arange(b*65536,(b+1)*65536,dtype=np.uint32);a=np.stack([(t>>(3*i))&7 for i in range(8)],axis=1).astype(np.int32);word=(a@G)%8;ok=np.all(word!=4,axis=1);excluded+=int((~ok).sum());zeros=(word[ok]==0).sum(axis=1);hist+=np.bincount(zeros,minlength=21)
 if b%32==31:print('batch',b+1,flush=True)
print('hist',hist.tolist(),'excluded',excluded,flush=True)
from fractions import Fraction
expected=[45056,245760,130560,430080,130720,107520,59424,3840,9032,0,0,0,0,0,0,0,0,0,0,0,1]
assert hist.tolist()==expected and excluded==15615223
weighted=sum(int(count)*2**(20+k) for k,count in enumerate(hist));assert weighted==18546205655040
moment=Fraction(21743271936*weighted,6**4*8**20);assert moment==Fraction(34545,128)
covered=Fraction(16**2,1)/moment;assert covered==Fraction(32768,34545)
print('PASS second moment',moment,'coverage lower bound',covered,'=',float(covered))

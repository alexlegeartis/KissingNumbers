"""Split paired tier candidates and tier seed deletions by hemisphere."""
import numpy as np
R='research/mixed_tier17_20261104/';F=np.fromfile(R+'mixed17_equator_frontier.bin',np.int32).reshape(-1,6);n=len(F);oldw=np.fromfile(R+'mixed17_switch_weights.bin',np.int32);E=np.fromfile(R+'mixed17_switch_edges.bin',np.int32).reshape(-1,2);base=2*n;new=[]
for a,b in E:
 a=int(a);b=int(b)
 if b<n:new.append((a,b));new.append((a+n,b+n));continue
 j=b-n
 if j<480:new.append((a,base+j));new.append((a+n,base+j))
 elif j<512:new.append((a,base+j));new.append((a+n,base+j+32))
 elif j<524:new.append((a,base+j+32))
 else:new.append((a+n,base+j+32))
weights=np.r_[oldw[:n]//2,oldw[:n]//2,oldw[n:n+480],np.full(64,16),np.full(24,16)].astype(np.int32);assert len(weights)==base+568
np.asarray(new,np.int32).tofile(R+'mixed17_asymmetric_edges.bin');weights.tofile(R+'mixed17_asymmetric_weights.bin');print('per hemisphere',n,'old atoms',568,'edges',len(new))

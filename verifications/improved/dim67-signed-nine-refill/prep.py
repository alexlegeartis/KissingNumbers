from pathlib import Path
import hashlib,json,numpy as np,sys,importlib.util,struct
D=Path(__file__).resolve().parent.parent
if not (D/'dim61-68-95-96-code-layers').is_dir(): D=D/'mid_49_69_20261121'/'verifications'/'improved'
S=D/'dim61-68-95-96-code-layers/source'
P=D/'dim61-67-compatible-sign-families/data/subcode67.json'
def blob(p):
 q=p.read_bytes();return hashlib.sha1(f'blob {len(q)}\0'.encode()+q).hexdigest()
def input67(expanded=None):
 assert hashlib.sha256((S/'dim61_96/shortening_support_search.json').read_bytes()).hexdigest()=='5422ed8f2b3be75267a57ea6fccc88f403a80e5bbe60dc64adc8edda51b61d88'
 assert hashlib.sha256(P.read_bytes()).hexdigest()=='90a421e678dd2f61f9fd668c374152145a6937dcfa4cbd4f5107825a44b5fe93'
 ext=Path(expanded) if expanded else S/'shared_parity/expanded67_68.npz'
 if not ext.is_file():ext=Path(__file__).with_name('expanded67_68.npz')
 assert blob(ext)=='92b815f8cea3197578fae34336095f418abd754a'
 O=np.load(ext,allow_pickle=False)
 assert len(O['new_supports'])==1792
 p=Path(S/'dim61_96/cyclic73.weight16.bin')
 assert hashlib.sha256(p.read_bytes()).hexdigest()=='5d240938e26a7a477de62bc255439e743cb33120bb13b986a70c6b1fea0b5239'
 R=[int(x)|(int(y)<<64)for x,y in np.fromfile(p,'<u8').reshape(-1,2)]
 U=json.loads((S/'dim61_96/shortening_support_search.json').read_text())['dimensions']['67']
 mark=(21,37,38,49,57,63,67,70,72);rem=sum(1<<i for i in U['deleted_coordinates'])
 coords=[i for i in range(1,73) if i not in mark]+[0]+[i for i in mark if not rem>>i&1]
 assert len(set(coords))==67
 V=[sum(1<<j for j,i in enumerate(coords)if x>>i&1)for x in R if not x&rem]
 assert len(V)==40428
 J=json.loads(P.read_text())
 assert J['total']==355095446
 W=V+list(map(int,O['new_supports']))+[int(o['support'])for o in J['new_classes']]
 assert len(W)==len(set(W)) and all(k.bit_count()==16 for k in W)
 p=S/'continuation49_72/netto.py'
 assert hashlib.sha256(p.read_bytes()).hexdigest()=='c3e73dd365bc432d625bfeb106096d19cdf090772b9717b5f2b4a097496bb855'
 q=importlib.util.spec_from_file_location('netto_source',p);net=importlib.util.module_from_spec(q);q.loader.exec_module(net)
 T=[sum(1<<i for i in b)for b in net.blocks68() if max(b)<67]
 assert len(T)==11792
 return W,T
def save(P,X):
 with open(P,'wb')as f:
  for x in X:f.write(struct.pack('<QQ',x&((1<<64)-1),x>>64))
if __name__=='__main__':
 W,T=input67();save('stock67.bin',W);save('netto67.bin',T);print('stock',len(W),len(T))

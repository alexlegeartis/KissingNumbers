from pathlib import Path
import hashlib,json,numpy as np,sys,importlib.util,struct
D=Path(__file__).resolve().parent.parent
if not (D/'dim61-68-95-96-code-layers').is_dir(): D=D/'mid_49_69_20261121'/'verifications'/'improved'
S=D/'dim61-68-95-96-code-layers/source'
P=D/'dim61-62-66-affine-reconstruction/dim66/certificate.json'
A=D/'dim61-66-nonary-affine-caps/certificate.json'
def blob(p):
 q=p.read_bytes();return hashlib.sha1(f'blob {len(q)}\0'.encode()+q).hexdigest()
def input66():
 assert hashlib.sha256((S/'dim61_96/shortening_support_search.json').read_bytes()).hexdigest()=='5422ed8f2b3be75267a57ea6fccc88f403a80e5bbe60dc64adc8edda51b61d88'
 assert hashlib.sha256(P.read_bytes()).hexdigest()=='26a393f7e6202f84fb5a945b95bc623b4e1dde8711d8ce384f7bc9ed396fc501'
 assert hashlib.sha256(A.read_bytes()).hexdigest()=='d2f7743669ca5ef39a6f7a58fdffc859bbc84550c815fcd89a6d958cb8b3f6de'
 p=Path(S/'dim61_96/cyclic73.weight16.bin')
 assert hashlib.sha256(p.read_bytes()).hexdigest()=='5d240938e26a7a477de62bc255439e743cb33120bb13b986a70c6b1fea0b5239'
 R=[int(x)|(int(y)<<64)for x,y in np.fromfile(p,'<u8').reshape(-1,2)]
 U=json.loads((S/'dim61_96/shortening_support_search.json').read_text())['dimensions']['66']
 mark=(21,37,38,49,57,63,67,70,72);rem=sum(1<<i for i in U['deleted_coordinates'])
 coords=[i for i in range(1,73) if i not in mark]+[0]+[i for i in mark if not rem>>i&1]
 assert len(set(coords))==66
 V=[sum(1<<j for j,i in enumerate(coords)if x>>i&1)for x in R if not x&rem]
 assert len(V)==34028
 J=json.loads(P.read_text())
 W=V+[int(o['support'])for o in J['new_classes']+J['family']['new']]
 assert len(W)==len(set(W)) and all(k.bit_count()==16 for k in W)
 p=S/'continuation49_72/netto.py'
 assert hashlib.sha256(p.read_bytes()).hexdigest()=='c3e73dd365bc432d625bfeb106096d19cdf090772b9717b5f2b4a097496bb855'
 q=importlib.util.spec_from_file_location('netto_source',p);net=importlib.util.module_from_spec(q);q.loader.exec_module(net)
 T=[sum(1<<i for i in b)for b in net.blocks68() if max(b)<66]
 assert len(T)==11088
 U=json.loads(A.read_text())['supports']['66']; assert len(U)==138
 return W,T,U
def save(P,X):
 with open(P,'wb')as f:
  for x in X:f.write(struct.pack('<QQ',x&((1<<64)-1),x>>64))
if __name__=='__main__':
 W,T,U=input66();save('stock66.bin',W);save('netto66.bin',T);save('old9_66.bin',U);print('stock',len(W),len(T))


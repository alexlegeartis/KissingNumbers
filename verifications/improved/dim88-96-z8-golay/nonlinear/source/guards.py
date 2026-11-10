"""Production structural guards shared with semantic corruption tests."""
import itertools

def validate_source(data,source,generated_nr,parts,full):
 G=data['octacode_generator'];assert len(G)==4 and all(len(row)==8 and all(type(x)==int and 0<=x<4 for x in row) for row in G)
 gray=[0,1,3,2]
 explicit={sum(gray[sum(c[i]*G[i][j] for i in range(4))%4]<<(2*j) for j in range(8)) for c in itertools.product(range(4),repeat=4)}
 assert explicit==set(generated_nr)==set(data['NR_words']) and len(explicit)==len(data['NR_words'])==256
 reps=data['translate_representatives'];assert len(reps)==8 and all(type(x)==int and 0<=x<65536 for x in reps)
 expected=[{w^t for w in explicit} for t in reps]
 assert len(data['partition'])==8 and all(expected[i]==set(parts[i])==set(data['partition'][i]) and len(data['partition'][i])==256 for i in range(8))
 H=set.union(*expected);assert len(H)==2048
 span={0}
 for x in data['span_basis']:
  assert type(x)==int and 0<=x<65536
  span|={y^x for y in span}
 assert H==span and len(span)==2048
 assert len(source)==len(set(source))==4096 and set(source)==set(full)
 assert all(type(x)==int and 0<=x<(1<<24) for x in source)

def validate_pattern(key,pat,equal):
 sizes=list(map(int,key.split(',')));weights=[i for i,n in enumerate(sizes) for _ in range(n)] if equal else [2]*sizes[0]+[4]*sizes[1]
 assert pat['family'] in ('linear','NR')
 assert type(pat.get('translation',0))==int and 0<=pat.get('translation',0)<1<<24
 if pat['family']=='linear':assert pat.get('translation',0)==0
 if equal and len({w%2 for w in weights})>1:
  assert pat['count']==0
  return None
 groups=pat['positions'];flat=sum(groups,[])
 assert list(map(len,groups))==weights and len(flat)==len(set(flat)) and all(type(x)==int and 0<=x<24 for x in flat)
 return groups

def count_pattern(key,pat,equal,linear,nr):
 groups=validate_pattern(key,pat,equal)
 if groups is None:return 0
 base=linear if pat['family']=='linear' else nr;translation=pat.get('translation',0);masks=[sum(1<<i for i in g) for g in groups]
 if equal:n=sum(len({len(g)-2*((w^translation)&m).bit_count() for g,m in zip(groups,masks)})==1 for w in base)
 else:n=sum(all(2*((w^translation)&m).bit_count()==len(g) for g,m in zip(groups,masks)) for w in base)
 assert n==pat['count']
 return n

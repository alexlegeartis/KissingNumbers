def need(ok,msg):
 if not ok:raise ValueError(msg)
def span(rows):
 a=[0]
 for h in rows:a += [x^h for x in a]
 return a
def expand(s,w):
 return sum(((w>>j)&1)<<i for j,i in enumerate(i for i in range(s.bit_length()) if s>>i&1))

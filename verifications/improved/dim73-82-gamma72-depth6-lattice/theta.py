from sympy import Matrix,divisor_sigma
from fractions import Fraction as F
from math import comb
N=32
def mul(a,b):
 v=[0]*N
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b[:N-i]):
    if y:v[i+j]+=x*y
 return v
one=[1]+[0]*(N-1)
def power(a,k):
 R=one
 for _ in range(k):R=mul(R,a)
 return R
A=[1]+[24*(divisor_sigma(i,1)-2*(divisor_sigma(i//2,1)if i%2==0 else 0)) for i in range(1,N)]
B=[-3]+[240*(divisor_sigma(i,3)-4*(divisor_sigma(i//2,3)if i%2==0 else 0)) for i in range(1,N)]
def eta(a,b):
 R=one
 for i in range(1,N):
  v=one[:];v[i]=-1
  ex=a+(b if i%2==0 else 0)
  if ex>=0:R=mul(R, power(v,ex))
  else:
   w=one[:]
   for j in range(1,N//i+1):
    if i*j<N:w[i*j]=comb(-ex+j-1,j)
   R=mul(R,w)
 return R

# Untwisted Γ72 full original shell:
E4=[1]+[240*divisor_sigma(i,3)for i in range(1,N)]
Del=[0]+eta(24,0)[:-1]
P=power(E4,9)
for j in (1,2,3):
 U=mul(power(E4,9-3*j),power(Del,j))
 assert U[j]==1
 P=[x-P[j]*y for x,y in zip(P,U)]
assert P[:5]==[1,0,0,0,6218175600]
# Γ0(2) quotient is not needed for isotropic shifts. Use both cusp jets.
X=[];Y=[]
for b in range(10):
 a=18-2*b
 R=mul(power(A,a),power(B,b))
 X.append(R);Y.append([(-1)**(a+b)*z for z in R])
C=Matrix(X).T;D=Matrix(Y).T
assert C[:10,:].det()!=0
K=C[:4,:].col_join(D[:6,:])
v=K.inv()*Matrix([1]+[0]*9)
S=D*v/2**18
assert S[:7]==[0,0,0,0,0,0,136704]
print('PASS isotropic depth 6 count',S[6])

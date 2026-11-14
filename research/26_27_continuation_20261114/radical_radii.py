"""Integer sign routines for the bounded discovery catalogues."""
import numpy as np

def qsign(a,b,r):
    a=np.asarray(a,np.int64);b=np.asarray(b,np.int64)
    return np.where(b==0,np.sign(a),np.where(b>0,np.where(a>=0,1,np.sign(r*b*b-a*a)),np.where(a<=0,-1,np.sign(a*a-r*b*b))))

def bsign(a,b,c,d,r,s):
    a=np.asarray(a,np.int64);b=np.asarray(b,np.int64);c=np.asarray(c,np.int64);d=np.asarray(d,np.int64)
    sx=qsign(a,b,r);sy=qsign(c,d,r);u=a*a+r*b*b-s*c*c-r*s*d*d;v=2*a*b-2*s*c*d;diff=qsign(u,v,r)
    return np.where(sy==0,sx,np.where(sx==0,sy,np.where(sx==sy,sx,sx*diff)))

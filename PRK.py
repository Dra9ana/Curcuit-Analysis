import sympy
from sympy import *
from sympy import symbols, ones, Eq, S
from sympy.integrals.transforms import laplace_transform
from sympy import I, Abs, conjugate
from sympy.plotting import plot

V2,V3,V4,U,R1,R2,R3,C1,C2,t,s,x=symbols("V2,V3,V4,U,R1,R2,R3,C1,C2,t,s,x");
e1=(U-V2)/R1-(V2-V3)/R3-V2/R2-C2*s*(V2-V4)
e2=(V2-V3)/R3-C1*s*V3
e3=x-V4/U
solution=nonlinsolve([e1,e2,e3],(V2,x,U)).subs([(V3,V4),(R1,2),(R2,2),(R3,1),(C2,sqrt(2)),(C1,1/sqrt(2))])
H=simplify(next(iter(solution))[1])
solveset(H,s);
w=Symbol('w',real=True)
A0=abs(H.subs(s,w*sqrt(-1)))
A=simplify(A0)
plot(A,(w,0,10))
Aref=simplify(limit(A,w,0))
w3dB=solveset(Eq(A,Aref/sqrt(2)),w,domain=Interval(0,sympy.oo))
f=inverse_laplace_transform(H/s,s,t)
print(f)
plot(f,(t,-20,20))
g=inverse_laplace_transform(H,s,t)
print(g)
plot(g)

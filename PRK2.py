import sympy
from sympy import *
from sympy import symbols, ones, Eq, S
from sympy.integrals.transforms import laplace_transform
from sympy import I, Abs, conjugate
from sympy.plotting import plot

R1,t,s,Ug,Tau, zc,zt,I1,I2,U1,U2=symbols("R1,t,s,Ug,tau, zc,zt,I1,I2,U1,U2");

#Ug=laplace_transform(Heaviside(t),t,s)[0]
e1=zc*I1+zc*I2*exp(-s*Tau)+U2*exp(-s*Tau)-U1;
e2=zc*I2+zc*I1*exp(-s*Tau)+U1*exp(-s*Tau)-U2;
e3=U1+R1*I1-Ug;


solution=linsolve([e1,e2,e3],(U1,U2,I1,I2)).subs([(R1,zc),(zc,50),(Tau,0.001),(I2,0)])
Ut=simplify(next(iter(solution))[1])
U11=simplify(next(iter(solution))[0])
solution1=nonlinsolve([e1,e2,e3],(U1,I1,U2,I2)).subs([(R1,zc),(zc,50),(t,0.001),(Ug,0)])
I22=simplify(next(iter(solution1))[3])
zt=U2/I22

Ut=Ut.subs(Ug,laplace_transform(Heaviside(t),t,s)[0])
pprint(Ut)
plot(inverse_laplace_transform(Ut,s,t),(t,-0.01,0.01))
U11=U11.subs(Ug,laplace_transform(Heaviside(t),t,s)[0])
plot(inverse_laplace_transform(U11,s,t),(t,-0.01,0.01))

clc
clear all
syms ms mp
syms t theta(t) a(t) b(t) phi(t)
syms l R g
syms u1 u2
tdot=diff(theta,t)
tdd=diff(tdot,t)
adot=diff(a,t)
add=diff(adot,t)
bdot=diff(b,t)
bdd=diff(bdot,t)
pdot=diff(phi,t)
pdd=diff(pdot,t)
Lx=0.5*ms*(-R*tdot)^2+0.5*mp*((-R*tdot+(adot-tdot)*l*cos(a-theta))^2+((adot-tdot)*l*sin(a-theta))^2)+0.5*(2/3*ms*R^2)*(-tdot)^2+1/2*(mp*l^2/12+mp*(l/2)^2)*(adot-tdot)^2-mp/l*((-R*tdot+(adot-tdot)*l*cos(a-theta))^2+((adot-tdot)*l*sin(a-theta))^2)-mp*g*l*cos(a-theta)
Ly=0.5*ms*(R*pdot)^2 + 0.5*mp*((R*pdot+(bdot-pdot)*l*cos(b-phi))^2 +((bdot-pdot)*l*sin(b-phi))^2)+0.5*(2/3*ms*R^2)*(-pdot)^2+0.5*((mp*l^2)/12+mp*(l/2)^2)*(bdot-pdot)^2-(mp*(R*pdot+(bdot-pdot)*l*cos(b-phi))^2 +((bdot-pdot)*l*sin(b-phi))^2)/l-mp*g*l*cos(b-phi) 
%% Lx, q1
%M11, M12 and V11 from here

%c=dLx/d(thetadot)
c=(mp*(2*(R + l*cos(a - theta))*(R*tdot - l*cos(a - theta)*(adot - tdot)) - l^2*sin(a - theta)^2*(2*adot - 2*tdot)))/2 - (mp*(2*(R + l*cos(a - theta))*(R*tdot - l*cos(a - theta)*(adot - tdot)) - l^2*sin(a - theta)^2*(2*adot - 2*tdot)))/l - (l^2*mp*(2*adot - 2*tdot))/6 + (5*R^2*ms*tdot)/3
%d=dLx/d(theta)
d=(mp*(2*l^2*sin(a - theta)*cos(a - theta)*(adot - tdot)^2 + 2*l*sin(a - theta)*(adot - tdot)*(R*tdot - l*cos(a - theta)*(adot - tdot))))/l - (mp*(2*l^2*sin(a - theta)*cos(a - theta)*(adot - tdot)^2 + 2*l*sin(a - theta)*(adot - tdot)*(R*tdot - l*cos(a - theta)*(adot - tdot))))/2 - g*l*mp*sin(a - theta)
%Euler Lagrange
simplify(diff(c,t)-d)
%% Lx, q2
%M21, M22 and V22 from here

%e=dLx/d(adot)
e=(mp*(l^2*sin(a - theta)^2*(2*adot - 2*tdot) - 2*l*cos(a - theta)*(R*tdot - l*cos(a - theta)*(adot - tdot))))/2 - (mp*(l^2*sin(a - theta)^2*(2*adot - 2*tdot) - 2*l*cos(a - theta)*(R*tdot - l*cos(a - theta)*(adot - tdot))))/l + (l^2*mp*(2*adot - 2*tdot))/6
%f=dLx/d(a)
f=(mp*(2*l^2*sin(a - theta)*cos(a - theta)*(adot - tdot)^2 + 2*l*sin(a - theta)*(adot - tdot)*(R*tdot - l*cos(a - theta)*(adot - tdot))))/2 - (mp*(2*l^2*sin(a - theta)*cos(a - theta)*(adot - tdot)^2 + 2*l*sin(a - theta)*(adot - tdot)*(R*tdot - l*cos(a - theta)*(adot - tdot))))/l + g*l*mp*sin(a - theta)
%Euler Lagrange
simplify(diff(e,t)-f)

%% Ly, q3
%M33, M34 and V3 from here

%g=dLy/d(bdot)
g=(mp*(l^2*sin(b - phi)^2*(2*bdot - 2*pdot) + 2*l*cos(b - phi)*(R*pdot + l*cos(b - phi)*(bdot - pdot))))/2 - (l^2*sin(b - phi)^2*(2*bdot - 2*pdot) + 2*l*mp*cos(b - phi)*(R*pdot + l*cos(b - phi)*(bdot - pdot)))/l + (l^2*mp*(2*bdot - 2*pdot))/6
%h=dLy/d(b)
h=(mp*(2*l^2*sin(b - phi)*cos(b - phi)*(bdot - pdot)^2 - 2*l*sin(b - phi)*(bdot - pdot)*(R*pdot + l*cos(b - phi)*(bdot - pdot))))/2 - (2*l^2*sin(b - phi)*cos(b - phi)*(bdot - pdot)^2 - 2*l*mp*sin(b - phi)*(bdot - pdot)*(R*pdot + l*cos(b - phi)*(bdot - pdot)))/l + g*l*mp*sin(b - phi)
%Euler Lagrange
simplify(diff(g,t)-h)
%% Ly, q4
%M43, M44 and V44 from here

%i=dLy/d(pdot)
i=(mp*(2*(R - l*cos(b - phi))*(R*pdot + l*cos(b - phi)*(bdot - pdot)) - l^2*sin(b - phi)^2*(2*bdot - 2*pdot)))/2 - (2*mp*(R - l*cos(b - phi))*(R*pdot + l*cos(b - phi)*(bdot - pdot)) - l^2*sin(b - phi)^2*(2*bdot - 2*pdot))/l - (l^2*mp*(2*bdot - 2*pdot))/6 + (5*R^2*ms*pdot)/3
%j=dLy/d(p)
j=(2*l^2*sin(b - phi)*cos(b - phi)*(bdot - pdot)^2 - 2*l*mp*sin(b - phi)*(bdot - pdot)*(R*pdot + l*cos(b - phi)*(bdot - pdot)))/l - (mp*(2*l^2*sin(b - phi)*cos(b - phi)*(bdot - pdot)^2 - 2*l*sin(b - phi)*(bdot - pdot)*(R*pdot + l*cos(b - phi)*(bdot - pdot))))/2 - g*l*mp*sin(b - phi)
%Euler Lagrange
simplify(diff(i,t)-j)
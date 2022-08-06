clear;
% 
%Dutch and Spiral Approximation
%
b11=a11-a12*a21/a22
b12=a13-a12*a23/a22 
b13=a14; b21=a31-a32*a21/a22
b22=a33-a32*a23/a22
b23=0
b31=-a21/a22
b32=-a23/a22
b33=0
%
B=[b11 b12 b13;b21 b22 b23; b31 b32 b33] 

[u2 v2] = eig(B)
%
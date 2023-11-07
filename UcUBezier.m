%Hello! this will plot Bezier curve for n control points 
%This is a replacement of the program 'Parametic Cubic Bezier Curve'
%submitted before ...Bezier for any number of points ...enjoy 
function [y,x] = UcUBezier
%n=4;
% n=input('Enter no. of points  ');
%w=input('Press 1 for entry through mouse or 2 for keyboard repectively-->');
w=2;
%n1=n-1;
if w==1
    axis([0 100 0 100])
    [p]=ginput(4);
end
if w==2
    %[p]=[0 1;0.4 1; 0.5 0;0.6 0; 1 0];%
    [p]=[0 1; 0.3 0.9; 0.4 0; 1 0];%
    %input('Enter co-ordinates of points within brackets ->[x1 y1;x2 y2;x3 y3;...;xn yn] ');
end
n=length(p);
n1=n-1;
sigma=zeros(1,n-1);
    i=[0:1:n1];
    sigma(i+1)=factorial(n1)./(factorial(i).*factorial(n1-i));
% for    i=0:1:n1
%     sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
% end
l=[];
UB=zeros(1,n);
for u=0:0.002:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end
P=l*p;
%line(P(:,1),P(:,2))
%line(p(:,1),p(:,2))
x=P(:,2);
y=P(:,1);
end
% books for reference on the subject of cad/cam author Roger Adams  ,
% Ibrahim Zeid 
%Prepared by Mechanical engg. student NIT Allahabad , India
% for any questions feel free to mail me at slnarasimhan89@gmail.com
function [P,V,T] = LSex3(p0,v0)
%LS Summary of this function goes here
%   Detailed explanation goes here

if nargin==0
 p0=[0;1;0];v0=[1;1;1];
end


k=2000;

v=PrTpM(v0,hg(p0));

[P,V,T]=AB5se(p0,v);

while 1
   [P,V,T]= AB5st(P,V,T);
    if length(P)== k
       P1=P(1:1,1:k);
       P2=P(2:2,1:k);
       P3=P(3:3,1:k);
       %plot3(P1,P2,P3, 'Color', [255, 153, 51]/255); hold on; grid on; 
       plot3(P1,P2,P3, 'Color', [255, 0, 0]/255); hold on; grid on; 
       xlabel('x');
       ylabel('y');
       zlabel('z');
       %xlim([-3 3]); ylim([-3 3]);
       return 
    end
end
end

function [out] = hf(x)
out=[x(1)^2+x(2)^2+x(3)^2;x(1)^2-x(3)];
end

function [out]= hg(x)
out=[[2*x(1);2*x(2);2*x(3)],[2*x(1);0;-1]];
end
function [out] = hH(x)
    out=[2*eye(3),[2;0;0],[0;0;0],[0;0;0]];
end

function [proj] = PrTpM(v,J)
%PROJTPM Summary of this function goes here
%   Detailed explanation goes here

[n,m]=size(J);
[Q,~]=qr(J);

Q=Q(1:n,1:m);
Q=Q*Q.'; %overrides Q with proj NpM as a matrix

proj=v-Q*v;
proj=proj/norm(proj);
end

function [P,V,T] = AB5se(p0,v0,h)
%RK5SETUP Summary of this function goes here
%   Detailed explanation goes here

if nargin==2
  h=0.01;
end

%Euler 
p1=p0+h*v0;
t0=Mfeq(p0,v0);
v1=(v0+h*t0);

%AdamBash2
temp=3/2*v1-1/2*v0;
p2=(p1+h*temp);
t1=Mfeq(p1,v1);
temp=3/2*t1-1/2*t0;
v2=v1+h*temp;

%AB3
temp=23/12*v2-4/3*v1+5/12*v0;
p3=(p2+h*temp);
t2=Mfeq(p2,v2);
temp=23/12*t2-4/3*t1+5/12*t0;
v3=v2+h*temp;
%AB4
temp=55/24*v3-59/24*v2+37/24*v1-3/8*v0;
p4=(p3+h*temp);
t3=Mfeq(p3,v3);
temp=55/24*t3-59/24*t2+37/24*t1-3/8*t0;
v4=v3+h*temp;

P=round([p0,p1,p2,p3,p4],10);
V=round([v0,v1,v2,v3,v4],10);
T=round([t0,t1,t2,t3],10);
end

function [P,V,T] = AB5st(P,V,T,h)
%RK5 Summary of this function goes here
%   Detailed explanation goes here

n=3;

if nargin==3
  h=0.01;
end

[i,j]=size(P);

t=Mfeq(P(1:n,j),V(1:n,j));

p_x=1901/720*V(1:n,j)-2774/720*V(1:n,j-1)+2616/720*V(1:n,j-2)-1274/720*V(1:n,j-3)+251/720*V(1:n,j-4);
v_x=1901/720*t-2774/720*T(1:n,j-1)+2616/720*T(1:n,j-2)-1274/720*T(1:n,j-3)+251/720*T(1:n,j-4);
p_new=round(P(1:n,j)+h*p_x,10);
v_new=round(V(1:n,j)+h*v_x,10);
v_new=round(v_new/norm(v_new),10);

P=[P,p_new];
T=[T,t];
V=[V,v_new];

end


function [w]=Mfeq(q,v)
    eJ=hg(q);
    eH=hH(q);
    
    [n,m]=size(eJ);
    [Q,R]=qr(eJ);
    
    R=R(1:m,1:m);
    R=R.'; %why the heck is it R transposed?
    
    H=zeros([m 1]);
    for i=1:m
        H(i)=v.'*eH(1:n,(i-1)*n+1:i*n)*v;
    end
    t=R\H;
    w=-Q(1:n,1:m)*t;
end

function [w] = Logarithmic2(p,q)

global h;
h=0.01;

clf('reset');
if nargin==0
   p=[-1;2;0.2];
   q=[0.5;0.5;2]; 
   
   
   fsurf(@(x,y) 1/(x^2+y^2) ,[-2 2]);
   hold on; grid on;
   plot3(q(1),q(2),q(3),'-s','color', [0.93 0.7 0.1]);
   xlim([-1.1 0.8]);
   ylim([0 2.1]);
   zlim([0.2 3]);
   xlabel('x');
   ylabel('y');
   zlabel('z');
   drawnow;
end

%setup
%f(x)=(x(1)^2+x(2)^2)*x(3)-1;

v=q-p;
v=PrTpM(v,hg(p));
%v=round(v/norm(v),10);

drawcounter=1;

while 1
    %outer loop -> correct v
        %setup for new cycle
        [n,~]=size(p);
        p0=zeros(n,n+1);
        v0=p0;
        p0(1:n,1)=p;
        p0(1:n,2:n+1)=zeros(n);
        v0(1:n,1)=PrTpM(v,hg(p));
        v0(1:n,2:n+1)=eye(3);
        [P,V,T]=AB5se(p0,v0);
        while length(P)<1/h
            %inner loop -> geodesic for specific v
            [P,V,T]= AB5st(P,V,T);
        end
        
        
        k=length(P);
        s=P{k}(1:n,1);
        
        round([s,v0(1:n,1)],10)
        
        P1=zeros(1,k);
        P2=zeros(1,k);
        P3=zeros(1,k);
        for i=1:k
            P1(i)=P{i}(1,1);
            P2(i)=P{i}(2,1);
            P3(i)=P{i}(3,1);
        end
        
        if drawcounter==1
            col=[1 1 1];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1'})
            drawcounter=drawcounter+1;
        elseif drawcounter==2
            col=[0 1 1];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2'})
            drawcounter=drawcounter+1;
        elseif drawcounter==3
            col=[1 0 1];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2','v3'})
            drawcounter=drawcounter+1;
        elseif drawcounter==4
            col=[1 1 0];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2','v3','v4'})
            drawcounter=drawcounter+1;
        elseif drawcounter==5
            col=[0.635 0.08 0.19];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2','v3','v4','v5'})
            drawcounter=drawcounter+1;
        elseif drawcounter==6
            col=[0.3 0.745 0.93];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2','v3','v4','v5','v6'})
            drawcounter=drawcounter+1;
        elseif drawcounter==7
            col=[1 0 0];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2','v3','v4','v5','v6','v7'})
            drawcounter=drawcounter+1;
        elseif drawcounter==8
            col=[0 1 0];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2','v3','v4','v5','v6','v7','v8'})
            drawcounter=drawcounter+1;
        else
            col=[0.5 0.2 0.56];
            plot3(P1,P2,P3,'Color', col);
            %plot3(P1(k),P2(k),P3(k),'-s', 'Color', col);
            legend({'M','Target','v1','v2','v3','v4','v5','v6','v7','v8','v9'})

        end
        
        drawnow;

        if norm(s-q)<0.005
            w=v;
            return
        else 
            cor=P{k}(1:n,2:n+1);
            w=s-q;
            cor=cor\w;
            v=v-cor;
        end
        
    end

end

function [out] =hg(x)
    out=[2*x(1)*x(3);2*x(2)*x(3);(x(1)^2+x(2)^2)];
end

function [out]=hH(x)
    out=[[2*x(3);0;2*x(1)],[0;2*x(3);2*x(2)],[2*x(1);2*x(2);0]];
end

function [out] = HH(x)
    temp1=zeros(3);
    temp2=zeros(3);
    temp3=zeros(3);
    temp1(3,1)=2;
    temp1(1,3)=2;
    temp2(2,3)=2;
    temp2(3,2)=2;
    temp3(1,1)=2;
    out={temp1,temp2,temp3};
    
end

function [P,V,T] = AB5se(p0,v0,h)
%AB5SETUP Setup for a 5 step ODE solver
%   Detailed explanation goes here

if nargin==2
  global h;
end

%Euler 
p1=p0+h*v0;
t0=Mflog(p0,v0);
v1=(v0+h*t0);

%AdamBash2
temp=3/2*v1-1/2*v0;
p2=round(p1+h*temp,10);
t1=Mflog(p1,v1);
temp=3/2*t1-1/2*t0;
v2=v1+h*temp;

%AB3
temp=23/12*v2-4/3*v1+5/12*v0;
p3=round(p2+h*temp,10);
t2=Mflog(p2,v2);
temp=23/12*t2-4/3*t1+5/12*t0;
v3=v2+h*temp;
%AB4
temp=55/24*v3-59/24*v2+37/24*v1-3/8*v0;
p4=round(p3+h*temp,10);
t3=Mflog(p3,v3);
temp=55/24*t3-59/24*t2+37/24*t1-3/8*t0;
v4=v3+h*temp;

Pr=round([p0,p1,p2,p3,p4],10);
Vr=round([v0,v1,v2,v3,v4],10);
Tr=round([t0,t1,t2,t3],10);

[n,~]=size(Pr);
m=n+1;
P={Pr(1:n,1:m),Pr(1:n,m+1:2*m),Pr(1:n,2*m+1:3*m),Pr(1:n,3*m+1:4*m),Pr(1:n,4*m+1:5*m)};
V={Vr(1:n,1:m),Vr(1:n,m+1:2*m),Vr(1:n,2*m+1:3*m),Vr(1:n,3*m+1:4*m),Vr(1:n,4*m+1:5*m)};
T={Tr(1:n,1:m),Tr(1:n,m+1:2*m),Tr(1:n,2*m+1:3*m),Tr(1:n,3*m+1:4*m)};
end

function [P,V,T] = AB5st(P,V,T,h)
%AB5 step
%   Detailed explanation goes here
if nargin==3
    global h;
end

[~,j]=size(P);

t=Mflog(P{j},V{j});

p_x=1901/720*V{j}-2774/720*V{j-1}+2616/720*V{j-2}-1274/720*V{j-3}+251/720*V{j-4};
v_x=1901/720*t-2774/720*T{j-1}+2616/720*T{j-2}-1274/720*T{j-3}+251/720*T{j-4};
p_new=round(P{j}+h*p_x,10);
v_new=round(V{j}+h*v_x,10);
%v_new=round(v_new/norm(v_new),10);

P{j+1}=p_new;
T{j}=t;
V{j+1}=v_new;

end


function [w]=Mflog(qin,vin)
    
    
    %get points and matrix
    [n,m]=size(qin);
    q=qin(1:n,1);
    v=vin(1:n,1);
    qm=qin(1:n,2:m);
    vm=vin(1:n,2:m);
    
    Hess=HH(q);
    
    %geodesic equation of point p to point q
    eJ=hg(q);
    eH=hH(q);
    
    [n,m]=size(eJ);
    [Q,R]=qr(eJ);
    
    R=R(1:m,1:m);
    R=R.'; 
    
    H=zeros([m 1]);
    for i=1:m
        H(i)=v.'*eH*v;
    end
    t=R\H;
    w=-Q(1:n,1:m)*t;
    
    %equation for pmatrix and vmatrix
    %we assume from here on that f=f_1 !

    no=norm(eJ)^2;
    
    H1=zeros([n 1]);
    for i=1:n
        H1(i)=v.'*Hess{i}*v;
    end
    
    tJ1=diag(H1);
    tJ2=diag(eH*eJ);
    
    J1=1/no*tJ1-2/no^2*tJ2*H;
    J2=zeros([n]);
    for i=1:n
        J2(1:n,i)=eJ;
    end
    
    J3=H/no*eH;
    Jgamma=-J1*J2-J3;
    
    K1=-eH*v*eJ.';
    
    Jdgamma=2/no*K1;
        
    temp=(Jgamma.'*qm+Jdgamma.'*vm);
    
    
    w=[w,temp];
end

function [proj] = PrTpM(v,J)
%PROJTPM Projection to the tangent space, with the normal space spanned by
%the coloumns of J
%   Detailed explanation goes here

[n,m]=size(J);
[Q,R]=qr(J);

Q=Q(1:n,1:m);
Q=Q*Q.'; %overrides Q with proj NpM as a matrix

proj=v-Q*v;
%proj=proj/norm(proj);
end
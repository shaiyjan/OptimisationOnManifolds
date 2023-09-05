function [out1] = Ladder1(p0,v,W)
%PARATRANS Summary of this function goes here
%   Detailed explanation goes here

global h;
h=0.01;
target=200;
geodsteps=50;

if nargin==0
   p0=[1;0;0];
   v=[0;1;0];
   W={[0;0;1]};
   clf;
   fsurf(@(x,y) sqrt(1-(x^2+y^2)), [-1 1], 'EdgeColor', 'none');
   hold on; grid on;
   xlabel('x');
   ylabel('y');
   zlabel('z');
   view(150,10);
   fsurf(@(x,y) -sqrt(1-(x^2+y^2)) ,[-1 1],'EdgeColor', 'none');
   %fsurf(@(x,y) 1.4 ,[-1 1]);
   zlim([-0.15 0.3]);
   
   th=0:pi/50:2*pi;
               
   for height=[-0.15,0.3,0.4]
       r=sqrt(1-height^2);
       xun=r*cos(th);
       yun=r*sin(th);
       zun=zeros(1,length(th))+height;
       plot3(xun,yun,zun,'color',0.5*[1,1,1]);
   end
   drawnow;
end

v=PrTpM(v,hg(p0));
v=round(v/norm(v),10);

[P,V,T]=AB5se(p0,v);
steps=5;
[~,m]=size(W);

Q=cell(1,m);
ParVec=W;
Norms=zeros(1,m);
for i=1:m
   temp=W{i};
   Norms(i)=norm(temp); 
end
ParPoi=cell(1,m);
tcount=1;
PT=0;

%a=P(k-steps); is p0
for i=1:m
    Q{i}=LS(p0,ParVec{i},1/(10*h));
end
while 1
    [P,V,T]= AB5st(P,V,T);
    
 if length(P)==target %random cancel  criteria
    PT=2;
 elseif steps==geodsteps
    PT=1;
 end
    
    if PT~=0
        k=length(P);
        
        %drawing
        P1=zeros(1,steps);
        P2=zeros(1,steps);
        P3=zeros(1,steps);
        for i=1:steps
            P1(i)=P{k-i}(1,1);
            P2(i)=P{k-i}(2,1);
            P3(i)=P{k-i}(3,1);
        end
        plot3(P1,P2,P3,'Color', 'k');
        %if (-1)^tcount==1
        %    plot3(P1,P2,P3,'Color', 'm');
        %else
        %    yel=[1 1 0];
        %    plot3(P1,P2,P3,'Color', 'y');
        %end
        
        drawnow
        %till here drawing

        
        
        mp=P(k-round(steps/2));
        for i=1:m
             Q{i}=Log(Q{i},mp{1},2);
             %hf(Q{i})
        end

        if PT==1
            PT=0;
            steps=0;
        elseif PT==2
            for i=1:m
                ParVec{i}=ParVec{i};
            end
            b=P(k);
            temp=-Log(b{1},Q{i},1);
             Q{i}=(-1)^tcount*temp*Norms(i)/norm(temp);
            out1=Q;
            if nargin==0
               cplot=fsurf(@(x,y) -sqrt(1-(x^2+y^2)) ,[0 0], 'EdgeColor', 'none');
               mplot=plot3(100,100,100,'Color','m');
               rplot=plot3(100,100,100,'Color','r');
               bplot=plot3(100,100,100,'Color','k');
               oplot=plot3(100,100,100,'Color',[0.9 0.6 0.2]);
               legend([cplot,bplot,mplot,rplot,oplot],'Sphere','Geodesic v','Geodesic w','Log center','Log transport');
               xlim([-1 1]);
               
            end
            return
        end
        
    tcount=tcount+1;    
    end
    steps=steps+1;
    
end

end

function [out] = hf(x)
    out=norm(x)-1;
end

function [out] = hg(x)
    out=2*x;
end 

function [out]= hH(x)
    out=2*eye(3);
end

function [out]=HH(x)
    out={zeros(3),zeros(3),zeros(3)};
end

function [out] = LS(p0,v,n)
%LS Calculate the geodesic with starting points p0 and initial velocity v
%in TpM, for up to n steps/ at t=n*h
[P,V,T]=AB5se(p0,v);
n=round(n);
while length(P)<n
    [P,V,T]= AB5st(P,V,T);
end
%drawing
   k=length(P);
   P1=zeros(1,k);
   P2=zeros(1,k);
   P3=zeros(1,k);
   for i=1:k
      P1(i)=P{i}(1,1);
      P2(i)=P{i}(2,1);
      P3(i)=P{i}(3,1);
   end
   plot3(P1,P2,P3,'Color', 'm');
   drawnow
%drawing

out=P{length(P)};
end

function [out] = Log(p,q,lcase)
%Log Solves the logarithmic problem, finding v in TpM such that the
%geodesic with starting point p and initial velocity v hits q at t=1
global  h;

v=q-p;
while 1
    [n,~]=size(p);
    p0=zeros(n,n+1);
    v0=p0;
    p0(1:n,1)=p;
    p0(1:n,2:n+1)=zeros(n);
    v0(1:n,1)=round(PrTpM(v,hg(p)),10);
    v0(1:n,2:n+1)=eye(3);
    [P,V,T]=AB5se(p0,v0);
    while length(P)<1/h
       %inner loop -> geodesic for specific v
       [P,V,T]= AB5st(P,V,T);
    end
        
    k=length(P);
    s=P{k}(1:n,1);
        
    if norm(s-q)<0.005
        break
    else 
        cor=P{k}(1:n,2:n+1);
        w=s-q;
        cor=cor\w;
        v=v-cor;
    end
end

    if lcase==1
        out=v;
        %drawing
        k=length(P);
        P1=zeros(1,k);
        P2=zeros(1,k);
        P3=zeros(1,k);
        for i=1:k
           P1(i)=P{i}(1,1);
           P2(i)=P{i}(2,1);
           P3(i)=P{i}(3,1);
        end
        orange=[0.9 0.6 0.2];
        plot3(P1,P2,P3,'Color', orange);
        drawnow
        %drawing
        return
    elseif lcase==2
        while length(P)<2/h
            %inner loop -> geodesic for specific v
            [P,V,T]= AB5st(P,V,T);
        end
        %drawing
        k=length(P);
        P1=zeros(1,k);
        P2=zeros(1,k);
        P3=zeros(1,k);
        for i=1:k
            P1(i)=P{i}(1,1);
            P2(i)=P{i}(2,1);
            P3(i)=P{i}(3,1);
        end
        red=[1 0 0];
        plot3(P1,P2,P3,'Color', red);
        drawnow
        %drawing end
        
        out=P{length(P)};
        out=out(1:n,1);
        return
    end
end

function [P,V,T] = AB5se(p0,v0,h)
%RK5Se  Setup to a 5 step ODE solver, with Euler and Adam-Bashford 2 to 4 step algorithm
%p0 starting point, v0 initial velocity, h steplength
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

[n,m]=size(p0);

%if m~=1
%    m=n+1;
%else
%    m=1;
%end

P={Pr(1:n,1:m),Pr(1:n,m+1:2*m),Pr(1:n,2*m+1:3*m),Pr(1:n,3*m+1:4*m),Pr(1:n,4*m+1:5*m)};
V={Vr(1:n,1:m),Vr(1:n,m+1:2*m),Vr(1:n,2*m+1:3*m),Vr(1:n,3*m+1:4*m),Vr(1:n,4*m+1:5*m)};
T={Tr(1:n,1:m),Tr(1:n,m+1:2*m),Tr(1:n,2*m+1:3*m),Tr(1:n,3*m+1:4*m)};
end

function [P,V,T] = AB5st(P,V,T,h)
%Adam-Bashford 5 step
%
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
%evaluation of geodesic ODE with option to log ODE
    
    [n,~]=size(qin);
    
    q=qin(1:n,1);
    v=vin(1:n,1);
    
    Hess=HH(q);
    
    %init points
    eJ=hg(q);
    eH=hH(q);
    
    [n,m]=size(eJ);
    
    
    
    %calculations
    [Q,R]=qr(eJ);
    
    R=R(1:m,1:m);
    R=R.'; 
    H=zeros([m 1]);
    for i=1:m
        H(i)=v.'*eH*v;
    end
    t=R\H;
    w=-Q(1:n,1:m)*t;
    
    [n,m]=size(qin);
    
    if m>1 
        % only entered for log - matrix valued ODE calc.
        
        %init matrix
        qm=qin(1:n,2:m);
        vm=vin(1:n,2:m);

        %calculations
        no=norm(eJ)^2;
    
        H1=zeros([n 1]);
        for i=1:n
            H1(i)=v.'*Hess{i}*v;
        end
    
        tJ1=diag(H1);
        tJ2=diag(eH*eJ);
    
        J1=1/no*tJ1-2/no^2*tJ2*H;
        J2=zeros(n);
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
    
end

function [proj] = PrTpM(v,J)
%PROJTPM Calc projection to TpM
%   v in R^n  J is transpose of jacobian of gradients

[n,m]=size(J);
[Q,~]=qr(J);

Q=Q(1:n,1:m);
Q=Q*Q.'; %overrides Q with proj NpM as a matrix

proj=v-Q*v;
%proj=proj/norm(proj);
end
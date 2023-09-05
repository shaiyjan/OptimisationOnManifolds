% This is the "perfect line search" version.


function tester()
global SET
cd('/home/gc1355/Downloads/Cutest')
addpath('/home/gc1355/Dropbox/Diss/Data/Definite_LBFGSMCutest_ALL_V3.m')
addpath('/home/gc1355/Downloads/Cutest/Gould_Ref/cutest/src/matlab')

MainFolder="/home/gc1355/Dropbox/Diss/Data/V4-"+string(datetime('today'));
mkdir(MainFolder);
for SET=["Circ","4Circ","CircInt","Stiefel"]
    destination=MainFolder+"/"+SET;
    mkdir(destination);
    for t={'ACOP','DEGE','NASH','AVGA','NET1'}
        file="/home/gc1355/Downloads/Cutest/"+t{1}+"/OUTSDIF.d";
        exc="/home/gc1355/Downloads/Cutest";
        copyfile(file,exc)
        folder=destination+"/"+t{1};
        mkdir(folder);
        try
            setup=cutest_setup();
        catch 
            cutest_terminate
            setup=cutest_setup();
        end
        fid=fopen(folder+"/String_"+SET+"_"+setup.name(1:4)+".txt","wt");
        fprintf(fid,"\\[\\begin{array}{|c|c|c|c|c|c|}\\hline \n");
        fprintf(fid,"\\multicolumn{3}{|c|}{\\text{Name:"+setup.name+"} } &\\multicolumn{3}{|c|}{\\text{Dim}: "+setup.n+"} \\\\ \\hline \n");
        fprintf(fid,"L & \\nabla L & f & It & Rest & LogCor\\\\ \\hline \n");
        for i=4:4
            for j=[0,setup.n]
                fprintf(SET+"-"+t{1}+"-"+i+"-"+j+"\n")
                [out,name]=LBFGSMCutest(10^(-i),j,setup);
                save(folder+"/"+name+".mat",'-struct','out');
                createString(out,fid);
            end
        end
        fprintf(fid,"\\end{array}\\]");
        cutest_terminate 
    end
end
end

function createString(data,fid)
    errorpow=floor(log10(data.error));
    error=num2str(data.error*10^(-errorpow),'%5.2f');
    fprintf(fid," \\multicolumn{1}{|c|}{D="+data.listlength+"} &\\multicolumn{2}{|c|}{h=10^{"+log10(data.stepsize)+"}} &  \\multicolumn{3}{|c|}{f_\\epsilon="+error+"\\cdot10^{"+errorpow+"}} \\\\ \\hline \n");
    
    ControlVal=inf;
    try
        if data.listlength==0
            StrRest="-";
            StrCor="-";
        else
            StrRest=data.restarts1;
            StrCor=data.retries1;
        end
        ControlVal=data.gradf1;
        GradPow=floor(log10(data.gradnorm1));
        GradDig=num2str(data.gradnorm1*10^(-GradPow),'%5.4f');
        Func=num2str(data.gradf1,'%5.4f');
        fprintf(fid, Func +" & "+ GradDig+"\\cdot10^{"+GradPow +"} &" + data.gradeval1 + " & " + data.gradit1 + " & " + StrRest + " & " + StrCor + "\\\\ \n");
    catch
    end
    try
        if data.listlength==0
            StrRest="-";
            StrCor="-";
        else
            StrRest=data.restarts2;
            StrCor=data.retries2;
        end
        ControlVal=data.gradf2;
        GradPow=floor(log10(data.gradnorm2));
        GradDig=num2str(data.gradnorm2*10^(-GradPow),'%5.4f');
        FuncDif=data.gradf1-data.gradf2;
        FuncPow=floor(log10(FuncDif));
        FuncDig=num2str(-FuncDif*10^(-FuncPow),'%5.4f');
        fprintf(fid,"\\ast"+ FuncDig+"\\cdot10^{"+FuncPow +"} & "+ GradDig+"\\cdot10^{"+GradPow +"} &" + data.gradeval2 + " & " + data.gradit2 + " & " + StrRest + " & " + StrCor + "\\\\ \n");
    catch
    end
    try
        if data.listlength==0
            StrRest="-";
            StrCor="-";
        else
            StrRest=data.restarts3;
            StrCor=data.retries3;
        end
        GradPow=floor(log10(data.gradnorm3));
        GradDig=num2str(data.gradnorm3*10^(-GradPow),'%5.4f');
        FuncDif=data.gradf2-data.gradf3;
        FuncPow=floor(log10(FuncDif));
        FuncDig=num2str(-FuncDif*10^(-FuncPow),'%5.4f');
        fprintf(fid,"\\ast"+ FuncDig+"\\cdot10^{"+FuncPow +"} & "+ GradDig+"\\cdot10^{"+GradPow +"} &" + data.gradeval3 + " & " + data.gradit3 + " & " + StrRest + " & " + StrCor + "\\\\ \n");
    catch
    end
    try
        if data.listlength==0
            StrRest="-";
            StrCor="-";
        else
            StrRest=data.restartsstep;
            StrCor=data.retriesstep;
        end
        GradPow=floor(log10(data.gradnormstep));
        GradDig=num2str(data.gradnormstep*10^(-GradPow),'%5.4f');
        if GradPow==0
            GradStr=""+GradDig;
        else
            GradStr=""+GradDig+"\\cdot10^{"+GradPow +"}";
        end
        if ControlVal==inf
        Func=num2str(data.gradfstep,'%5.4f%');
        fprintf(fid, Func +" & "+  GradStr+ "&" + data.gradevalstep + " & " + data.graditstep + " & " + StrRest + " & " + StrCor + "\\\\ \n");    
        else
        FuncDif=ControlVal-data.gradfstep;
        FuncPow=floor(log10(FuncDif));
        FuncDig=num2str(-FuncDif*10^(-FuncPow),'%5.4f');
        fprintf(fid,"\\ast"+ FuncDig+"\\cdot10^{"+FuncPow +"} & "+ GradStr+ "&" + data.gradevalstep + " & " + data.graditstep + " & " + StrRest + " & " + StrCor + "\\\\ \n");  
        end
    catch
    end
    fprintf(fid,"\\hline \n");
end

function[out,name] = LBFGSMCutest(s,leng,setup)
global SET row col maxstep listlength restart h minv dgf evaluations dim Gamma Delta rto m retries iterations;

h=s;
listlength=leng;
maxstep=15000;
gradval=10^(-2);

restart=1;
evaluations=1;
retries=0;

con=0;
rto=20;
iterations=0;
restarts=-1;
Gamma={};
Delta={};


out.algo="LBFGS3_"+SET;
out.name=setup.name;
out.dim=setup.n;
out.stepsize=h;
out.listlength=listlength;
dim=setup.n;
if SET=="Circ" || SET=="4Circ"
    m=1;
elseif SET=="CircInt"
    m=3;
elseif SET=="Stiefel"
    if dim==72
       row=9;
       col=8;
    elseif dim==8
        row=4;
        col=2;
    elseif dim==48
        row=8;
        col=6;
    elseif dim==20
        row=5;
        col=4;
    end
    m=col*(col+1)/2;
end


p=stpoint(setup.x);
minv=cutest_obj(p);
name=out.algo+"_"+setup.name(1:4)+"_"+string(h)+"_"+string(listlength);
while 1
   v=-round(PrTpM(cutest_grad(p),hg(p)),rto);
   dgf=norm(v);
   %Stopping criteria
   
   try
   if iterations==5000
         out.gradfstep=minv;
         out.gradnormstep=dgf;
         out.gradevalstep=evaluations;
         out.retriesstep=retries;
         out.restartsstep=restarts;
         out.graditstep=iterations;
         out.error=norm(hf(p));
         fprintf('Failed on steps \n' )
         %con=con+1;
         return
   elseif norm(hf(p))>= 1
         out.fail='Failed on distance /n';
         return
   elseif con==2
       if dgf<gradval*10^(-2) 
         out.gradf3=minv;
         out.crit3=gradval*10^(-2);
         out.gradnorm3=dgf;
         out.gradeval3=evaluations;
         out.retries3=retries;
         out.restarts3=restarts;
         out.gradit3=iterations;
         out.error=norm(hf(p));
         %fprintf('End by gradient \n' )
         %con=con+1;
         return
       end      
   elseif con==1
       if dgf<gradval*10^(-1)
         out.gradf2=minv;
         out.crit2=gradval*10^(-1);
         out.gradnorm2=dgf;
         out.gradeval2=evaluations;
         out.restarts2=restarts;
         out.retries2=retries;
         out.gradit2=iterations;
         con=con+1;  
       end      
   elseif con==0
       if dgf<gradval
         out.gradf1=minv;
         out.crit1=gradval;
         out.gradnorm1=dgf;
         out.gradeval1=evaluations;
         out.restarts1=restarts;
         out.retries1=retries;
         out.gradit1=iterations;
         con=con+1;    
       end
   end
   catch
   end
   iterations=iterations+1;
   
   %steepest descent
   if restart==1
       restarts=restarts+1;
       desc=v;
       c=dgf;
       dgf=-dgf;
       if c>1
          desc=desc/c;
       else
          c=1;
       end
       [p,vold,k]=LS(p,desc,maxstep);
       if listlength==1
           Delta={c*vold};
       else
           Delta={h*k*vold};
       end
       
       
       
       gold=c*vold;
       Gamma={};
       if listlength==0 %Variable to ensure steepest descent if steepest=1
           restart=1;
       else
           restart=0;
       end
    continue
   end
   %if not steepest do
   
   if listlength==1
       Gamma{1}=gold;
   else
       Gamma{end+1}=-v+gold;
   end
   
   
   desc1=descend(v,Delta,Gamma);
   desc=PrTpM(desc1,hg(p));
   dgf=round(desc.'*(-v),rto); % -v=grad
   
   if dgf>0
       restart=1;
       iterations=iterations-1;
       continue
   else
       c=norm(desc);
       if c>1
        desc=desc/c;
       end
   end
   
   %output
   
   if length(Delta)==listlength
       Gamma=Gamma(2:end);
       Delta=Delta(2:end);
   end
   
   Gamma{end+1}=v;
   
   [p,vold,k]=Lad(p,desc,maxstep);
   if listlength==1
      Delta={c*vold};
   else
      Delta{end+1}=h*k*vold;
   end
   
   gold=Gamma{end};
   Gamma=Gamma(1:end-1);
end

end

function [out] = stpoint(p)
global dim s1 s2 s3 SET row col
if SET=="Circ" || SET=="4Circ"
    s1=0;
    if norm(p)==0
        out=zeros([dim 1])+1;
        out=out/norm(out);
    else
        out=p/norm(p)  ;
    end
    s1=hf(out);  
elseif SET=="CircInt"
    s1=0;
    s2=0;
    s3=0;
    if norm(p)==0
        out=zeros([dim 1])+1;
        out=out/norm(out);
    else
        out=p/norm(p)  ;
    end
    s=hf(out);
    s1=s(1);
    s2=s(2);
    s3=s(3);
elseif SET=="Stiefel"
    if dim==48
        out=[
        [1;-1;0;0;0;0;0;0]/sqrt(2),
        [0;0;1;-1;0;0;0;0]/sqrt(2),
        [1;1;1;1;1;1;1;1]/sqrt(8),
        [0;0;0;0;1;-1;0;0]/sqrt(2),
        [0;0;0;0;0;0;1;-1]/sqrt(2),
        [1;1;-1;-1;0;0;0;0]/sqrt(4)        
        ];
    else
        M=zeros(row,col);
        M(1:col,1:col)=eye(col);
        out=Vec(M);    
    end
    
end  
end

function [out]=Mat(x)
    global row col,
    out=zeros(row,col);
    for i=(1:col)-1
        out(1:row,i+1)=x(i*row+1:(i+1)*row);
    end
end

function [out]=Vec(x)
global dim col row;
    out=zeros(dim,1);
    for j=(1:col)-1
        out(j*row+1:(j+1)*row)=x(1:row,j+1);
    end
end

function [out]= hL(x)
global evaluations
evaluations=evaluations+1;
out=cutest_obj(x);
end

function [out] = hf(x)
global SET s1 s2 s3 dim col
if SET=="Circ"
    out=x.'*x-s1;
elseif SET=="4Circ"
    out=-s1;
    for i=1:dim
        out=out+x(i)^4;
    end
elseif SET=="CircInt"
    out=zeros([3 1]);
    out(1)=x(1)^4+x(2)^4-s1;
    out(2)=norm(x)^2-s2;
    out(3)=x(1)^3*x(dim-1)^3-x(dim)-s3;
elseif SET=="Stiefel"
    out=Mat(x).'*Mat(x)-eye(col);
end
end

function [out] = hg(x)
global SET dim m row col
if SET=="Circ"
    out=2*x;
elseif SET=="4Circ"
    out=zeros([dim 1]);
    for i=1:dim
        out(i)=4*x(i)^3; 
    end    
elseif SET=="CircInt"
    out=zeros([dim m]);

    out(1,1)=4*x(1)^3;
    out(2,1)=4*x(2)^3;
    out(1:dim,2)=2*x;
    out(1,3)=3*x(1)^2*x(dim-1)^3;
    out(dim-1,3)=3*x(1)^3*x(dim-1)^2;
    out(dim,3)=-1;
elseif SET=="Stiefel"
    out=zeros([dim m]);
    M=Mat(x);
    for i=(1:col)
        for j=1:i
            temp=zeros(row,col);
            if j==i
                temp(1:row,i)=2*M(1:row,i);
            else
                temp(1:row,i)=M(1:row,j);
                temp(1:row,j)=M(1:row,i);
            end
            out(1:dim,(i-1)*i/2+j)=Vec(temp);
        end
    end
end 
end

function [out] = hH(x)
global SET dim m col row
if SET=="Circ"
    out{1}=2*eye(dim);
elseif SET=="4Circ"
    v=zeros([dim 1]);
    for i=1:dim
        v(i)=12*x(i)^2;
    end
    out{1}=diag(v);
elseif SET=="CircInt"
    H=cell([m 1]);
    Mat=zeros(dim);
    H{1}=Mat;
    H{3}=Mat;
    H{1}(1,1)=12*x(1)^2;
    H{1}(2,2)=12*x(2)^2;
    H{2}=2*eye(dim);
    H{3}(1,1)=6*x(1)*x(dim-1)^3;
    H{3}(1,dim-1)=9*x(1)^2*x(dim-1)^2;
    H{3}(dim-1,1)=9*x(1)^2*x(dim-1)^2;
    H{3}(dim-1,dim-1)=6*x(1)^3*x(dim-1);
    out=H;
elseif SET=="Stiefel"
    H=cell([m,1]);
    for i=1:col
        for j=1:i
            k=(i-1)*i/2+j;
            M=zeros(dim);
            if i==j
                M(row*(i-1)+1:row*i,row*(i-1)+1:row*i)=2*eye(row);
            else
                M(row*(i-1)+1:row*i,row*(j-1)+1:row*j)=eye(row);
                M(row*(j-1)+1:row*j,row*(i-1)+1:row*i)=eye(row);
            end
            H{k}=M;
        end
    end
    out=H;   
end
end

function [out]=HessCalc(x,i,j)
global SET dim 
out=zeros(dim);
if SET=="Circ" || SET=="Stiefel"    
elseif SET=="4Circ"
    out=zeros(dim);
    out(j,j)=24*x(j);
elseif SET=="CircInt"
    if i==1 && j==1
        out(1,1)=24*x(1);
    elseif i==1 && j==2
        out(2,2)=24*x(2);
    elseif i==3 && j==1
        out(1,1)=6*x(dim-1)^3;
        out(1,dim-1)=18*x(1)*x(dim-1)^2;
        out(dim-1,1)=out(1,dim-1);
        out(dim-1,dim-1)=18*x(1)^2*x(dim-1);
    elseif i==3 && j==dim-1
        out(1,1)=18*x(1)*x(dim-1)^2;
        out(1,dim-1)=18*x(1)^2*x(dim-1);
        out(dim-1,1)=out(1,dim-1);
        out(dim-1,dim-1)=6*x(1)^3;
    end
end



    %Hessian of d/dx_j of f_i

    
end

function [q]=descend(z,Delta,Gamma)
global listlength
if listlength==1
    beta=z.'*(z-Gamma{1})/norm(Gamma{1});
    q=z+beta*Delta{1};
else
    q=z;
    l=length(Delta);
    alpha=zeros(1,l);
    initscale=1;
    for i=flip(1:l)
        alpha(i)=Delta{i}.'*q/(Gamma{i}.'*Delta{i});
        q=q-alpha(i)*Gamma{i};
    end
    q=initscale^(-1)*q;
    for i=1:l
        beta=Gamma{i}.'*q/(Gamma{i}.'*Delta{i});
        q=q+(alpha(i)-beta)*Delta{i};
    end
end
end

function [out]= c(P,V)
global   minv maxstep restart iterations;
%variables

%intern variables
l=length(P);
CompVal=inf;
for i=1:l
    k=hL(P{i});
    if i==l
        out=i;
        CompVal=k;
        break
    elseif k<CompVal 
        CompVal=k;
    else
        out=i-1;
        break
    end
end
minv=CompVal;
if out==1 && restart==1
    iterations=5000;
elseif out==1
    restart=1;
end
    
    
if  0.8*maxstep <= out
   maxstep=3*maxstep;
   %fprintf('0.8max< %d , maxstep x3 \n',out);
elseif 0.2*maxstep <= out
   %fprintf('0.2max<= %d < 0.8max , maxstep 1x \n',out);
else
   %fprintf('10< %d<0.2max , maxstep/2 \n',out);
   maxstep=floor(maxstep/2);
end
end

function [out,out2] =Jacob(p,v)
%calculate the Jacobian of gamma''=f(gamma,gamma') in regard to gamma
global dim h m

A=hg(p); %is a jacobian


[Q,R]=qr(A);
R=R(1:m,1:m).';
Rin=inv(R);
J=zeros(dim);

H=hH(p);

w=zeros([m 1]);
for i=1:m
    w(i)=v.'*H{i}*v;
end
w=R\w;
for i=1:dim 
    %dQ,dR
    dA=zeros([dim m]);
    temp=zeros([dim 1]);
    temp(i)=1;
    for j=1:m
       dA(1:dim,j)=H{j}*temp; 
    end
    Q0=Q(1:dim,1:m);
    QP=Q(1:dim,m+1:dim);
    dK=QP.'*dA*Rin;
    A=Q0.'*dA*Rin;
    dRR=zeros(m);
    temp=A+A.';
    for s=1:m
       for t=s+1:m
          dRR(s,t)=temp(s,t);
       end
    end
    for s=1:m
       dRR(s,s)=A(s,s)/2; 
    end
    dOm=A-dRR;
    dQ=Q0*dOm+QP*dK;
    dR=dRR*R;
    %dQ till here
    t1=zeros([m 1]);
    for j=1:m
        t1(j)=v.'*HessCalc(p,i,j)*v;
    end
    t1=Q0*R*t1;
    t2=R\ (dR*w);
    J(1:dim,i)=(dQ*w-Q0*t2+t1);
    
end
out=J;

M=zeros([m dim]);
for i=1:m
   temp=H{i}*v;
   M(i,1:dim)=temp.';
end
out2=2*Q0*Rin*M;

end

function [out,out2,out3] = Lad(p0,v,n)
%PARATRANS Summary of this function goes here
%   Detailed explanation goes here
global Gamma Delta restart;
%variables 
logsteps=20;
tleng=100;
%internvariables
counter=1;
stop=0;

[P,V,T]=AB5se(p0,v);
while length(P)<n
    [P,V,T]= AB5st(P,V,T);
end  

k=c(P,V);
out=P{k};
out2=V{k};
out3=k;

if restart==1
   return 
end

[~,m]=size(Gamma);
Q=cell(1,m);
ParVec=Gamma;
Norm=zeros(length(Gamma));
for i=1:length(Gamma)
    Norm(i)=norm(Gamma{i});
end

for i=1:m
    Q{i}=LS(p0,ParVec{i}/Norm(i),logsteps,1);
end

while 1
    tsteps=tleng*counter;
    if tsteps>=k
        tsteps=k;
        stop=1;
    end

for i=1:m
    mid=floor((tsteps+tleng*(counter-1))/2);
	Q{i}=Log(Q{i},P{mid},2);
    if Q{i}==inf
            out=p0;
            return
    end
end  

if stop==1
    for i=1:m
        Q{i}=(-1)^counter*Norm(i)*Log(P{k},Q{i},1); % which sign? overwriting points which direction
        if Q{i}==inf
            out=p0;
            return
        end
    end
    Gamma=Q;
    break
end
counter=counter+1; 
end

[~,m]=size(Delta);
Q=cell(1,m);
ParVec=Delta;
counter=1;
stop=0;
for i=1:length(Delta)
    Norm(i)=norm(Delta{i});
end

for i=1:m
    Q{i}=LS(p0,ParVec{i}/Norm(i),logsteps,1);
end

while 1
    tsteps=tleng*counter;
    if tsteps>=k
        tsteps=k;
        stop=1;
    end
    
    for i=1:m
        mid=floor((tsteps+tleng*(counter-1))/2);
        Q{i}=Log(Q{i},P{mid},2);
        if Q{i}==inf
            out=p0;
            return
        end
    end  

    if stop==1
        for i=1:m
            Q{i}=(-1)^counter*Norm(i)*Log(P{k},Q{i},1); % which sign? overwriting points which direction
            if Q{i}==inf
                out=p0;
                return
            end
        end
        Delta=Q;
        return
    end
counter=counter+1; 
end

end

function [out,out2,out3] = LS(p0,v,n,check)
%LS Calculate the geodesic with starting points p0 and initial velocity v
%in TpM, for up to n steps/ at t=n*h



[P,V,T]=AB5se(p0,v);
while length(P)<n
    [P,V,T]= AB5st(P,V,T);
end    

try
    if check==1
        out=P{end};
        out2=V{end};
    end
catch
            k=c(P,V);
            out=P{k};
            out2=V{k};
            out3=k;
end

end

function [P,V,T] = AB5se2(p0,v0)
%RK5Se  Setup to a 5 step ODE solver, with Euler and Adam-Bashford 2 to 4 step algorithm
%p0 starting point, v0 initial velocity, h steplength
global rto
h=10^(-3);
%Euler 
p1=p0+h*v0;
t0=Mflog(p0,v0);
v1=round(v0+h*t0,rto);

%AdamBash2
temp=3/2*v1-1/2*v0;
p2=round(p1+h*temp,rto);
t1=Mflog(p1,v1);
temp=3/2*t1-1/2*t0;
v2=v1+h*temp;

%AB3
temp=23/12*v2-4/3*v1+5/12*v0;
p3=round(p2+h*temp,rto);
t2=Mflog(p2,v2);
temp=23/12*t2-4/3*t1+5/12*t0;
v3=v2+h*temp;
%AB4
temp=55/24*v3-59/24*v2+37/24*v1-3/8*v0;
p4=round(p3+h*temp,rto);
t3=Mflog(p3,v3);
temp=55/24*t3-59/24*t2+37/24*t1-3/8*t0;
v4=v3+h*temp;

Pr=round([p0,p1,p2,p3,p4],rto);
Vr=round([v0,v1,v2,v3,v4],rto);
Tr=round([t0,t1,t2,t3],rto);

[n,m]=size(p0);

P={Pr(1:n,1:m),Pr(1:n,m+1:2*m),Pr(1:n,2*m+1:3*m),Pr(1:n,3*m+1:4*m),Pr(1:n,4*m+1:5*m)};
V={Vr(1:n,1:m),Vr(1:n,m+1:2*m),Vr(1:n,2*m+1:3*m),Vr(1:n,3*m+1:4*m),Vr(1:n,4*m+1:5*m)};
T={Tr(1:n,1:m),Tr(1:n,m+1:2*m),Tr(1:n,2*m+1:3*m),Tr(1:n,3*m+1:4*m)};
end

function [P,V,T] = AB5st2(P,V,T)
%Adam-Bashford 5 step
global rto;
h=10^(-3);
[~,j]=size(P);

t=Mflog(P{j},V{j});

p_x=1901/720*V{j}-2774/720*V{j-1}+2616/720*V{j-2}-1274/720*V{j-3}+251/720*V{j-4};
v_x=1901/720*t-2774/720*T{j-1}+2616/720*T{j-2}-1274/720*T{j-3}+251/720*T{j-4};
p_new=round(P{j}+h*p_x,rto);
v_new=round(V{j}+h*v_x,rto);

P{j+1}=p_new;
T{j}=t;
V{j+1}=v_new;
end

function [out] = Log(p,q,lcase)
%Log Solves the logarithmic problem, finding v in TpM such that the
%geodesic with starting point p and initial velocity v hits q at t=1
global dim  rto retries h;
hint=0.001;
v=q-p;
failcounter=0;
while 1
    v=round(PrTpM(v,hg(p)),rto);
    [P,V,T]=AB5se2(p,v);
    while length(P)<1/hint
    [P,V,T]= AB5st2(P,V,T);
    end
    if norm(P{end}-q)<h/2
       % fprintf('log succesfull \n');
       break 
    else
        retries=retries+1;
        if failcounter==10
             out=inf;
             return
        else
            failcounter=failcounter+1;
        end
        
        %fprintf('log failed, retry \n');
        p0(1:dim,1)=p;
        p0(1:dim,2:dim+1)=zeros(dim);
        v0(1:dim,1)=v;
        v0(1:dim,2:dim+1)=eye(dim);
        [PM,VM,TM]=AB5se(p0,v0);
        for i=length(PM)+1:length(P)
            [PM,VM,TM]= AB5st(PM,VM,TM);
            X=PM{i};
            X(1:dim,1)=P{i};
            PM{i}=X;
            X=VM{i};
            X(1:dim,1)=V{i};
            VM{i}=X;
        end       
        M=PM{end}; %Calculate Corrector by logsequence
        M=M(1:dim,2:dim+1);
        w=P{end}-q;
        cor=M\w;
        v=v-cor;
       
    end
end
if lcase==1
   out=v/norm(v);
   return
elseif lcase==2
   while length(P)<2/hint
      %inner loop -> geodesic for specific v
      [P,V,T]= AB5st2(P,V,T);
   end
   out=P{end};
   out=out(1:dim,1);
   return
end
end

function [P,V,T] = AB5se(p0,v0)
%RK5Se  Setup to a 5 step ODE solver, with Euler and Adam-Bashford 2 to 4 step algorithm
%p0 starting point, v0 initial velocity, h steplength
global h rto
%Euler 
p1=p0+h*v0;
t0=Mflog(p0,v0);
v1=round(v0+h*t0,rto);

%AdamBash2
temp=3/2*v1-1/2*v0;
p2=round(p1+h*temp,rto);
t1=Mflog(p1,v1);
temp=3/2*t1-1/2*t0;
v2=v1+h*temp;

%AB3
temp=23/12*v2-4/3*v1+5/12*v0;
p3=round(p2+h*temp,rto);
t2=Mflog(p2,v2);
temp=23/12*t2-4/3*t1+5/12*t0;
v3=v2+h*temp;
%AB4
temp=55/24*v3-59/24*v2+37/24*v1-3/8*v0;
p4=round(p3+h*temp,rto);
t3=Mflog(p3,v3);
temp=55/24*t3-59/24*t2+37/24*t1-3/8*t0;
v4=v3+h*temp;

Pr=round([p0,p1,p2,p3,p4],rto);
Vr=round([v0,v1,v2,v3,v4],rto);
Tr=round([t0,t1,t2,t3],rto);

[n,m]=size(p0);

P={Pr(1:n,1:m),Pr(1:n,m+1:2*m),Pr(1:n,2*m+1:3*m),Pr(1:n,3*m+1:4*m),Pr(1:n,4*m+1:5*m)};
V={Vr(1:n,1:m),Vr(1:n,m+1:2*m),Vr(1:n,2*m+1:3*m),Vr(1:n,3*m+1:4*m),Vr(1:n,4*m+1:5*m)};
T={Tr(1:n,1:m),Tr(1:n,m+1:2*m),Tr(1:n,2*m+1:3*m),Tr(1:n,3*m+1:4*m)};
end

function [P,V,T] = AB5st(P,V,T)
%Adam-Bashford 5 step
global h rto;
[~,j]=size(P);

t=Mflog(P{j},V{j});

p_x=1901/720*V{j}-2774/720*V{j-1}+2616/720*V{j-2}-1274/720*V{j-3}+251/720*V{j-4};
v_x=1901/720*t-2774/720*T{j-1}+2616/720*T{j-2}-1274/720*T{j-3}+251/720*T{j-4};
p_new=round(P{j}+h*p_x,rto);
v_new=round(V{j}+h*v_x,rto);

P{j+1}=p_new;
T{j}=t;
V{j+1}=v_new;
end

function [w]=Mflog(qin,vin)
%evaluation of geodesic ODE with option to log ODE
[n,m]=size(qin);
q=qin(1:n,1);
v=vin(1:n,1);
if m==1  
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
        H(i)=v.'*eH{i}*v; % changed from eH
    end
    t=R\H;
    w=-Q(1:n,1:m)*t;
else
    qm=qin(1:n,2:m);
    vm=vin(1:n,2:m);
    [Jgamma,Jdgamma]=Jacob(q,v);
    temp=(Jgamma.'*qm+Jdgamma.'*vm);
    w=[zeros([n,1]),temp];
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
end
function [fitness]=ErPSOGainFlatness(par)
K=2;                
I1=601;            
H1=50;              
   
%lamdas=1610;
%m=(lamdas-1510)/0.8;

p1=par.p1;           %6个泵浦光功率    
p2=par.p2;
p3=par.p3;
p4=par.p4;
p5=par.p5;               
p6=par.p6;

v1=par.v1;           %6个泵浦光波长 
v2=par.v2;
v3=par.v3;
v4=par.v4;
v5=par.v5;
v6=par.v6;

%n=6;
p7=0.00001;     
v7=198.675;


%l=m+n;


for j=1:106
    if(j==1)
        v(j)=v1*1e12;
        p(j)=p1;
    elseif(j==2)
        v(j)=v2*1e12;
        p(j)=p2;
    elseif (j==3)
        v(j)=v3*1e12;
        p(j)=p3;
    elseif (j==4)
        v(j)=v4*1e12;
        p(j)=p4;
    elseif(j==5)
        v(j)=v5*1e12;
        p(j)=p5;
    elseif(j==6)
        v(j)=v6*1e12;
        p(j)=p6;
    elseif(j==7)
        v(j)=v7*1e12;
        p(j)=p7;
    else
        v(j)=(v7-0.124*(j-7))*1e12;
        p(j)=p7;
    end
end

%%泵浦波、信号波损耗系数%%
for j=1:6
    alf(j)=0.231045/(1000*4.343);%泵浦波损耗系数(单位1/m)  %/(1000*4.343)将单位dB/km换算成1/m
end
for j=7:106
    alf(j)=0.190272/(1000*4.434);%信号波损耗系数（单位1/m）
end

%%增益系数函数%%
gain=zeros(106,106);
x=zeros(106,106);
for j=1:106
    for k=1:106
        b=1e12;
        x(j,k)=abs(v(j)-v(k));
        x(j,k)=x(j,k)/b;
        c=exp(-(x(j,k)-1.686)*(x(j,k)-1.686)/(2*1.023146*1.023146));
        d=exp(-(x(j,k)-3.039)*(x(j,k)-3.039)/(2*1.873365*1.873365));
        e=exp(-(x(j,k)-5.553)*(x(j,k)-5.553)/(2*3.19248*3.19248));
        f1=exp(-(x(j,k)-10.821)*(x(j,k)-10.821)/(2*5.575755*5.575755));
        s=exp(-(x(j,k)-13.71)*(x(j,k)-13.71)/(2*3.61371*3.61371));
        t=exp(-(x(j,k)-14.664)*(x(j,k)-14.664)/(2*0.625194*0.625194));
        o=exp(-(x(j,k)-17.82)*(x(j,k)-17.82)/(2*1.41888*1.41888));
        q=exp(-(x(j,k)-19.743)*(x(j,k)-19.743)/(2*4.68895*4.68895));
        gain(j,k)=3.68*1e-14*(0.077*c+0.149*d+0.171*e+0.615*f1+0.655*s+0.185*t+0.207*o+0.126*q);
    end
end
%%光纤有效截面积%%
Aeff=zeros(1,106);
x1=zeros(1,106);
for j=1:106
    x1(j)=v(j);%x1为频率，单位Thz
    x1(j)=x1(j)/b;
    Aeff(j)=(-0.5541*x1(j)+184.02)*1e-12;%光纤有效截面积和频率之间的关系表达式
end
%%非线性耦合方程PiPj的系数%%
g=zeros(106,106);
for j=1:106
    for k=1:106
        if(v(k)>v(j))
                g(j,k)=gain(j,k)/(K*Aeff(j));
        end
        if(v(k)<v(j))
        g(j,k)=-v(j)*gain(j,k)/(v(k)*K*Aeff(j));
        end
    end
end
%%将初始功率值赋给y(i,k),i代表迭代次数，此时i=1%%
y1=zeros(I1,106);%I1是迭代点数
y1(1,:)=p(1:106);
f1=zeros(I1-1,7,106);
%%龙格-库塔求前三项%%
for i=1:I1-1
   for j=1:106
        sum1=0;
        sum2=0;
        sum3=0;
        sum4=0;
        sum5=0;
        sum6=0;
        sum7=0;
        sum8=0;
       for k=1:j-1
            sum1=sum1+g(j,k)*y1(i,j)*y1(i,k);%当频率zk>zj时，将k=1:j-1所有项求和
        end
        for k=j+1:106
             sum2=sum2+g(j,k)*y1(i,j)*y1(i,k);%当频率zk<zj时，将k=j+1:54所有项求和
        end
        f1(i,1,j)=(-alf(j)*y1(i,j))+sum1+sum2;
        for k=1:j-1
            sum3=sum3+g(j,k)*(y1(i,j)+H1*f1(i,1,j)/2)*y1(i,k);
        end
        for k=j+1:106
            sum4=sum4+g(j,k)*(y1(i,j)+H1*f1(i,1,j)/2)*y1(i,k);
        end
        f1(i,2,j)=(-alf(j)*(y1(i,j)+H1*f1(i,1,j)/2))+sum3+sum4;
        for k=1:j-1
            sum5=sum5+g(j,k)*(y1(i,j)+H1*f1(i,2,j)/2)*y1(i,k);
        end
        for k=j+1:106
            sum6=sum6+g(j,k)*(y1(i,j)+H1*f1(i,2,j)/2)*y1(i,k);
        end
        f1(i,3,j)=(-alf(j)*(y1(i,j)+H1*f1(i,2,j)/2))+sum5+sum6;
        for k=1:j-1
            sum7=sum7+g(j,k)*(y1(i,j)+H1*f1(i,3,j))*y1(i,k);
        end
        for k=j+1:106
            sum8=sum8+g(j,k)*(y1(i,j)+H1*f1(i,3,j))*y1(i,k);
        end
        f1(i,4,j)=(-alf(j)*(y1(i,j)+H1*f1(i,3,j)))+sum7+sum8;
        y1(i+1,j)=y1(i,j)+(f1(i,1,j)+2*f1(i,2,j)+2*f1(i,3,j)+f1(i,4,j))*H1/6;
   end
end
    
%%龙格-库塔法解方程--无泵浦光时，信号光的功率变化%%
y2=zeros(I1,106);
y2(1,7:106)=p(7:106);
f2=zeros(I1-1,7,j);
for i=1:I1-1
    for j=7:106
        sum1=0;
        sum2=0;
        sum3=0;
        sum4=0;
        sum5=0;
        sum6=0;
        sum7=0;
        sum8=0;
        for k=7:j-1
            sum1=sum1+g(j,k)*y2(i,j)*y2(i,k);%当频率zk>zj时，将k=1:j-1所有项求和
        end
        for k=j+1:106
            sum2=sum2+g(j,k)*y2(i,j)*y2(i,k);%当频率zk<zj时，将k=j+1:54所有项求和
        end
        f2(i,1,j)=(-alf(j)*y2(i,j))+sum1+sum2;
        for k=7:j-1
            sum3=sum3+g(j,k)*(y2(i,j)+H1*f2(i,1,j)/2)*y2(i,k);
        end
        for k=j+1:106
            sum4=sum4+g(j,k)*(y2(i,j)+H1*f2(i,1,j)/2)*y2(i,k);
        end
        f2(i,2,j)=(-alf(j)*(y2(i,j)+H1*f2(i,1,j)/2))+sum3+sum4;
        for k=7:j-1
            sum5=sum5+g(j,k)*(y2(i,j)+H1*f2(i,2,j)/2)*y2(i,k);
        end
        for k=j+1:106
            sum6=sum6+g(j,k)*(y2(i,j)+H1*f2(i,2,j)/2)*y2(i,k);
        end
        f2(i,3,j)=(-alf(j)*(y2(i,j)+H1*f2(i,2,j)/2))+sum5+sum6;
        for k=7:j-1
            sum7=sum7+g(j,k)*(y2(i,j)+H1*f2(i,3,j))*y2(i,k);
        end
        for k=j+1:106
            sum8=sum8+g(j,k)*(y2(i,j)+H1*f2(i,3,j))*y2(i,k);
        end
        f2(i,4,j)=(-alf(j)*(y2(i,j)+H1*f2(i,3,j)))+sum7+sum8;
        y2(i+1,j)=y2(i,j)+(f2(i,1,j)+2*f2(i,2,j)+2*f2(i,3,j)+f2(i,4,j))*H1/6;
    end
end

lamda=zeros(1,100);
gain_=zeros(1,100);
for j=1:100
    lamda(j)=(3*1e17)/v(j+6);%用波长表示,lamda=c/v,这里单位是nm，所以是17次方
    gain_(j)=10*log10(y1(I1,j+6)/y2(I1,j+6));%开关增益dB
end
 gain=0;
 for j=1:100
     gain=gain+gain_(j);
 end
 gain1=gain/100;
aa=max(gain_(1,:));
bb=min(gain_(1,:));
fitness=aa-bb;

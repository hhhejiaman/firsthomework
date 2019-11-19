clear all           %%�źŹⲨ��1530nm-1570nm�����0.8nm%%     %%���ϵ����4.343��Ϊ4.434%%
clc
K=2;                %Keff=2,��ʱ���ֹ���źŹ��ƫ��̬���
I1=601;             %Ҫ�����ĵ����������޸�
H1=50;             %����-�������ļ��㲽��
% H2=1000;          %��Admas���ļ��㲽��   %���˳���L=H1*(I1-1)  
p1=1.1070;            %���ֲ�1406nm�Ĺ���
p2=0.6606;            %���ֲ�1413nm�Ĺ���
p3=0.2520;            %���ֲ�1445nm�Ĺ���
p4=0.1001;            %���ֲ�1492nm�Ĺ���
p5=0.1326;
p6=0.0636;
p7=0.00001;            %�źŹ�Ĺ���

v1=227.5554;           %���ֲ�1406nm��Ƶ��
v2=228.1089;           %���ֲ�1413nm��Ƶ��
v3=213.1949;           %���ֲ�1445nm��Ƶ��
v4=201.1772;           %���ֲ�1492nm��Ƶ��
v5=210.5423;           %���ֲ�1356nm��Ƶ��
v6=203.7131;             %�źŲ�1510nm��Ƶ��
v7=198.675;
%%����Ƶ�ʺ͹���  %k=1:54,54Ϊ���ֲ����źŲ�����·����kȡһ��ֵ��������Ӧ��Ƶ�ʺ͹��ʣ�Ƶ�ʰ���������%%
for j=1:106  
    if(j==1)
        v(j)=v1*1e12;%���ֲ�1427nm��Ƶ��,��λHZ
        p(j)=p1;
    elseif(j==2)
        v(j)=v2*1e12;%���ֲ�1445nm��Ƶ��
        p(j)=p2;
    elseif (j==3)
        v(j)=v3*1e12;%���ֲ�1466nm��Ƶ��
        p(j)=p3;
    elseif (j==4)
        v(j)=v4*1e12;%�źŲ�1530nm��Ƶ��
        p(j)=p4;
    elseif (j==5)
        v(j)=v5*1e12;%�źŲ�1530nm��Ƶ��
        p(j)=p5;
    elseif (j==6)
        v(j)=v6*1e12;%�źŲ�1530nm��Ƶ��
        p(j)=p6;
    elseif (j==7)
        v(j)=v7*1e12;%�źŲ�1530nm��Ƶ��
        p(j)=p7;
    else
        v(j)=(v7-0.124*(j-7))*1e12;%�����źŲ�Ƶ��
        p(j)=p7;
    end
end

%%���ֲ����źŲ����ϵ��%%
for j=1:6
alf(j)=0.231045/(1000*4.434);%���ֲ����ϵ��(��λ1/m)  %/(1000*4.343)����λdB/km�����1/m
end
for j=7:106
alf(j)=0.190272/(1000*4.434);%�źŲ����ϵ������λ1/m��
end

%%����ϵ������%%
gain=zeros(106,106);
x=zeros(106,106);
for j=1:106
    for k=1:106
        b=1e12;
        x(j,k)=abs(v(j)-v(k));%xΪƵ��������λΪThz
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
%%������Ч�����%%
Aeff=zeros(1,106);
x1=zeros(1,106);
for j=1:106
    x1(j)=v(j);%x1ΪƵ�ʣ���λThz
    x1(j)=x1(j)/b;
    Aeff(j)=(-0.5541*x1(j)+184.02)*1e-12;
end
%%��������Ϸ���PiPj��ϵ��%%
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
%%����ʼ����ֵ����y(i,k),i���������������ʱi=1%%
y1=zeros(I1,106);
y1(1,:)=p(1:106);
f1=zeros(I1-1,7,106);
%%����-������ǰ����%%
for i=1:I1-1
    for j=1:105
        sum1=0;
        sum2=0;
        sum3=0;
        sum4=0;
        sum5=0;
        sum6=0;
        sum7=0;
        sum8=0;
        for k=1:j-1
            sum1=sum1+g(j,k)*y1(i,j)*y1(i,k);%��Ƶ��zk>zjʱ����k=1:j-1���������
        end
        for k=j+1:106
             sum2=sum2+g(j,k)*y1(i,j)*y1(i,k);%��Ƶ��zk<zjʱ����k=j+1:54���������
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
%%����-�������ⷽ��--�ޱ��ֹ�ʱ���źŹ�Ĺ��ʱ仯%%
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
            sum1=sum1+g(j,k)*y2(i,j)*y2(i,k);%��Ƶ��zk>zjʱ����k=1:j-1���������
        end
        for k=j+1:106
            sum2=sum2+g(j,k)*y2(i,j)*y2(i,k);%��Ƶ��zk<zjʱ����k=j+1:54���������
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
  

% %%�����������Ϸ���
% L=zeros(1,I1);
 %for j=1:I1
 %   L(j)=0.05*(j-1);
%end
%plot(L,y1(:,1:3)*1000);
%axis([0 30 0 20])
%axis([0 30 0 1000]);
%title('һ�ױ��ֹ⹦���ع��˵ı仯');
%xlabel('L/km');
%ylabel('pump power/mw')

lamda=zeros(1,100);
gain_=zeros(1,100);
for j=1:100
    lamda(j)=(3*1e17)/v(j+6);%�ò�����ʾ
    gain_(j)=10*log10(y1(I1,j+6)/y2(I1,j+6));%��������dB
end
%aa=max(gain_(1,:));
%bb=min(gain_(1,:));
%cc=aa-bb;
plot(lamda,gain_,'r-');
axis([1510 1610 10 25]);
title('����������');
xlabel('lamda/nm');
ylabel('gain/dB');



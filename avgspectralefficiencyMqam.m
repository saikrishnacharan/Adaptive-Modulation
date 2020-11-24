clc;
close all;
clear all;
m_dB=15:1.5:30;% average snr range
m=10.^(m_dB/10);%linear average snr
N=length(m);
M=zeros(1,N);%for Pb=0.001
M1=zeros(1,N);%for Pb=0.000001
n=1e3;%no of snr values generated in each snr
Pb=0.001;%Pb value
Pb1=0.000001;%pb value
K=-1.5/log(5*Pb);%loss function depending on the Pb value
K1=-1.5/log(5*Pb1);%loss function depending on the Pb1 value
for k=1:1:N
    tem=0;
    for q=1:1:1e2%genrating snr values 1e2 time and average the spectral efficiency  
        H=exprnd((m(k)),1,n);
        y=sort(H);
        P=(1/m(k)).*exp(-y./m(k));
        for i=n:-1:1
            Pt=0;
            p=0;
            for j=n:-1:n+1-i
                Pt=Pt+(P(j)/y(j));
                p=p+P(j);
            end
            gm=p/(1+Pt);
            if y(j)>= gm
                break
            end
        end
        yk=gm/K;
        f=@(y) log2(y./yk).*(1./m(k)).*exp(-y./m(k));
        tem=tem+integral(f,yk,Inf);
    end
    M(k)=sum(tem)/1e2;
end
%similar as above for Pb1 value 
for k=1:1:N
    tem=0;
    for q=1:1:1e2
        H=exprnd(m(k),1,n);
        y=sort(H);
        P=(1/m(k)).*exp(-y./m(k));
        for i=n:-1:1
            Pt=0;
            p=0;
            for j=n:-1:n+1-i
                Pt=Pt+(P(j)/y(j));
                p=p+P(j);
            end
            gm=p/(1+Pt);
            if y(j)>= gm
                break
            end
        end
        yk=gm/K1;
        f=@(y) log2(y./yk).*(1./m(k)).*exp(-y./m(k));
        tem=tem+integral(f,yk,Inf);
    end
    M1(k)=sum(tem)/1e2;
end
figure(1)
plot(m_dB,M);hold on;grid on;%plotiing for Pb=0.001
plot(m_dB,M1);%plotting for 0.000001
title("Average spectral efficiency in Rayleigh fading");
xlabel("Average SNR in dB");
ylabel("Spectral efficiency(bps/Hz)");
legend('Pb is 10-3','Pb is 10-6');

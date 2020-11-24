clc;
clear all;
%no of receiver and transmitter  antennas
Nr=1;
N=1e5;
%linear snr
snr_dB=15:1:30;
snr=10.^(snr_dB/10);
snr=2*snr;
temper=zeros(1,length(snr));
    for b=1:1:length(snr)
        tem=0;
        temt=0;
        temtt=0;
        for q=1:1:N
            H=sqrt(1/2).*(randn(Nr,Nr)+1i*randn(Nr,Nr));
            S=svd(H);%getting singular values
            K=(S).^2;%getting squares of the singular values
            r=rank(H);%finding rank of the matrix;
            for i=r:-1:1
                Pt=0;
                for j=1:1:i
                    temp=snr(b)*K(j);
                    Pt=Pt+(1/temp);
                end
                gm=i/(1+Pt);
                if snr(b)*K(i)>= gm
                    break
                end
            end
            for k=1:1:i
                P(k)=1/gm-1/(snr(b)*K(k));
                C(k)=log2(snr(b)*K(k)/gm);
            end
            for w=1:1:r
                CC(w)=log2(1+((snr(b)*K(w))/r));
            end
            %CC=log2(1+snr(b)*K(1));
            tem=tem+sum(CC);
        end
        temper(1,b)=tem/N;
    end

    title("capacity vs snr");
    xlabel("SNR in dB");
    ylabel("Capacity");grid on;
    plot(snr_dB,10*log10(temper(1,:)),'-s');

legend('Water filling theorem');
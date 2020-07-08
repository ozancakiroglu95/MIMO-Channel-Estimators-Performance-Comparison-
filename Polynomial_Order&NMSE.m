close all;
clear;
clc;

M = 50;
theta = 30*pi/180;
deviation = 10*pi/180;

R1 = Mas_MIMO_covariance_matrix(M, theta, deviation);
R2 = Mas_MIMO_covariance_matrix(M, 40*pi/180, deviation);
a = 2/trace(R1+eye(M));
A_peach = zeros(M,M,30);

SNR1dB = 10;
SNR1 = 10.^(SNR1dB/10);

SNR2dB = 10;
SNR2 = 10.^(SNR2dB/10);
SNR2 = 0;


for l1 = 1:50

    a = 2/trace(SNR1*R1+SNR2*R2+eye(M));
    A_peach(:,:,l1) = zeros(M,M);
    
    for l2 = 1:l1
    
        A_peach(:,:,l1) = A_peach(:,:,l1) + a*((eye(M)-a*(SNR1*R1+SNR2*R2+eye(M)))^(l2-1));
    
    end
    
    for i = 1:l1+1

        for j = 1:l1+1

            A(i,j)  =  trace(R1*sqrt(SNR1)*(SNR1*R1+SNR2*R2+eye(M))^(i+j-1)*R1*sqrt(SNR1));
            b(i)    =  trace(R1*sqrt(SNR1)*(SNR1*R1+SNR2*R2+eye(M))^(i-1)*R1*sqrt(SNR1));
             
        end

    end
   
    w = A\transpose(b);
    NMSE_WPeach(l1) = (trace(R1) + w'*A*w - transpose(b')*w - w'*transpose(b))/trace(R1);
    NMSE_WPeach(l1) = 10*log(abs(NMSE_WPeach(l1)));
    
    NMSE_Peach(l1) = (trace(R1) + trace(R1*sqrt(SNR1)*A_peach(:,:,l1)*(SNR1*R1+SNR2*R2+eye(M))*A_peach(:,:,l1)'*R1*sqrt(SNR1)) - 2*(trace(sqrt(SNR1)*R1*A_peach(:,:,l1)'*R1*sqrt(SNR1))))/trace(R1);
    NMSE_Peach(l1) = 10*log(abs(NMSE_Peach(l1)));
    
    NMSE_MMSE = real(trace(R1 - SNR1*R1*((SNR1*R1+SNR2*R2+eye(M))\R1)))/trace(R1);
    NMSE_MMSE = 10*log(abs(NMSE_MMSE));
    
    A_EWMMSE = (sqrt(SNR1)/(SNR1+SNR2+1))*eye(M);
    NMSE_EWMMSE = (trace(R1) + trace(A_EWMMSE*(SNR1*R1+SNR2*R2+eye(M))*A_EWMMSE') - 2*real(trace(A_EWMMSE'*R1))*sqrt(SNR1))/trace(R1);
    NMSE_EWMMSE = 10*log(abs(NMSE_EWMMSE));
    
    A_LS = eye(M)/sqrt(SNR1);
    NMSE_LS = (trace(R1) + trace(A_LS*(SNR1*R1+SNR2*R2+eye(M))*A_LS') - 2*real(trace(A_LS'*R1))*sqrt(SNR1))/trace(R1);
    NMSE_LS = 10*log(abs(NMSE_LS));
    
end

figure;
hold on; box on;


plot(1:50,NMSE_LS*ones(1,50),'Color',[0.6350 0.0780 0.1840],'LineWidth',3);
plot(1:50,NMSE_EWMMSE*ones(1,50),'Color',[0.9290 0.6940 0.1250],'LineWidth',3);
plot(1:50,NMSE_MMSE*ones(1,50),'Color',[0 0.4470 0.7410],'LineWidth',3);
plot(1:50,NMSE_Peach,'Color',[0.6789 0.4470 0.7410],'LineWidth',3)
plot(1:50,NMSE_WPeach,'LineWidth',3)

ylim([-42 -5])

title("One user without interferer",'FontSize', 15);
xlabel("Polynomial Order (L)",'FontSize', 15);
ylabel("Normalized MSE (dB)",'FontSize', 15);

dim = [.28 .1 .8 .8];

set(gca,'Color',[255/256 228/256 181/256])
legend('LS','EW-MMSE','MMSE','Peach','Weighted Peach','Location','NorthEast',"Color","w");
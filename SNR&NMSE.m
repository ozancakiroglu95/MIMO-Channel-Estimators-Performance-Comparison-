close all;
clear;

M = 50;
deviation = 10*pi/180;
theta1 = 30*pi/180;
theta2 = 40*pi/180;

SNR1dB = -10:1:20;
SNR1 = 10.^(SNR1dB/10);

SNR2dB = SNR1dB;
SNR2 = 10.^(SNR2dB/10);
SNR2 = 0*SNR2;

R1 = Mas_MIMO_covariance_matrix(M, theta1, deviation);
R2 = Mas_MIMO_covariance_matrix(M, theta2, deviation);

for s = 1:length(SNR1dB)
    
    a = 2/trace(SNR1(s)*R1+SNR2(s)*R2+eye(M));
    A_peach = zeros(M,M);
    
    for l = 1:20
    
        A_peach = A_peach + a*((eye(M)-a*(SNR1(s)*R1+SNR2(s)*R2+eye(M)))^(l-1));
    
    end
    
    NMSE_Peach(s) = (trace(R1) + trace(R1*sqrt(SNR1(s))*A_peach*(SNR1(s)*R1+SNR2(s)*R2+eye(M))*A_peach'*R1*sqrt(SNR1(s))) - 2*(trace(sqrt(SNR1(s))*R1*A_peach'*R1*sqrt(SNR1(s)))))/trace(R1);
    
    for i = 1:20+1

        for j = 1:20+1

            A(i,j)  =  trace(R1*sqrt(SNR1(s))*(SNR1(s)*R1+SNR2(s)*R2+eye(M))^(i+j-1)*R1*sqrt(SNR1(s)));
            b(i)    =  trace(R1*sqrt(SNR1(s))*(SNR1(s)*R1+SNR2(s)*R2+eye(M))^(i-1)*R1*sqrt(SNR1(s)));
             
        end

    end
   
    w = A\transpose(b);
    NMSE(s) = (trace(R1) + w'*A*w - transpose(b')*w - w'*transpose(b))/trace(R1);
    
    NMSE_MMSE(s) = real(trace(R1 - SNR1(s)*R1*((SNR1(s)*R1+SNR2(s)*R2+eye(M))\R1)))/trace(R1);
   
    A_EWMMSE = (sqrt(SNR1(s))/(SNR1(s)+SNR2(s)+1))*eye(M);
    NMSE_EWMMSE(s) = (trace(R1) + trace(A_EWMMSE*(SNR1(s)*R1+SNR2(s)*R2+eye(M))*A_EWMMSE') - 2*real(trace(A_EWMMSE'*R1))*sqrt(SNR1(s)))/trace(R1);
    
    A_LS = eye(M)/sqrt(SNR1(s));
    NMSE_LS(s) = (trace(R1) + trace(A_LS*(SNR1(s)*R1+SNR2(s)*R2+eye(M))*A_LS') - 2*real(trace(A_LS'*R1))*sqrt(SNR1(s)))/trace(R1);

end


figure;
hold on; box on;

plot(SNR1dB,NMSE_LS,'Color',[0.6350 0.0780 0.1840],'LineWidth',3);
plot(SNR1dB,NMSE_EWMMSE,'Color',[0.9290 0.6940 0.1250],'LineWidth',3);
plot(SNR1dB,NMSE_MMSE,'Color',[0 0.4470 0.7410],'LineWidth',3);
plot(SNR1dB,NMSE_Peach,'Color',[0.6789 0.4470 0.7410],'LineWidth',3);
plot(SNR1dB,NMSE,'LineWidth',3)

title("One User without interferer",'FontSize', 15)
xlabel('SNR[dB]');
ylabel('NMSE');
set(gca, 'YScale', 'log')
set(gca,'Color',[0.4 0.6 0.7])
legend('LS','EW-MMSE','MMSE','Peach','Weighted Peach','Location','NorthEast');
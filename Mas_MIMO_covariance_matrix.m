% Toeplitz matrix implementation idea was taken from emilbjornson Github implementation since
% this style provides least computational time compared to other implementation styles and easiness 
% for further performance evaluations.

function R = Mas_MIMO_covariance_matrix(M,theta,ASD)

firstRow = zeros(M,1);

for column = 1:M   
    distance = (column-1)/4;
    F = @(Delta)exp(1i*2*pi*distance*sin(theta+Delta)).*exp(-Delta.^2/(2*ASD^2))/(sqrt(2*pi)*ASD); 
    firstRow(column) = integral(F,-20*ASD,20*ASD);       
end
 
R = toeplitz(firstRow);



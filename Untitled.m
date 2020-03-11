clear all;
close all;
clc;


ITER = 200;
K = 10; %No. of users
Mv = 20:30:500; %No. of BS antennas
R = 1;
pu = zeros(ITER,length(Mv));
D = Dmatrix(K);
beta = diag(D);

for i = 1:ITER
    for M = 1:length(Mv)
        for k = 1:K
            nr = 2^R-1;
            dr = ((M-1)*beta(k))-(nr*(sum(beta)-beta(k)));
            pu_temp(k) = nr/dr;
        end
        pu(i,M) = pu_temp(5)/sum(pu_temp);
    end
end

pu_plot = log10(sum(pu)/ITER);
figure;
plot(Mv,pu_plot);
grid on;
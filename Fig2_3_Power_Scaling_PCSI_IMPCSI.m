clear all;
close all;
clc;

ITER = 200;
K = 10; %No. of users
Mv = 20:30:500; %No. of BS antennas
M = 10
Eu_dB = 20;
%Eu_dB = 5;
Eu = 10^(Eu_dB/10);

figure;

for i = 1:2
    
    [rate_MRC_PCSI,bound_MRC_PCSI,rate_ZF_PCSI,bound_ZF_PCSI,rate_MMSE_PCSI,bound_MMSE_PCSI] = PCSI(Eu,i);
    [rate_MRC_IPCSI,bound_MRC_IPCSI,rate_ZF_IPCSI,bound_ZF_IPCSI,rate_MMSE_IPCSI,bound_MMSE_IPCSI] = IPCSI(Eu,i);

    p1 = plot(Mv,rate_MRC_PCSI,'rsquare--','LineWidth',1.2);
    hold on;
    %plot(Mv,bound_MRC_PCSI,'square','MarkerEdgeColor','r','LineWidth',1);
    %hold on;
    p2 = plot(Mv,rate_ZF_PCSI,'b*--','LineWidth',1.2);
    hold on;
    %plot(Mv,bound_ZF_PCSI,'*','MarkerEdgeColor','b','LineWidth',1);
    %hold on;
    p3 = plot(Mv,rate_MMSE_PCSI,'o--','color',[0 0.5 0],'LineWidth',1.2);
    hold on;
    %plot(Mv,rate_MMSE_PCSI,'o','MarkerEdgeColor',[0 0.5 0],'LineWidth',1);
    %hold on;
    p4 = plot(Mv,rate_MRC_IPCSI,'rsquare-','LineWidth',1.2);
    hold on;
    %plot(Mv,bound_MRC_IPCSI,'square','MarkerEdgeColor','r','LineWidth',1);
    %hold on;
    p5 = plot(Mv,rate_ZF_IPCSI,'b*-','LineWidth',1.2);
    hold on;
    %plot(Mv,bound_ZF_IPCSI,'*','MarkerEdgeColor','b','LineWidth',1);
    %hold on;
    p6 = plot(Mv,rate_MMSE_IPCSI,'o-','color',[0 0.5 0],'LineWidth',1.2);
    hold on;
    %plot(Mv,rate_MMSE_IPCSI,'o','MarkerEdgeColor',[0 0.5 0],'LineWidth',1);
    
end

grid on;
grid minor;
%legend('Perfect CSI, MRC','Imperfect CSI, MRC','Perfect CSI, ZF','Imperfect CSI, ZF','Perfect CSI, MMSE','imperfect CSI, MMSE');
legend([p1 p4 p2 p5 p3 p6], 'Perfect CSI, MRC','Imperfect CSI, MRC','Perfect CSI, ZF','Imperfect CSI, ZF','Perfect CSI, MMSE','Imperfect CSI, MMSE');
xlabel('Number of BS Antennas (M)');
ylabel('Spectral Efficiency (bits/s/Hz)'); 
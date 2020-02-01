f_opt= 4.0645582;
ASGARD_fc_fl_rcv1= importdata('../IPALM/results/DLRCSGR2_FL_outer_rcv1tau_47236',' ');
DLRCSGR_eps_fl_rcv1= importdata('../IPALM/results/DLRCSGR3_FL_outer_rcv1tau_1',' ');
SMART_CD_fl_rcv1= importdata('../IPALM/results/PDCD_FL_outer_rcv1tau_1',' ');
KATYUSHA_fl_rcv1= importdata('../IPALM/results/ALM_FL_outer_rcv1tau_259',' ');
s1= size(ASGARD_fc_fl_rcv1,1);
s2= size(DLRCSGR_eps_fl_rcv1,1);
s3= size(SMART_CD_fl_rcv1,1);
s4= size(KATYUSHA_fl_rcv1,1);
plot(ASGARD_fc_fl_rcv1(1:3:s1,2),log10((abs(ASGARD_fc_fl_rcv1(1:3:s1,3)- f_opt)/f_opt)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_fl_rcv1(1:15:s2,2),log10((abs(DLRCSGR_eps_fl_rcv1(1:15:s2,3)- f_opt)/f_opt)),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_fl_rcv1(:,2),log10((abs(DLRCSGR_m_fl_rcv1(:,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%hold on
plot(SMART_CD_fl_rcv1(1:10:s3,2),log10((abs(SMART_CD_fl_rcv1(1:10:s3,3)- f_opt)/f_opt)),'--x','LineWidth',2);
hold on
plot(KATYUSHA_fl_rcv1(1:6:s4,2),log10((abs(KATYUSHA_fl_rcv1(1:6:s4,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%hold on
%plot(ADMM_fl_rcv1(:,2),ADMM_fl_rcv1(:,6),'LineWidth',5);
%hold on
%plot(CVX_fl_rcv1(:,3),log10(abs(CVX_fl_rcv1(:,1)- f_opt)),'--+','LineWidth',5);
ylim([-4 4]);
xlim([0 1800]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('rcv1');
%legend('ASGARD-fc','DLRCSGR-eps','DLRCSGR-m','SMART-CD');
legend('ASGARD-DL','IPALM-APPROX','SMART-CD','IPALM-KATYUSHA');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/fl_rcv1.eps','epsc');

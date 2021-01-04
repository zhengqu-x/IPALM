f_ub= 2.10983690e+00;
f_lb= 2.10979594e+00;
ASGARD_fc_bp_rcv2= importdata('../IPALM/results/ALM2_BP_outer_rcv2tau_47236',' ');
DLRCSGR_eps_bp_rcv2= importdata('../IPALM/results/ALM3_BP_outer_rcv2tau_1',' ');
SMART_CD_bp_rcv2= importdata('../IPALM/results/PDCD_I_BP_outer_rcv2tau_1',' ');
ADMM_bp_rcv2= importdata('../IPALM/results/ADMM_BP_rcv2_beta_1',' ');
s1= size(ASGARD_fc_bp_rcv2,1);
s2= size(DLRCSGR_eps_bp_rcv2,1);
s3= size(SMART_CD_bp_rcv2,1);
s4= size(ADMM_bp_rcv2,1);
figure(1);
plot(ASGARD_fc_bp_rcv2(1:8:s1,2),log10(abs(ASGARD_fc_bp_rcv2(1:8:s1,5)- f_ub)/f_ub),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_rcv2(1:6:s2,2),log10(abs(DLRCSGR_eps_bp_rcv2(1:6:s2,5)- f_lb)/f_lb),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_rcv2(:,2),log10(abs(DLRCSGR_m_bp_rcv2(:,5)- f_opt)/f_opt),'--o','LineWidth',2);
%hold on
plot(ADMM_bp_rcv2(1:20:s4,2),log10(abs(ADMM_bp_rcv2(1:20:s4,3)- f_ub)/f_ub),'--+','LineWidth',2);
hold on
plot(SMART_CD_bp_rcv2(1:40:s3,2),log10((abs(SMART_CD_bp_rcv2(1:40:s3,4)- f_ub)/f_ub)),'--x','LineWidth',2);
%hold on
%plot(ADMM_bp_rcv2(:,2),ADMM_bp_rcv2(:,6),'LineWidth',5);
%hold on
%plot(CVX_bp_rcv2(:,4),log10(abs(CVX_bp_rcv2(:,2)- f_opt)/f_opt),'--+','LineWidth',2);
ylim([-4 5]);
xlim([0 2200]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('rcv1');
%legend('ASGARD-fc','ASGARD-nurc');
%legend('ASGARD-fc','DLRCSGR-eps','DLRCSGR-m','SMART-CD','CVX');
legend({'ASGARD-DL','IPALM-APPROX','LADMM','SMART-CD'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(1),'myplots/bp_rcv1_obj.eps','epsc');
figure(2);
plot(ASGARD_fc_bp_rcv2(1:8:s1,2),log10(ASGARD_fc_bp_rcv2(1:8:s1,4)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_rcv2(1:6:s2,2),log10(DLRCSGR_eps_bp_rcv2(1:6:s2,4)),'--s','LineWidth',2);
hold on
plot(ADMM_bp_rcv2(1:20:s4,2),log10(ADMM_bp_rcv2(1:20:s4,4)),'--+','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_rcv2(:,2),log10(DLRCSGR_m_bp_rcv2(:,4)),'--o','LineWidth',2);
%hold on
plot(SMART_CD_bp_rcv2(1:40:s3,2),log10(SMART_CD_bp_rcv2(1:40:s3,3)),'--x','LineWidth',2);
xlim([0 2200]);
ylim([-7 4]);
xlabel('time');
ylabel('log|Ax- b|');
title('rcv1');
%legend('ASGARD-fc','ASGARD-nurc');
legend({'ASGARD-DL','IPALM-APPROX','LADMM','SMART-CD'},'interpreter','latex','Fontsize',10);
%legend({'ASGARD-fc','DLRCSGR-eps','DLRCSGR-m','CVX'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(2),'myplots/bp_rcv1_infeas.eps','epsc');

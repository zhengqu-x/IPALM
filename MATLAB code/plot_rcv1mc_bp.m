f_opt= 1.64446161;
ASGARD_fc_bp_rcv1mc= importdata('../IPALM/results/ALM2_BP_outer_rcv1mctau_47236',' ');
DLRCSGR_eps_bp_rcv1mc= importdata('../IPALM/results/ALM3_BP_outer_rcv1mctau_1',' ');
SMART_CD_bp_rcv1mc= importdata('../IPALM/results/PDCD_I_BP_outer_rcv1mctau_1',' ');
s1= size(ASGARD_fc_bp_rcv1mc,1);
s2= size(DLRCSGR_eps_bp_rcv1mc,1);
s3= size(SMART_CD_bp_rcv1mc,1);
figure(1);
plot(ASGARD_fc_bp_rcv1mc(1:10:s1,2),log10(abs(ASGARD_fc_bp_rcv1mc(1:10:s1,5)- f_opt)/f_opt),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_rcv1mc(1:10:s2,2),log10(abs(DLRCSGR_eps_bp_rcv1mc(1:10:s2,5)- f_opt)/f_opt),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_rcv1mc(:,2),log10(abs(DLRCSGR_m_bp_rcv1mc(:,5)- f_opt)/f_opt),'--o','LineWidth',2);
%hold on
plot(SMART_CD_bp_rcv1mc(1:120:s3,2),log10((abs(SMART_CD_bp_rcv1mc(1:120:s3,4)- f_opt)/f_opt)),'--x','LineWidth',2);
%hold on
%plot(ADMM_bp_rcv1mc(:,2),ADMM_bp_rcv1mc(:,6),'LineWidth',5);
%hold on
%plot(CVX_bp_rcv1mc(:,4),log10(abs(CVX_bp_rcv1mc(:,2)- f_opt)/f_opt),'--+','LineWidth',2);
ylim([-4 5]);
xlim([0 3000]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('rcv1mc');
%legend('ASGARD-fc','ASGARD-nurc');
%legend('ASGARD-fc','DLRCSGR-eps','DLRCSGR-m','SMART-CD','CVX');
legend({'ASGARD-DL','IPALM-APPROX','SMART-CD'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(1),[pwd '/my plots/bp_rcv1mc_obj.eps']);
figure(2);
plot(ASGARD_fc_bp_rcv1mc(1:10:s1,2),log10(ASGARD_fc_bp_rcv1mc(1:10:s1,4)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_rcv1mc(1:10:s2,2),log10(DLRCSGR_eps_bp_rcv1mc(1:10:s2,4)),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_rcv1mc(:,2),log10(DLRCSGR_m_bp_rcv1mc(:,4)),'--o','LineWidth',2);
%hold on
plot(SMART_CD_bp_rcv1mc(1:120:s3,2),log10(SMART_CD_bp_rcv1mc(1:120:s3,3)),'--x','LineWidth',2);
xlim([0 3000]);
ylim([-7 4]);
xlabel('time');
ylabel('log|Ax- b|');
title('rcv1mc');
%legend('ASGARD-fc','ASGARD-nurc');
legend({'ASGARD-DL','IPALM-APPROX','SMART-CD'},'interpreter','latex','Fontsize',10);
%legend({'ASGARD-fc','DLRCSGR-eps','DLRCSGR-m','CVX'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(2),[pwd '/my plots/bp_rcv1mc_infeas.eps']);

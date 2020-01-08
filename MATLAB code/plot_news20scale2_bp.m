f_opt= 481.3869;
ASGARD_fc_bp_news20scale2= importdata('../IPALM/results/ALM2_BP_outer_news20scale2tau_62061',' ');
DLRCSGR_eps_bp_news20scale2= importdata('../IPALM/results/ALM3_BP_outer_news20scale2tau_1',' ');
SMART_CD_bp_news20scale2= importdata('../IPALM/results/PDCD_I_BP_outer_news20scale2tau_1',' ');
s1= size(ASGARD_fc_bp_news20scale2,1);
s2= size(DLRCSGR_eps_bp_news20scale2,1);
s3= size(SMART_CD_bp_news20scale2,1);
figure(1);
plot(ASGARD_fc_bp_news20scale2(1:10:s1,2),log10(abs(ASGARD_fc_bp_news20scale2(1:10:s1,5)- f_opt)/f_opt),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_news20scale2(1:15:s2,2),log10(abs(100*DLRCSGR_eps_bp_news20scale2(1:15:s2,5)- f_opt)/f_opt),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_news20scale2(1:3:49,2),log10(abs(DLRCSGR_m_bp_news20scale2(1:3:49,5)- f_opt)/f_opt),'--o','LineWidth',2);
%hold on
%plot(ADMM_bp_news20scale2(:,2),ADMM_bp_news20scale2(:,6),'LineWidth',5);
%hold on
plot(SMART_CD_bp_news20scale2(1:40:s3,2),log10((abs(SMART_CD_bp_news20scale2(1:40:s3,4)- f_opt)/f_opt)),'--x','LineWidth',2);
%plot(CVX_bp_news20scale2(:,4),log10(abs(CVX_bp_news20scale2(:,2)- f_opt)/f_opt),'--+','LineWidth',2);
ylim([-6 4]);
xlim([0 1500]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('news20scale2');
%legend('ASGARD-fc','ASGARD-nurc');
legend({'ASGARD-DL','IPALM-APPROX','SMART-CD'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(1),[pwd '/my plots/bp_news20scale2_obj.eps']);
figure(2);
plot(ASGARD_fc_bp_news20scale2(1:10:s1,2),log10(ASGARD_fc_bp_news20scale2(1:10:s1,4)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_news20scale2(1:15:s2,2),log10(DLRCSGR_eps_bp_news20scale2(1:15:s2,4)),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_news20scale2(1:3:49,2),log10(DLRCSGR_m_bp_news20scale2(1:3:49,4)),'--o','LineWidth',2);
plot(SMART_CD_bp_news20scale2(1:40:s3,2),log10(SMART_CD_bp_news20scale2(1:40:s3,3)),'--x','LineWidth',2);
xlim([0 1500]);
ylim([-7 4]);
xlabel('time');
ylabel('log|Ax- b|');
title('news20scale2');
%legend('ASGARD-fc','ASGARD-nurc');
legend({'ASGARD-DL','IPALM-APPROX','SMART-CD'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(2),[pwd '/my plots/bp_news20scale2_infeas.eps']);

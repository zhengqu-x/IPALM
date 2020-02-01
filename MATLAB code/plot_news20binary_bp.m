f_ub= 20.6456;
f_lb= 20.64553;
ASGARD_fc_bp_news20binary= importdata('../IPALM/results/ALM2_BP_outer_news20binary_bptau_1355191',' ');
DLRCSGR_eps_bp_news20binary= importdata('../IPALM/results/ALM3_BP_outer_news20binary_bptau_1',' ');
SMART_CD_bp_news20binary= importdata('../IPALM/results/PDCD_I_BP_outer_news20binary_bptau_1',' ');
s1= size(ASGARD_fc_bp_news20binary,1);
s2= size(DLRCSGR_eps_bp_news20binary,1);
s3= size(SMART_CD_bp_news20binary,1);
figure(1);
plot(ASGARD_fc_bp_news20binary(1:7:s1,2),log10(abs(ASGARD_fc_bp_news20binary(1:7:s1,5)- f_ub)/f_ub),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_news20binary(1:4:s2,2),log10(abs(DLRCSGR_eps_bp_news20binary(1:4:s2,5)- f_lb)/f_lb),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_news20binary(:,2),log10(abs(DLRCSGR_m_bp_news20binary(:,5)- f_opt)/f_opt),'--o','LineWidth',2);
%hold on
plot(SMART_CD_bp_news20binary(1:2:s3,2),log10((abs(SMART_CD_bp_news20binary(1:2:s3,4)- f_ub)/f_ub)),'--x','LineWidth',2);
%hold on
%plot(ADMM_bp_news20binary(:,2),ADMM_bp_news20binary(:,6),'LineWidth',5);
%hold on
%plot(CVX_bp_news20binary(:,4),log10(abs(CVX_bp_news20binary(:,2)- f_opt)/f_opt),'--+','LineWidth',2);
ylim([-2 2]);
xlim([0 2200]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('news20binary');
%legend('ASGARD-fc','ASGARD-nurc');
%legend('ASGARD-fc','DLRCSGR-eps','DLRCSGR-m','SMART-CD','CVX');
legend({'ASGARD-DL','IPALM-APPROX','SMART-CD'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(1),'myplots/bp_news20binary_obj.eps','epsc');
figure(2);
plot(ASGARD_fc_bp_news20binary(1:7:s1,2),log10(ASGARD_fc_bp_news20binary(1:7:s1,4)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_bp_news20binary(1:4:s2,2),log10(DLRCSGR_eps_bp_news20binary(1:4:s2,4)),'--s','LineWidth',2);
hold on
%plot(DLRCSGR_m_bp_news20binary(:,2),log10(DLRCSGR_m_bp_news20binary(:,4)),'--o','LineWidth',2);
%hold on
plot(SMART_CD_bp_news20binary(1:2:s3,2),log10(SMART_CD_bp_news20binary(1:2:s3,3)),'--x','LineWidth',2);
xlim([0 2200]);
ylim([-3 2]);
xlabel('time');
ylabel('log|Ax- b|');
title('news20binary');
%legend('ASGARD-fc','ASGARD-nurc');
legend({'ASGARD-DL','IPALM-APPROX','SMART-CD'},'interpreter','latex','Fontsize',10);
%legend({'ASGARD-fc','DLRCSGR-eps','DLRCSGR-m','CVX'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 350 400]);
saveas(figure(2),'myplots/bp_news20binary_infeas.eps','epsc');

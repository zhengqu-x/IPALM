f_opt= 292.858;
ASGARD_fc_l1l_news20binary= importdata('../IPALM/results/DLRCSGR2_LAD_outer_news20binarytau_1355191',' ');
DLRCSGR_eps_l1l_news20binary= importdata('../IPALM/results/DLRCSGR3_LAD_outer_news20binarytau_1',' ');
SMART_CD_l1l_news20binary= importdata('../IPALM/results/PDCD_LAD_outer_news20binarytau_1',' ');
ADMM_l1l_news20binary= importdata('../IPALM/results/ADMM_LAD_news20binary_beta_1',' ');
s1= size(ASGARD_fc_l1l_news20binary,1);
s2= size(DLRCSGR_eps_l1l_news20binary,1);
s3= size(SMART_CD_l1l_news20binary,1);
s4= size(ADMM_l1l_news20binary);
plot(ASGARD_fc_l1l_news20binary(:,2),log10((abs(ASGARD_fc_l1l_news20binary(:,3)- f_opt)/f_opt)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_l1l_news20binary(1:2:s2,2),log10((abs(DLRCSGR_eps_l1l_news20binary(1:2:s2,3)- f_opt)/f_opt)),'--s','LineWidth',2);
hold on
plot(ADMM_l1l_news20binary(1:2:s4,2),log10(abs(ADMM_l1l_news20binary(1:2:s4,3)- f_opt)/f_opt),'--+','LineWidth',2);
hold on
plot(SMART_CD_l1l_news20binary(:,2),log10((abs(SMART_CD_l1l_news20binary(:,3)- f_opt)/f_opt)),'--x','LineWidth',2);
%plot(DLRCSGR_m_l1l_news20binary(1:2:52,2),log10((abs(DLRCSGR_m_l1l_news20binary(1:2:52,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%hold on
%plot(ADMM_l1l_news20binary(:,2),ADMM_l1l_news20binary(:,6),'LineWidth',5);
%hold on
%plot(CVX_l1l_news20binary(:,3),log10((abs(CVX_l1l_news20binary(:,1)- f_opt)/f_opt)),'--+','LineWidth',2);
ylim([-2 4]);
xlim([0 3000]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('news20binary');
legend({'ASGARD-DL','IPALM-APPROX','LADMM','SMART-CD'},'interpreter','latex','Fontsize',10);
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/lasso_news20binary.eps','epsc');

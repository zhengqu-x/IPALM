f_opt= 4.81386755;
ASGARD_fc_l1l_news20scale2= importdata('../IPALM/results/DLRCSGR2_LAD_outer_news20scale2tau_62061',' ');
DLRCSGR_eps_l1l_news20scale2= importdata('../IPALM/results/DLRCSGR3_LAD_outer_news20scale2tau_1',' ');
SMART_CD_l1l_news20scale2= importdata('../IPALM/results/PDCD_LAD_outer_news20scale2tau_1',' ');
s1= size(ASGARD_fc_l1l_news20scale2,1);
s2= size(DLRCSGR_eps_l1l_news20scale2,1);
s3= size(SMART_CD_l1l_news20scale2,1);
plot(ASGARD_fc_l1l_news20scale2(1:2:s1,2),log10((abs(ASGARD_fc_l1l_news20scale2(1:2:s1,3)- f_opt)/f_opt)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_l1l_news20scale2(1:20:s2,2),log10((abs(DLRCSGR_eps_l1l_news20scale2(1:20:s2,3)- f_opt)/f_opt)),'--s','LineWidth',2);
hold on
plot(SMART_CD_l1l_news20scale2(1:40:s3,2),log10((abs(SMART_CD_l1l_news20scale2(1:40:s3,3)- f_opt)/f_opt)),'--x','LineWidth',2);
%plot(DLRCSGR_m_l1l_news20scale2(1:2:52,2),log10((abs(DLRCSGR_m_l1l_news20scale2(1:2:52,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%hold on
%plot(ADMM_l1l_news20scale2(:,2),ADMM_l1l_news20scale2(:,6),'LineWidth',5);
%hold on
%plot(CVX_l1l_news20scale2(:,3),log10((abs(CVX_l1l_news20scale2(:,1)- f_opt)/f_opt)),'--+','LineWidth',2);
ylim([-5 5]);
xlim([0 1000]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('news20scale2');
legend('ASGARD-DL','IPALM-APPROX','SMART-CD');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/lasso_news20scale2.eps','epsc');

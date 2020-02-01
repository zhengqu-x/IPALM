f_opt= 7.372570420000000e+03;
ASGARD_fc_svm_realsim= importdata('../IPALM/results/DLRCSGR2_SVM_outer_realsimtau_20959',' ');
DLRCSGR_eps_svm_realsim= importdata('../IPALM/results/DLRCSGR3_SVM_outer_realsimtau_1',' ');
SMART_CD_svm_realsim= importdata('../IPALM/results/PDCD_SMSVM_SVM_outer_realsimtau_1',' ');
KATYUSHA_svm_realsim= importdata('../IPALM/results/ALM_SVM_outer_realsimtau_268',' ');
s1= size(ASGARD_fc_svm_realsim,1);
s2= size(DLRCSGR_eps_svm_realsim,1);
s3= size(SMART_CD_svm_realsim,1);
s4= size(KATYUSHA_svm_realsim,1);
plot(ASGARD_fc_svm_realsim(1:4:s1,2),log10((abs(ASGARD_fc_svm_realsim(1:4:s1,3)- f_opt)/f_opt)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_svm_realsim(1:9:s2,2),log10((abs(DLRCSGR_eps_svm_realsim(1:9:s2,3)- f_opt)/f_opt)),'--s','LineWidth',2);
hold on
plot(SMART_CD_svm_realsim(1:2:s3,2),log10((abs(SMART_CD_svm_realsim(1:2:s3,5)- f_opt)/f_opt)),'--x','LineWidth',2);
%plot(DLRCSGR_m_svm_realsim(1:2:52,2),log10((abs(DLRCSGR_m_svm_realsim(1:2:52,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%hold on
plot(KATYUSHA_svm_realsim(1:7:s4,2),log10((abs(KATYUSHA_svm_realsim(1:7:s4,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%plot(ADMM_svm_realsim(:,2),ADMM_svm_realsim(:,6),'LineWidth',5);
%hold on
%plot(CVX_svm_realsim(:,3),log10((abs(CVX_svm_realsim(:,1)- f_opt)/f_opt)),'--+','LineWidth',2);
ylim([-4 4]);
xlim([0 3600]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('realsim');
legend('ASGARD-DL','IPALM-APPROX','SMART-CD','IPALM-KATYUSHA');
%legend('ASGARD-DL','IPALM-APPROX','SMART-CD');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/svm_realsim.eps','epsc');

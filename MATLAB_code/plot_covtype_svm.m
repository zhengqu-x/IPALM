f_opt= 3.379019290232700e+05;
f_lb= 3.376823216715768e+05;
ASGARD_fc_svm_covtypem= importdata('../IPALM/results/DLRCSGR2_SVM_outer_covtypetau_55',' ');
DLRCSGR_eps_svm_covtypem= importdata('../IPALM/results/DLRCSGR3_SVM_outer_covtypetau_1',' ');
SMART_CD_svm_covtypem= importdata('../IPALM/results/PDCD_SMSVM_SVM_outer_covtypetau_1',' ');
KATYUSHA_svm_covtypem= importdata('../IPALM/results/ALM_SVM_outer_covtypetau_762',' ');
ADMM_svm_covtypem= importdata('../IPALM/results/ADMM_SVM_covtype_beta_2',' ');
s1= size(ASGARD_fc_svm_covtypem,1);
s2= size(DLRCSGR_eps_svm_covtypem,1);
s3= size(SMART_CD_svm_covtypem,1);
s4= size(KATYUSHA_svm_covtypem,1);
s5= size(ADMM_svm_covtypem,1);
plot(ASGARD_fc_svm_covtypem(1:2:s1,2),log10((abs(ASGARD_fc_svm_covtypem(1:2:s1,3)- f_opt)/f_opt)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_svm_covtypem(1:s2,2),log10((abs(DLRCSGR_eps_svm_covtypem(1:s2,3)- f_opt)/f_opt)),'--s','LineWidth',2);
hold on
plot(ADMM_svm_covtypem(1:2:s5,2),log10(abs(ADMM_svm_covtypem(1:2:s5,3)- f_opt)/f_opt),'--+','LineWidth',2);
hold on
plot(SMART_CD_svm_covtypem(1:s3,2),log10((abs(SMART_CD_svm_covtypem(1:s3,5)- f_opt)/f_opt)),'--x','LineWidth',2);
%plot(DLRCSGR_m_svm_covtypem(1:2:52,2),log10((abs(DLRCSGR_m_svm_covtypem(1:2:52,3)- f_opt)/f_opt)),'--o','LineWidth',2);
hold on
plot(KATYUSHA_svm_covtypem(1:s4,2),log10((abs(KATYUSHA_svm_covtypem(1:s4,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%plot(ADMM_svm_covtypem(:,2),ADMM_svm_covtypem(:,6),'LineWidth',5);
%hold on
%plot(CVX_svm_covtypem(:,3),log10((abs(CVX_svm_covtypem(:,1)- f_opt)/f_opt)),'--+','LineWidth',2);
ylim([-4.5 2]);
xlim([0 1000]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('covtype');
legend('ASGARD-DL','IPALM-APPROX','LADMM','SMART-CD','IPALM-KATYUSHA');
%legend('ASGARD-DL','IPALM-APPROX','SMART-CD');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/svm_covtype.eps','epsc');

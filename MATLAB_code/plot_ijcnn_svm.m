f_opt= 8.568546470949960e+3;
f_lb= 8.568544079948495e+03;
ASGARD_fc_svm_ijcnn1m= importdata('../IPALM/results/DLRCSGR2_SVM_outer_ijcnn1tau_23',' ');
DLRCSGR_eps_svm_ijcnn1m= importdata('../IPALM/results/DLRCSGR3_SVM_outer_ijcnn1tau_1',' ');
SMART_CD_svm_ijcnn1m= importdata('../IPALM/results/PDCD_SMSVM_SVM_outer_ijcnn1tau_1',' ');
KATYUSHA_svm_ijcnn1m= importdata('../IPALM/results/ALM_SVM_outer_ijcnn1tau_223',' ');
ADMM_svm_ijcnn1m= importdata('../IPALM/results/ADMM_SVM_ijcnn1_beta_0',' ');
s1= size(ASGARD_fc_svm_ijcnn1m,1);
s2= size(DLRCSGR_eps_svm_ijcnn1m,1);
s3= size(SMART_CD_svm_ijcnn1m,1);
s4= size(KATYUSHA_svm_ijcnn1m,1);
s5= size(ADMM_svm_ijcnn1m,1);
plot(ASGARD_fc_svm_ijcnn1m(1:5:s1,2),log10((abs(ASGARD_fc_svm_ijcnn1m(1:5:s1,3)- f_opt)/f_opt)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_svm_ijcnn1m(1:10:s2,2),log10((abs(DLRCSGR_eps_svm_ijcnn1m(1:10:s2,3)- f_opt)/f_opt)),'--s','LineWidth',2);
hold on
plot(ADMM_svm_ijcnn1m(1:30:s5,2),log10(abs(ADMM_svm_ijcnn1m(1:30:s5,3)- f_opt)/f_opt),'--+','LineWidth',2);
hold on
plot(SMART_CD_svm_ijcnn1m(1:20:s3,2),log10((abs(SMART_CD_svm_ijcnn1m(1:20:s3,5)- f_opt)/f_opt)),'--x','LineWidth',2);
%plot(DLRCSGR_m_svm_ijcnn1m(1:2:52,2),log10((abs(DLRCSGR_m_svm_ijcnn1m(1:2:52,3)- f_opt)/f_opt)),'--o','LineWidth',2);
hold on
plot(KATYUSHA_svm_ijcnn1m(1:10:s4,2),log10((abs(KATYUSHA_svm_ijcnn1m(1:10:s4,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%plot(ADMM_svm_ijcnn1m(:,2),ADMM_svm_ijcnn1m(:,6),'LineWidth',5);
%hold on
%plot(CVX_svm_ijcnn1m(:,3),log10((abs(CVX_svm_ijcnn1m(:,1)- f_opt)/f_opt)),'--+','LineWidth',2);
ylim([-7 2]);
xlim([0 1000]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('ijcnn1');
legend('ASGARD-DL','IPALM-APPROX','LADMM','SMART-CD','IPALM-KATYUSHA');
%legend('ASGARD-DL','IPALM-APPROX','SMART-CD');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/svm_ijcnn1.eps','epsc');

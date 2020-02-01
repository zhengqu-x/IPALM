f_opt= 5.7252151300236e+03;
f_lb= 5.725203673581847e+03;
ASGARD_fc_svm_a7am= importdata('../IPALM/results/DLRCSGR2_SVM_outer_a7atau_123',' ');
DLRCSGR_eps_svm_a7am= importdata('../IPALM/results/DLRCSGR3_SVM_outer_a7atau_1',' ');
SMART_CD_svm_a7am= importdata('../IPALM/results/PDCD_SMSVM_SVM_outer_a7atau_1',' ');
KATYUSHA_svm_a7am= importdata('../IPALM/results/ALM_SVM_outer_a7atau_126',' ');
s1= size(ASGARD_fc_svm_a7am,1);
s2= size(DLRCSGR_eps_svm_a7am,1);
s3= size(SMART_CD_svm_a7am,1);
s4= size(KATYUSHA_svm_a7am,1);
plot(ASGARD_fc_svm_a7am(1:5:s1,2),log10((abs(ASGARD_fc_svm_a7am(1:5:s1,3)- f_opt)/f_opt)),'-->','LineWidth',2);
hold on
plot(DLRCSGR_eps_svm_a7am(1:20:s2,2),log10((abs(DLRCSGR_eps_svm_a7am(1:20:s2,3)- f_opt)/f_opt)),'--s','LineWidth',2);
%hold on
plot(SMART_CD_svm_a7am(1:30:s3,2),log10((abs(SMART_CD_svm_a7am(1:30:s3,5)- f_opt)/f_opt)),'--x','LineWidth',2);
%plot(DLRCSGR_m_svm_a7am(1:2:52,2),log10((abs(DLRCSGR_m_svm_a7am(1:2:52,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%hold on
plot(KATYUSHA_svm_a7am(1:10:s4,2),log10((abs(KATYUSHA_svm_a7am(1:10:s4,3)- f_opt)/f_opt)),'--o','LineWidth',2);
%plot(ADMM_svm_a7am(:,2),ADMM_svm_a7am(:,6),'LineWidth',5);
%hold on
%plot(CVX_svm_a7am(:,3),log10((abs(CVX_svm_a7am(:,1)- f_opt)/f_opt)),'--+','LineWidth',2);
ylim([-6 2]);
xlim([0 1000]);
xlabel('time');
ylabel('log|F(x)- F^*|/F^*');
title('a7a');
legend('ASGARD-DL','IPALM-APPROX','SMART-CD','IPALM-KATYUSHA');
%legend('ASGARD-DL','IPALM-APPROX','SMART-CD');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/svm_a7a.eps','epsc');

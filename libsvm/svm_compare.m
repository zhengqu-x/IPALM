[b, A] = libsvmread('w8a.txt');
[b2, A2]= libsvmread('w8a.t');
tic
model = svmtrain(b, A );
toc
[pl, acc, dop] = svmpredict(b2, A2, model );
%acc is the accuracy of libsvm

n1= size(A,2);
n2= size(A2,2);
%copy the output x from the folder "/IPALM/results"
%x= importdata("L_Katyusha_x_SVM_outer_a7atau_126");
%x= importdata("L_Katyusha_x_SVM_outer_a8atau_150");
%x= importdata("L_Katyusha_x_SVM_outer_a9atau_180");
%x= importdata("L_Katyusha_x_SVM_outer_ijcnn1tau_223");
%x= importdata("L_Katyusha_x_SVM_outer_w6atau_131");
%x= importdata("L_Katyusha_x_SVM_outer_w7atau_157");
x= importdata("L_Katyusha_x_SVM_outer_w8atau_223");
n= min(n1,n2);
b_test= A2(:,1:n)*x(1:n,1)+ x(n1+ 1,1);
com= b2.*b_test;
accuracy= size(com(com>=0),1)/size(A2,1);
%accuracy is the accuracy of IPALM

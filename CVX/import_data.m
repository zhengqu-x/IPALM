[b, A]= libsvmread('matrix_news20binary');
save('news20binary_b.mat','b');
save('news20binary_A.mat','A');
[b1, A1]= libsvmread('matrix_news20binary_bp');
save('news20binary_bp_b.mat','b1');
save('news20binary_bp_A.mat','A1');

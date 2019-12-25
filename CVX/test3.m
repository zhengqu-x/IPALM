A = importdata('news20binary_A.mat'); b = importdata('news20binary_b.mat'); 
%M= importdata('Fused_news20scale2.mat');
n= size(A,2); 
m= size(A,1);
cvx_begin
    variable x(n)
    minimize(  0.01*norm(x,1)+ norm(A*x- b,1) ) 
cvx_end

A = importdata('news20binary_A.mat'); b = importdata('news20binary_b.mat'); 
M= importdata('Fused_news20binary.mat');
n= size(A,2); 
m= size(A,1);
cvx_begin
    variable x(n)
    minimize(  0.5*(A*x- b)'*(A*x- b)+ 0.01*norm(x,1)+ 0.01*norm(M*x,1) ) 
cvx_end

A = importdata('news20binary_A.mat'); b = importdata('news20binary_bp_b.mat'); 
%M= importdata('Fused_leu.mat');
n= size(A,2);
cvx_begin
    variable x(n)
    minimize(  0.01*norm(x,1) )
    subject to
                  A*x== b 
cvx_end
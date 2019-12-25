A = importdata('rcv1_A.mat'); b = importdata('rcv1_b.mat'); 
%M= importdata('Fused_leu.mat');
n= size(A,2);
cvx_begin
    variable x(n)
    minimize(  0.01*norm(x,1) )
    subject to
                  A*x== b 
cvx_end

A = importdata('qcqp1_A.mat'); 
%A1 = importdata('qcqp4_A1.mat'); 
%A2 = importdata('qcqp4_A2.mat'); 
%A3 = importdata('qcqp4_A3.mat'); 
%A4 = importdata('qcqp4_A4.mat'); 
%A5 = importdata('qcqp4_A5.mat'); 
%A= [A1;A2;A3;A4;A5];
b = importdata('qcqp1_b.mat');
m = size(A,1)/size(A,2)-1 ;
n = size(A,2); 
e= ones(n,1);
ne= -1*ones(n,1);
tic
cvx_begin
    variable x(n)
    minimize( 0.5*x'*A(1:n,:)*x+ b(1:n,1)'*x )
    subject to
        for i= 1:m
            0.5*x'*A(n*i+1:n*(i+1),:)*x+ b(n*i+1:n*(i+ 1),1)'*x<= 1
        end
        x<= e
        x>= ne
cvx_end
toc

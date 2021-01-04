n= 1000;
% the number of the constraints
m= 500;
A= sparse(n*(m+ 1),n);
b= zeros(n*(m+ 1),1);
for i= 1:m+1
    R = sprandsym(n,0.1,0.1);
    R = R'*R;
    c= randn(n,1);
    A((i- 1)*n+1:i*n,:)= R;
    b((i- 1)*n+1:i*n,:)= c;
end
% If m is too large, matlab cannot save the matrix A. We need seperate the
% matrix A into 5 parts.
A1= sparse(n*m/5,n);
A1= A(1:n*(m/5),:);
A2= sparse(n*m/5,n);
A2= A(n*(m/5)+ 1:2*n*(m/5),:);
A3= sparse(n*m/5,n);
A3= A(2*n*m/5+ 1:3*n*(m/5),:);
A4= sparse(n*m/5,n);
A4= A(3*n*m/5+ 1:4*n*(m/5),:);
A5= sparse(n*(m/5+ 1),n);
A5= A(4*n*m/5+ 1:n*(m+ 1),:);
save qcqp4_A1.mat A1;
save qcqp4_A2.mat A2;
save qcqp4_A3.mat A3;
save qcqp4_A4.mat A4;
save qcqp4_A5.mat A5;
save qcqp4_b.mat b;
% For small m, we can directlt save.
% save qcqp3_A.mat A;
% save qcqp3_b.mat b;
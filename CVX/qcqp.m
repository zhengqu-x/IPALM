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
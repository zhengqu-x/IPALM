Instruction for CVX

This folder contains the code and data used to solving the optimization problem by CVX. So readers need download and install CVX package (http://cvxr.com/cvx/).

We include in the folder 3 datasets news20scale2, rcv1, rcv1mc. 
Since news20binary and news20bianry_bp are too large to upload. Readers need download and install LIBSVM package (https://www.csie.ntu.edu.tw/~cjlin/libsvm/#download). Then run import_data.m.




To solve basis pursuit problem, 
Run test2.m
For different dataset, just change the dataset name to define A and b in test2.m

To solve least absolute deviation problem, 
Run test3.m
For different dataset, just change the dataset name to define A and b in test3.m

To solve fused lasso problem, 
Run test4.m
For different dataset, just change the dataset name to define A, b and M in test4.m

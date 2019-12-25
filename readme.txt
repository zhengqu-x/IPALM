General instruction

This small package contains a code which was used to implement algorithms and generate plots in the paper
F. Li and Z. Qu: An inexact proximal augmented Lagrangian framework with arbitrary linearly convergent inner solver for composite convex optimization.


The optimization problem is:
minimize {f(x)+ g(x)+ h(Mx)}
with:
1. f(x) is convex and differentiable on an open set containing dom(g)
2. g(x) is proper, convex and closed
3. h(x) is proper, convex and closed or an indicator function of a convex and closed set 
4. M is a linear operator (matrix)

Problem type:
1. Basis pursuit: f(x)= 0, g(x)= lambda1*||x||_1, h(x)= I_{x= b}
2. Least absolute deviation: f(x)= 0, g(x)= lambda1*||x||_1, h(x)= lambda3*||x||_1
3. Fused lasso: f(x)= 0.5*||Ax- b||_2^2, g(x)= lambda1*||x||_1, h(x)= lambda3*||x||_1
4. L_1 norm support vector machine: f(x)= 0, g(x,w)= lambda1*||x||_1, h(x)= max(0,1- x)

Implemented algorithms:
a. IPALM_APPROX
b. IPALM_Katyusha
c. SMART_CD
d. ASGARD_DL



The matrices should be entered either in LIBSVM format. 

We include in the folder datas/ 16 datasets news20scale2, rcv1, rcv2, rcv1mc, news20binary, news20binary_bp, ijcnn1, w4a, w6a, w7a, w8a, a7a, a8a, a9a, real-sim, covtype, downloaded or modified from https://www.csie.ntu.edu.tw/~cjlin/libsvm/, that are used in our experiments. Since real-sim is too large to upload, readers need to download it by themselves and save it as matrix_realsim. Datasets news20binary, news20binary_bp and covtype are zip files. So readers need unzip them before running the following commands.




To compile, please do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas

To run, please type:
>> ./main arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9 arg10

arg1: problem_type;
arg2: algorithm;
arg3: filename;
arg4: lambda1;
arg5: lambda3;
arg6: maximum number of outer iterations;
arg7: maximum running time;
arg8: (optional): beta_0, default value= 1; 
arg9: (optional): epsilon_0, default value=0;
arg10: (optional): p_N_2 (print or compute every p_N_2 iterations), default value=100;

For example, to solve the least absolute deviation problem with lambda1=0.01, lambda3= 1 in the paper, using IPALM_APPROX, beta_0= 1, epsilon_0= 100, p_N_2= 100, do: 

>>main 2 a news20scale2 0.01 1 100000 100 1 100 100 

If instead you want to solve the problem with SMART_CD, and beta_0= 1000, p_N_2= 50, do:

>>main 2 c news20scale2 0.01 1 100000 100 1000 50  




The results are saved in the folder results/.




This code comes as is, with no guarantee. We did not perform extensive testings. We welcome any suggestion to solve bugs or improve the code.

Reproduction of figures in the paper

To reproduce Figure 1(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a news20scale2 0.01 1 100000 1200 1 100 100 
>> ./main 2 c news20scale2 0.01 1 100000 1200 1000 50 
>> ./main 2 d news20scale2 0.01 1 100000 1200 0.01 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_news20scale2_lasso



To reproduce Figure 1(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a rcv1 0.01 1 100000 1200 1 100 100 
>> ./main 2 c rcv1 0.01 1 100000 1200 1000 50 
>> ./main 2 d rcv1 0.01 1 100000 1200 0.01 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_rcv1_lasso



To reproduce Figure 1(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a rcv1mc 0.01 1 100000 1200 1 100 100 
>> ./main 2 c rcv1mc 0.01 1 100000 1200 1000 50 
>> ./main 2 d rcv1mc 0.01 1 100000 1200 0.01 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_rcv1mc_lasso



To reproduce Figure 1(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a news20binary 0.01 1 100000 3600 0.1 100 50 
>> ./main 2 c news20binary 0.01 1 100000 3600 10 10 
>> ./main 2 d news20binary 0.01 1 100000 3600 0.1 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_news20binary_lasso



To reproduce Figure 2(a) and 3(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a news20scale2 0.01 1 100000 1800 0.01 1 50 
>> ./main 1 c news20scale2 1 1 100000 1800 1 50 
>> ./main 1 d news20scale2 1 1 100000 1800 1 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_news20scale2_bp



To reproduce Figure 2(b) and 3(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a rcv2 0.01 1 100000 2400 0.01 1 50 
>> ./main 1 c rcv2 0.01 1 100000 2400 1 50 
>> ./main 1 d rcv2 0.01 1 100000 2400 1 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_rcv1_bp



To reproduce Figure 2(c) and 3(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a rcv1mc 0.01 1 100000 3600 0.01 1 50 
>> ./main 1 c rcv1mc 0.01 1 100000 3600 1 50 
>> ./main 1 d rcv1mc 0.01 1 100000 3600 1 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_rcv1mc_bp



To reproduce Figure 2(d) and 3(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a news20binary_bp 0.01 1 100000 1200 0.01 1 50 
>> ./main 1 c news20binary_bp 0.01 1 100000 1200 0.1 50 
>> ./main 1 d news20binary_bp 0.01 1 100000 1200 0.1 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_news20binary_bp



To reproduce Figure 4(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a news20scale2 0.01 0.01 100000 2400 1 100 100 
>> ./main 3 b news20scale2 0.01 0.01 100000 2400 1 1 40 
>> ./main 3 c news20scale2 0.01 0.01 100000 2400 10 50 
>> ./main 3 d news20scale2 0.01 0.01 100000 2400 10 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_news20scale2_fl



To reproduce Figure 4(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a rcv1 0.01 0.01 100000 2400 1 100 100 
>> ./main 3 b rcv1 0.01 0.01 100000 2400 1 1 40 
>> ./main 3 c rcv1 0.01 0.01 100000 2400 10 50 
>> ./main 3 d rcv1 0.01 0.01 100000 2400 10 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_rcv1_fl




To reproduce Figure 4(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a rcv1mc 0.01 0.01 100000 2400 1 100 100 
>> ./main 3 b rcv1mc 0.01 0.01 100000 2400 1 1 40 
>> ./main 3 c rcv1mc 0.01 0.01 100000 2400 10 50 
>> ./main 3 d rcv1mc 0.01 0.01 100000 2400 10 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_rcv1mc_fl




To reproduce Figure 4(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a news20binary 0.01 0.01 100000 2400 1 100 100 
>> ./main 3 b news20binary 0.01 0.01 100000 2400 1 10 20 
>> ./main 3 c news20binary 0.01 0.01 100000 2400 10 10 
>> ./main 3 d news20binary 0.01 0.01 100000 2400 0.1 


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_news20binary_fl




To reproduce Figure 5(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a w4a 1 1 100000 200 1 1000 100 
>> ./main 4 b w4a 1 1 100000 200 10 100 20 
>> ./main 4 c w4a 1 1 100000 200 100 100 
>> ./main 4 d w4a 1 1 100000 200 1


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_w4a_svm



To reproduce Figure 5(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a a7a 1 1 100000 1200 1 1000 100 
>> ./main 4 b a7a 1 1 100000 1200 1 100 20 
>> ./main 4 c a7a 1 1 100000 1200 100 100 
>> ./main 4 d a7a 1 1 100000 1200 10


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_a7a_svm




To reproduce Figure 5(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a a8a 1 1 100000 1200 1 1000 100 
>> ./main 4 b a8a 1 1 100000 1200 1 100 20 
>> ./main 4 c a8a 1 1 100000 1200 100 100 
>> ./main 4 d a8a 1 1 100000 1200 10


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_a8a_svm




To reproduce Figure 5(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./ 4 a a9a 1 1 100000 1200 1 1000 100 
>> ./main 4 b a9a 1 1 100000 1200 1 100 20 
>> ./main 4 c a9a 1 1 100000 1200 100 100 
>> ./main 4 d a9a 1 1 100000 1200 10


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_a9a_svm





To reproduce Figure 5(e) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a w6a 1 1 100000 1200 1 1000 100 
>> ./main 4 b w6a 1 1 100000 1200 1 100 20 
>> ./main 4 c w6a 1 1 100000 1200 100 100 
>> ./main 4 d w6a 1 1 100000 1200 10


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_w6a_svm




To reproduce Figure 5(f) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a w7a 1 1 100000 1200 1 1000 100 
>> ./main 4 b w7a 1 1 100000 1200 1 100 20 
>> ./main 4 c w7a 1 1 100000 1200 100 100 
>> ./main 4 d w7a 1 1 100000 1200 10


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_w7a_svm




To reproduce Figure 5(g) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a ijcnn1 1 1 100000 1200 1 1000 100 
>> ./main 4 b ijcnn1 1 1 100000 1200 1 100 20 
>> ./main 4 c ijcnn1 1 1 100000 1200 100 100 
>> ./main 4 d ijcnn1 1 1 100000 1200 10


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_ijcnn_svm





To reproduce Figure 5(h) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a realsim 1 1 100000 3600 1 100 100 
>> ./main 4 b realsim 1 1 100000 3600 1 100 20 
>> ./main 4 c realsim 1 1 100000 3600 1 100 
>> ./main 4 d realsim 1 1 100000 3600 1


Then run
>>cd ..
>>cd MATLAB code
>>Matlab plot_realsim_svm


All the plots are saved in .eps format in the folder Matlab code/my plots/

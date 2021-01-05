Instruction for LIBSVM

This folder contains the code and data used to solve support vector machine problems by LIBSVM. Readers need download and install LIBSVM package (https://www.csie.ntu.edu.tw/~cjlin/libsvm/#download).

To generate Table 7 in our paper, 
First, 
copy the output x from the folder "../IPALM/results/" to current folder:
L_Katyusha_x_SVM_outer_a7atau_126;
L_Katyusha_x_SVM_outer_a8atau_150;
L_Katyusha_x_SVM_outer_a9atau_180;
L_Katyusha_x_SVM_outer_ijcnn1tau_223;
L_Katyusha_x_SVM_outer_w6atau_131;
L_Katyusha_x_SVM_outer_w7atau_157;
L_Katyusha_x_SVM_outer_w8atau_223.

Second, 
Run svm_compare.m
For different datasets, just change the dataset name to define [b, A], [b2, A2] in the 1st and 2nd line, and x in the 18th line.

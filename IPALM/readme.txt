=========
APPROX2.h
=========

This class implements the accelerated coordinate descent method for solving composite convex minoration problem, presented in the paper <<Accelerated, parallel and proximal coordinate descent, Fercoq & Richtarik, SIAM Review, 2015.>>





=========
Primal_Dual_LOOPLESS_Katyusha0.h
=========

This class implements the method loopless Katyusha with arbitrary sampling, in the paper <<L-SVRG and L-Katyusha with arbitrary sampling, Qian, Qu and Richtarik, arXiv:1906.01481.>>




=========
ALM_APG.h
=========


This class implements the method ASGARD_DL, in the paper 
<<An adaptive primal-dual framework for nonsmooth convex minimization, by Tran-Dinh, Alacaoglu, Fercoq and Cevher, Mathematical Programming Computation, 2019.>>





=========
ALM_APPROX.h
=========


This class implements the method IPALM_APPROX, in our paper.




=========
ALM_Katyusha.h
=========


This class implements the method IPALM_Katyusha, in our paper.



=========
ALM_I_APG.h
=========


This class implements the method ASGARD_DL, in the paper 
<<An adaptive primal-dual framework for nonsmooth convex minimization, by Tran-Dinh, Alacaoglu, Fercoq and Cevher, Mathematical Programming Computation, 2019.>>, for problem with equality constraints.



=========
ALM_I_APPROX.h
=========


This class implements the method IPALM_APPROX, in our paper, for problem with equality constraints.




=========
PDCD.h
=========

This class implements the method SMART_CD, in the paper << Smooth primal-dual coordinate de- scent algorithms for nonsmooth convex optimization, by Alacaoglu, Tran-Dinh, Fercoq and Cevher, Advances in Neural Information Processing Systems, 2017.>>




=========
PDCD_I.h
=========

This class implements the method SMART_CD, in the paper << Smooth primal-dual coordinate de- scent algorithms for nonsmooth convex optimization, by Alacaoglu, Tran-Dinh, Fercoq and Cevher, Advances in Neural Information Processing Systems, 2017.>>, for problem with equality constraints.




=========
PDCD_SMSVM.h
=========

This class implements the method SMART_CD, in the paper << Smooth primal-dual coordinate de- scent algorithms for nonsmooth convex optimization, by Alacaoglu, Tran-Dinh, Fercoq and Cevher, Advances in Neural Information Processing Systems, 2017.>>, for solving soft margin SVM problem.




=========
Basis_pursuit_a.h
=========

This class inherits the class ALM_I_APPROX. It solves basis pursuit problem by the method IPALM_APPROX.



=========
Basis_pursuit_c.h
=========

This class inherits the class PDCD_I. It solves basis pursuit problem by the method SMART_CD.



=========
Basis_pursuit_d.h
=========

This class inherits the class ALM_I_APG. It solves basis pursuit problem by the method ASGARD_DL.



=========
Fused_lasso_a.h
=========

This class inherits the class ALM_APPROX. It solves fused lasso problem by the method IPALM_APPROX.


=========
Fused_lasso_b.h
=========

This class inherits the class ALM_Katyusha. It solves fused lasso problem by the method IPALM_Katyusha.


=========
Fused_lasso_c.h
=========

This class inherits the class PDCD. It solves fused lasso problem by the method SMART_CD.



=========
Fused_lasso_d.h
=========

This class inherits the class ALM_APG. It solves fused lasso problem by the method ASGARD_DL.





=========
L_1_Lasso_a.h
=========

This class inherits the class ALM_APPROX. It solves Least absolute deviation problem by the method IPALM_APPROX.



=========
L_1_Lasso_c.h
=========

This class inherits the class PDCD. It solves Least absolute deviation problem by the method SMART_CD.



=========
L_1_Lasso_d.h
=========

This class inherits the class ALM_APG. It solves Least absolute deviation problem by the method ASGARD_DL.




=========
SMSVM_a.h
=========

This class inherits the class ALM_APPROX. It solves soft margin support vector machine problem by the method IPALM_APPROX.


=========
SMSVM_b.h
=========

This class inherits the class ALM_Katyusha. It solves soft margin support vector machine problem by the method IPALM_Katyusha.


=========
SMSVM_c.h
=========

This class inherits the class PDCD_SMSVM. It solves soft margin support vector machine problem by the method SMART_CD.



=========
SMSVM_d.h
=========

This class inherits the class ALM_APG. It solves soft margin support vector machine problem by the method ASGARD_DL.




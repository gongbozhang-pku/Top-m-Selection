# Source Codes of Top-m Selection Problem

This repository contains the source codes used in the following paper:

* Asymptotically Optimal Sampling Policy for Selecting Top-m Alternatives, Gongbo Zhang, Yijie Peng<sup>*</sup>, Jianghua Zhang, Enlu Zhou (2022) Submitted to *INFORMS Journal on Computing*.

### Disclaimer

LICENSE: Redistribution and use in source and binary forms, with or without modification, are permitted, under the terms of the [BSD license](./BSD_License.txt).

WARNING: These codes are written only for the purpose of demonstration and verification. While the correctness has been carefully checked, the quality such as standardability, clarity, generality, and efficiency has not been well considered.

## 1 Introduction

### 1.1 The `/Exp1-4` folder contains the MATLAB implementations of the numerical experiments in §5.1 of Zhang et al. (2022)

* The function `/Exp1-4/AOAPm.m` is the proposed policy of Zhang et al. (2022);
* The function `/Exp1-4/EAm.m` is a compared policy called equal allocation (EA);
* The function `/Exp1-4/OCBAmsa.m` is a modified compared policy from Chen et al. (2008), where a better separating parameter that improves performance of the original policy is sequentially implemented;
* The function `/Exp1-4/OCBAmjia.m` is a compared policy from Zhang et al. (2012,2015);
* The function `/Exp1-4/OCBAss.m` is a compared policy from Gao and Chen (2015);
* The function `/Exp1-4/OCBASS.m` is a compared policy from Gao and Chen (2016).

### 1.2 The `/Exp5` folder contains the MATLAB implementations of the numerical experiments in §5.2 of Zhang et al. (2022)

* The function `/Exp5/AOAPm.m` is the proposed policy of Zhang et al. (2022);
* The function `/Exp5/EAm.m` is a compared policy called equal allocation (EA);
* The function `/Exp5/OCBAmsa.m` is a modified compared policy from Chen et al. (2008), where a better separating parameter that improves performance of the original policy is sequentially implemented;
* The function `/Exp5/OCBAmjia.m` is a compared policy from Zhang et al. (2012,2015);
* The function `/Exp5/OCBAss.m` is a compared policy from Gao and Chen (2015);
* The function `/Exp5/OCBASS.m` is a compared policy from Gao and Chen (2016);
* The function `/Exp5/evacuation.m` is a simulator of the evacuation network;
* The function `/Exp5/Untitled.m` estimates the mean performance of each evacuation plan.

### 1.3 The `/App` folder contains the MATLAB implementations of the numerical experiments in the online Appendix of Zhang et al. (2022)

#### 1.3.1 The `/App/OCBA` folder contains codes for numerical experiments in §A.1 of the online Appendix
   
* The function `/App/OCBA/OCBAm.m` is a compared policy from Chen et al. (2008);

* The function `/App/OCBA/OCBAms.m` is a modified compared policy from Chen et al. (2008), where the original policy is implemeted in a sequential manner.

#### 1.3.2 The `/App/Sampling ratio` folder contains codes for numerical experiments in §A.3 of the online Appendix

* The function `/App/Sampling ratio/AOAPm.m` is the proposed policy of Zhang et al. (2022), where the output is the sampling ratio for each simulation budget.

#### 1.3.3 The `/App/trueparameter` folder contains codes for numerical experiments in §A.4.1 of the online Appendix
  
* The function `/App/trueparameter/OCBAss.m` is a modified compared policy from Gao and Chen (2015), where true parameters are used for implementation;

* The function `/App/trueparameter/OCBASSSS.m` is a modified compared policy from Gao and Chen (2016), where true parameters are used for implementation.

#### 1.3.4 The `/App/Non-normal` folder contains codes for numerical experiments in §A.4.2 of the online Appendix

* The function `/App/Non-normal/Bernoulli/AOAPm.m` is the proposed policy of Zhang et al. (2022), which is implemented under Bernoulli sampling distributions;

* The function `/App/Non-normal/Bernoulli/EAm.m` is a compared policy called equal allocation (EA), which is implemented under Bernoulli sampling distributions;
   
* The function `/App/Non-normal/Bernoulli/OCBAmsa.m` is a modified compared policy from Chen et al. (2008), which is implemented under Bernoulli sampling distributions, and a better separating parameter that improves performance of the original policy is sequentially implemented;

* The function `/App/Non-normal/Bernoulli/OCBAmjia.m` is a compared policy from Zhang et al. (2012,2015), which is implemented under Bernoulli sampling distributions;

* The function `/App/Non-normal/Bernoulli/OCBAss.m` is a compared policy from Gao and Chen (2015), which is implemented under Bernoulli sampling distributions;

* The function `/App/Non-normal/Bernoulli/OCBASS.m` is a compared policy from Gao and Chen (2016), which is implemented under Bernoulli sampling distributions.

* The function `/App/Non-normal/Exponential/AOAPm.m` is the proposed policy of Zhang et al. (2022), which is implemented under Exponential sampling distributions;

* The function `/App/Non-normal/Exponential/EAm.m` is a compared policy called equal allocation (EA), which is implemented under Exponential sampling distributions;
   
* The function `/App/Non-normal/Exponential/OCBAmsa.m` is a modified compared policy from Chen et al. (2008), which is implemented under Exponential sampling distributions, and a better separating parameter that improves performance of the original policy is sequentially implemented;

* The function `/App/Non-normal/Exponential/OCBAmjia.m` is a compared policy from Zhang et al. (2012,2015), which is implemented under Exponential sampling distributions;
   
* The function `/App/Non-normal/Exponential/OCBAss.m` is a compared policy from Gao and Chen (2015), which is implemented under Exponential sampling distributions;

* The function `/App/Non-normal/Exponential/OCBASS.m` is a compared policy from Gao and Chen (2016), which is implemented under Exponential sampling distributions.

* The function `/App/Non-normal/Queueing/queueing.m` is a simulator of the two-node tandem queueing system;

* The function `/App/Non-normal/Queueing/truevalue.m` estimates the mean performance of each worker allocation plan.

#### 1.3.5 The `/App/Top-m arms identification` folder contains codes for numerical experiments in §A.4.3 of the online Appendix

* The function `/App/Top-m arms identification/Normal/pSAR.m` is the proposed policy of Bubeck et al. (2013), which is implemented under Normal sampling distributions;

* The function `/App/Top-m arms identification/Normal/pSR.m` is the proposed policy of Audibert et al. (2010) and is then modified by Bubeck et al. (2013), which is implemented under Normal sampling distributions;

* The function `/App/Top-m arms identification/Normal/pGapE.m` is the proposed policy of Gabillon et al. (2011) and is then modified by Bubeck et al. (2013), which is implemented under Normal sampling distributions;

* The function `/App/Top-m arms identification/Bernoulli/pSAR.m` is the proposed policy of Bubeck et al. (2013), which is implemented under Bernoulli sampling distributions;

* The function `/App/Top-m arms identification/Normal/pSR.m` is the proposed policy of Audibert et al. (2010) and is then modified by Bubeck et al. (2013), which is implemented under Bernoulli sampling distributions;

* The function `/App/Top-m arms identification/Normal/pGapE.m` is the proposed policy of Gabillon et al. (2011) and is then modified by Bubeck et al. (2013), which is implemented under Bernoulli sampling distributions.

## 2 Installation

* The codes were written and run in MATLAB R2021a, on Windows 10 Home 64-bit OS, with Intel i5-11300H CPU @ 3.10 GHz, 16 GB RAM.

* To install the MATLAB codes, just copy the entire folder `/Exp1-4`, `/Exp5` and `/App`, respectively, into your MATLAB directory, or change the path of MATLAB to the folder `/Exp1-4`, `/Exp5` and `/App`, respectively.

## 3 Details on Numerical Experiments

### 3.1 Numerical Experiments in §5.1 of Zhang et al. (2022)

Get into folder `/Exp1-4`. Run `/Exp1-4/EAm.m` for *EA*, `/Exp1-4/OCBAmsa.m` for *OCBAm*, `/Exp1-4/OCBAmjia.m` for *OCBAm+*, `/Exp1-4/OCBAss.m` for *OCBAss*, `/Exp1-4/OCBASS.m` for *OCBASS*, and `/Exp1-4/AOAPm.m` for *AOAm*.

* Set corresponding input parameters for all policies in folder `/Exp1-4`;
* The policies can adjust ascending or descending;
* The policies can incoporate prior information or not;
* Update parameters in the same way when comparing;
* The independent macro experiments can parallelly run.
  
### 3.2 Numerical Experiments in §5.2 of Zhang et al. (2022)

Get into folder `/Exp5`. Run `/Exp5/EAm.m` for *EA*, `/Exp5/OCBAmsa.m` for *OCBAm*, `/Exp5/OCBAmjia.m` for *OCBAm+*, `/Exp5/OCBAss.m` for *OCBAss*, `/Exp5/OCBASS.m` for *OCBASS*, and `/Exp5/AOAPm.m` for *AOAm*.

* Set corresponding input parameters for all policies in folder `/Exp5`;
* The input parameter *truemu* of each policy is calculated by running `/Exp5/Untitled1.m`, which will call `/Exp5/evacuation.m` during the running;
* The independent macro experiments can parallelly run.

### 3.3 Numerical Experiments in §A.1 of the online Appendix

Get into folder `/App/OCBA`. Run `/Exp1-4/EAm.m` for *EA*, `/App/OCBA/OCBAm.m` for *OCBAm(two-stage)*, `/App/OCBA/OCBAms.m` for *OCBAm(sequential)*, `/Exp1-4/OCBAm.m` for *OCBAm*, `/Exp1-4/AOAPm.m` for *OCBASS*, and `/Exp1-4/AOAPm.m` for *AOAm*.

* Set corresponding input parameters for all policies in folder `/App/OCBA` and `/Exp1-4`.

### 3.4 Numerical Experiments in §A.3 of the online Appendix

Get into folder `/App/Sampling ratio`. Run `/AOAPm.m` for its sampling ratios.

* Set corresponding input parameters for the policy.

### 3.5 Numerical Experiments in §A.4.1 of the online Appendix

Get into folder `/App/trueparameter`. Run `/Exp1-4/EAm.m` for *EA*, `/Exp1-4/OCBAm.m` for *OCBAm*, `/Exp1-4/OCBAmjia.m` for *OCBAm+*, `/App/trueparameter/OCBAss.m` for *OCBAsst*, `/App/trueparameter/OCBASSS.m` for *OCBASSt*, and `/Exp1-4/AOAPm.m` for *AOAm*.

* Set corresponding input parameters for all policies in folder `/App/trueparameter` and `/Exp1-4`.

### 3.6 Numerical Experiments in §A.4.2 of the online Appendix

#### 3.6.1 Numerical Experiments in §A.4.2.1 of the online Appendix

Get into folder `/App/Non-normal/Bernoulli`. Run `/App/Non-normal/Bernoulli/EAm.m` for *EA*, `/App/Non-normal/Bernoulli/OCBAm.m` for *OCBAm*, `/App/Non-normal/Bernoulli/OCBAmjia.m` for *OCBAm+*, `/App/Non-normal/Bernoulli/OCBAss.m` for *OCBAss*, `/App/Non-normal/Bernoulli/OCBASSS.m` for *OCBASS*, and `/App/Non-normal/Bernoulli/AOAPm.m` for *AOAm*.

* Set corresponding input parameters for all policies in folder `/App/Non-normal/Bernoulli`;
* The independent macro experiments can parallelly run.

#### 3.6.2 Numerical Experiments in §A.4.2.2 of the online Appendix

Get into folder `/App/Non-normal/Exponential`. Run `/App/Non-normal/Exponential/EAm.m` for *EA*, `/App/Non-normal/Exponential/OCBAm.m` for *OCBAm*, `/App/Non-normal/Exponential/OCBAmjia.m` for *OCBAm+*, `/App/Non-normal/Exponential/OCBAss.m` for *OCBAss*, `/App/Non-normal/Exponential/OCBASSS.m` for *OCBASS*, and `/App/Non-normal/Exponential/AOAPm.m` for *AOAm*.

* Set corresponding input parameters for all policies in folder `/App/Non-normal/Exponential`;
* The independent macro experiments can parallelly run.

#### 3.6.3 Numerical Experiments in §A.4.2.3 of the online Appendix

Get into folder `/Exp5`. Run `/Exp5/EAm.m` for *EA*, `/Exp5/OCBAmsa.m` for *OCBAm*, `/Exp5/OCBAmjia.m` for *OCBAm+*, `/Exp5/OCBAss.m` for *OCBAss*, `/Exp5/OCBASS.m` for *OCBASS*, and `/Exp5/AOAPm.m` for *AOAm*.

* Set corresponding input parameters for all policies in folder `/Exp5`;
* The input parameter *truemu* of each policy is calculated by running `/App/Non-normal/Exponential/truevalue.m`, which will call `/App/Non-normal/Exponential/queueing.m` during the running;
* The independent macro experiments can parallelly run.

### 3.7 Numerical Experiments in §A.4.3 of the online Appendix

#### 3.7.1 Numerical Experiments in §A.4.3.1 and §A.4.3.2 of the online Appendix

Get into folder `/App/Top-m arms identification/Normal`. Run `/Exp1-4/EAm.m` for *EA*, `/Exp1-4/AOAPm.m` for *AOAm*, `/App/Top-m arms identification/Normal/pGapE.m` for *Gap-E*, `/App/Top-m arms identification/Normal/pSAR.m` for *SAR*, `/App/Top-m arms identification/Normal/pSR.m` for *SR*.

* Set corresponding input parameters for all policies in folder `/App/Top-m arms identification/Normal`;
* The independent macro experiments can parallelly run.

#### 3.7.2 Numerical Experiments in §A.4.3.3 of the online Appendix

Get into folder `/App/Top-m arms identification/Bernoulli`. Run `/App/Non-normal/Bernoulli/EAm.m` for *EA*, `/App/Non-normal/Bernoulli/AOAPm.m` for *AOAm*, `/App/Top-m arms identification/Bernoulli/pGapE.m` for *Gap-E*, `/App/Top-m arms identification/Bernoulli/pSAR.m` for *SAR*, `/App/Top-m arms identification/Bernoulli/pSR.m` for *SR*.

* Set corresponding input parameters for all policies in folder `/App/Top-m arms identification/Normal`;
* The independent macro experiments can parallelly run.

## References

* Chen CH, He D, Fu M, Lee LH (2008) Efficient Simulation Budget al.location for Selecting an Optimal Subset. INFORMS Journal on Computing 20(4): 579–595.
* Zhang S, Lee LH, Chew EP, Chen CH, Jen HY (2012) An Improved Simulation Budget al.location Procedure to Efficiently Select the Optimal Subset of Many Alternatives. 2012 IEEE International Conference on Automation Science and Engineering (CASE), 230–236 (IEEE).
* Zhang S, Lee LH, Chew EP, Xu J, Chen CH (2015) A Simulation Budget al.location Procedure for Enhancing the Efficiency of Optimal Subset Selection. IEEE Transactions on Automatic Control 61(1): 62–75.
* Gao S, Chen W (2015) A Note on the Subset Selection for Simulation Optimization. Proceedings of the 2015 Winter Simulation Conference, 3768–3776 (IEEE).
* Gao S, Chen W (2016) A New Budget al.location Framework for Selecting Top Simulated Designs. IIE Transactions 48(9): 855–863.
* Bubeck S, Wang T, Viswanathan N. Multiple identifications in multi-armed bandits[C] // In International Conference on Machine Learning (ICML). PMLR, 2013: 258-265.
* Gabillon V, Ghavamzadeh M, Lazaric A, et al. Multi-bandit best arm identification[J]. In Advances in Neural Information Processing Systems (NIPS), 2011, 24.
* Audibert J Y, Bubeck S, Munos R. Best arm identification in multi-armed bandits[C]// In Proceedings of the 23rd Annual Conference on Learning Theory (COLT). Citeseer, 2010: 41-53.
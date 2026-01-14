This repository contains the C++ source code and the speckle images (whose underlying displacement field is double-step) utlized in our paper "A criterion for assessing spatial resolution in digital image correlation: applications to conventional, Gaussian-weighted, and deep learning-based methods", published in Optics and Lasers in Engineering. 

Digital image correlation (DIC) is a full-field optical metrology widely used for measuring shape, motion, and deformation. Accurate characterization of spatial resolution is essential for assessing the performance of DIC methods. 

Inspired by the classical Rayleigh criterion, we considers a double-step displacement field and defines the spatial resolution as the separation at which two adjacent displacement discontinuities can just be resolved in the measured strain field. Based on this definition, theoretical models are developed to quantify the spatial resolution. 

(1)	For conventional DIC with a subset size of $2M+1$, the spatial resolution for the first-order shape function is approximately $2M+2$, while that for the second-order shape function is about $0.66(2M+1)$. 

(2)	For Gaussian-weighted DIC with a weighting radius of $R$, the spatial resolution is approximately $2R$ for the first-order shape function and $1.53 R$ for the second-order shape function. 

(3)	The proposed criterion is further applied to deep learning-based DIC methods (e.g., U-DICNet, DICTr, and USDICNet) to demonstrate its universality. 

In summary, we establish a physically grounded and quantitatively reliable framework for evaluating spatial resolution in conventional, weighted, and  learning-based DIC methods, with potential applications in performance assessment, algorithm optimization, and parameter selection.

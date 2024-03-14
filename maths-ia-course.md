---
layout: page
title: "The Mathematics of IA"
description: "10 lectures"
header-img: "img/boat-waves.jpg"
---

This series consists of 10 sessions, each lasting 2 hours, focused on the mathematics of machine learning. It outlines the primary concepts without delving into the intricacies of the proofs. Clicking on the title of each session provides access to the transcript. Additionally, there are basic notes available to guide you through the structure and progression of the content.

Here's your corrected and improved markdown:

### Course #1 - [Smooth Optimization](maths-ia-course-data/transcripts/1-smooth-optim.pdf)

**Content:**
- Introduction and motivation
- Gradients, Jacobians, Hessians
- Gradient descent and acceleration
- Stochastic Gradient Descent (SGD)

**Materials:**
- [Notebook on Regression](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_2_regression.ipynb)
- [Notebook on Classification](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_3_classification.ipynb)
- My course notes: [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Exercises Sheet](maths-ia-course-data/material/td-derivatives.pdf)

**Bibliography:**
- [Convex Optimization](http://stanford.edu/~boyd/cvxbook/), by Boyd and Vandenberghe
- [Introduction to Nonlinear Optimization: Theory, Algorithms, and Applications](https://sites.google.com/site/amirbeck314/books), by Amir Beck

### Course #2 - [From Smooth to Non-Smooth Optimization](maths-ia-course-data/transcripts/2-non-smooth-optim.pdf)

**Content:**
- Proofs of gradient descent and acceleration
- Linear models and regularization
- Ridge versus Lasso
- ISTA Algorithm

**Materials:**
- [Notebook on Linear Regression](https://nbviewer.jupyter.org/github/gpeyre/numerical-tours/blob/master/python/ml_2_regression.ipynb) (specifically the Lasso part)
- [Notebook on Interior Point Methods](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/optim_6_interior_points.ipynb)
- My course notes: [The Mathematical Tours of Signal Processing](https://mathematical-tours.github.io/book/)

**Bibliography:**
- [Course Notes on Convexity](https://who.rocq.inria.fr/Vincent.Duval/files/convexity.pdf) by Vincent Duval
- [Introduction to Nonlinear Optimization: Theory, Algorithms, and Applications](https://sites.google.com/site/amirbeck314/books), by Amir Beck

### Course #3 - [Lasso and Compressed Sensing](maths-ia-course-data/transcripts/3-compressed-sensing.pdf)

**Content:**
- Examples of non-smooth functionals (Lasso, TV regularization, constraints)
- Subgradient and proximal operators
- Forward-backward splitting, connection with FISTA
- ADMM, Douglas-Rachford (DR), Primal-Dual
- Compressive sensing theory

**Materials:**
- My course notes: [The Mathematical Tours of Signal Processing](https://mathematical-tours.github.io/book/),
- [Notebook on Douglas-Rachford Proximal Method](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/optim_4_condat_dr.ipynb)
- [Proximal Operators Repository](http://proximity-operator.net/) (including Python code)
- [Non-Smooth Optimization Slides](maths-ia-course-data/material/nonsmooth-optimization-slides.pdf)
- [Compressed Sensing Slides](maths-ia-course-data/material/compressed-sensing-slides.pdf)

**Bibliography:**
- *A Mathematical Introduction to Compressive Sensing* by Foucart, Simon and Rauhut, Holger (advanced)
- [Convex Optimization](http://stanford.edu/~boyd/cvxbook/), by Boyd and Vandenberghe
- [Proximal Algorithms](http://stanford.edu/~boyd/papers/prox_algs.html), by N. Parikh and S. Boyd


### Course #4 - [Kernel, Perceptron, CNN, and Transformers](maths-ia-course-data/transcripts/4-neural-networks.pdf)

**Content:**
- Transition from ridge regression to kernels
- Multilayer Perceptron (MLP)
- Convolutional Neural Networks (CNN)
- ResNet architecture
- Transformer models

**Materials:**
- [Slides on deep learning](maths-ia-course-data/material/course-ml-deep-learning.pdf)
- My course notes on [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Notebook on Multilayer Perceptron and Autograd](http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/python/ml_5_mlp.ipynb)

**Bibliography:**
- [The Elements of Statistical Learning](https://statweb.stanford.edu/~tibs/ElemStatLearn/), by Jerome H. Friedman, Robert Tibshirani, and Trevor Hastie
- [Machine Learning: A Probabilistic Perspective](http://www.cs.ubc.ca/~murphyk/MLbook/index.html), by Kevin Patrick Murphy (covers the theory of ML)


### Course #5 - [Deep Learning: Theory and Numerics](maths-ia-course-data/transcripts/5-neural-networks-theory.pdf)

**Content:**
- Review of MLP and its variants (CNN, ResNet)
- Theoretical framework of two-layer MLPs
- Gradient and Jacobians in neural networks
- Introduction to backpropagation

**Materials:**
- [Slides on deep learning](maths-ia-course-data/material/course-ml-deep-learning.pdf)
- [Slides on automatic differentiation](https://www.dropbox.com/s/4jr93tzz2uxgcwv/course-ml-autodiff.pdf?dl=0)
- My course notes on [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Notebook on deep learning](http://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_6_nn.ipynb)
- [Notebook on texture synthesis with deep networks](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_8_deep_texture_synthesis.ipynb)

**Bibliography:**
- [Christopher Olah's Blog](http://colah.github.io/)


### Course #6 - [Differential Programming](maths-ia-course-data/transcripts/6-autodiff.pdf)

**Content:**
- Recap on Gradient and Jacobian
- Forward and reverse mode automatic differentiation
- Introduction to PyTorch
- The adjoint method in computational mathematics

**Materials:**
- [Slides on automatic differentiation](maths-ia-course-data/material/course-ml-autodiff.pdf)
- My course notes on [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Code example: Multilayer perceptron and autograd](http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/python/ml_5_mlp.ipynb)

**Bibliography:**
- [Reverse-mode Automatic Differentiation: A Tutorial](https://rufflewind.com/2016-12-30/reverse-mode-automatic-differentiation)

### Course #7 - [Sampling and Diffusion Models](maths-ia-course-data/transcripts/7-diffusion-models.pdf)

**Content:**
- Refresher on Stochastic Gradient Descent (SGD)
- Introduction to Langevin dynamics
- Overview of diffusion models

**Materials:**
- [Numerical tour on diffusion models](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_12_diffusion_models.ipynb)
- Course notes on [Diffusion Models](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML-DiffusionModels.pdf)


### Course #8 - [LLM and Generative AI](maths-ia-course-data/transcripts/8-llm.pdf)

**Content:**
- Overview of different generative model concepts
- Introduction to generative models (VAE, GANs, U-Net, diffusion)
- Semi-supervised learning and next token prediction
- Tokenizers
- Transformer architectures, Flash attention
- State space models

**Materials:**
- [Introductory slides on Generative AI](https://speakerdeck.com/gpeyre/generative-ai-practice-and-theory)

**Bibliography:**
- [Andrej Karpathy's video on tokenization](https://www.youtube.com/watch?v=zduSFxRajkE)
- [Byte Pair Encoding](https://en.wikipedia.org/wiki/Byte_pair_encoding) (Wikipedia)
- Online [tokenizer demo by OpenAI](https://platform.openai.com/token

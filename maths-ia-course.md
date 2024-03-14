---
layout: page
title: "The Mathematics of IA"
description: "10 lectures"
header-img: "img/boat-waves.jpg"
---

This series consists of 10 sessions, each lasting 2 hours, focused on the mathematics of machine learning. It outlines the primary concepts without delving into the intricacies of the proofs. Clicking on the title of each session provides access to the transcript. Additionally, there are basic notes available to guide you through the structure and progression of the content.

### Course #1 - [Smooth Optimization](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/1-smooth-optim.pdf)

**Content:**
- Introduction and motivation
- Gradients, Jacobians, Hessians
- Gradient descent and acceleration
- Stochastic Gradient Descent (SGD)

**Materials:**
- [Notebook on Regression](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_2_regression.ipynb)
- [Notebook on Classification](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_3_classification.ipynb)
- My course notes: [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Exercises Sheet](https://mathematical-tours.github.io/maths-ia-course-data/material/td-derivatives.pdf)

**Bibliography:**
- [Convex Optimization](http://stanford.edu/~boyd/cvxbook/), by Boyd and Vandenberghe
- [Introduction to Nonlinear Optimization: Theory, Algorithms, and Applications](https://sites.google.com/site/amirbeck314/books), by Amir Beck

### Course #2 - [From Smooth to Non-Smooth Optimization](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/2-non-smooth-optim.pdf)

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

### Course #3 - [Lasso and Compressed Sensing](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/3-compressed-sensing.pdf)

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
- [Non-Smooth Optimization Slides](https://mathematical-tours.github.io/maths-ia-course-data/material/nonsmooth-optimization-slides.pdf)
- [Compressed Sensing Slides](https://mathematical-tours.github.io/maths-ia-course-data/material/compressed-sensing-slides.pdf)

**Bibliography:**
- *A Mathematical Introduction to Compressive Sensing* by Foucart, Simon and Rauhut, Holger (advanced)
- [Convex Optimization](http://stanford.edu/~boyd/cvxbook/), by Boyd and Vandenberghe
- [Proximal Algorithms](http://stanford.edu/~boyd/papers/prox_algs.html), by N. Parikh and S. Boyd


### Course #4 - [Kernel, Perceptron, CNN, and Transformers](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/4-neural-networks.pdf)

**Content:**
- Transition from ridge regression to kernels
- Multilayer Perceptron (MLP)
- Convolutional Neural Networks (CNN)
- ResNet architecture
- Transformer models

**Materials:**
- [Slides on deep learning](https://mathematical-tours.github.io/maths-ia-course-data/material/course-ml-deep-learning.pdf)
- My course notes on [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Notebook on Multilayer Perceptron and Autograd](http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/python/ml_5_mlp.ipynb)

**Bibliography:**
- [The Elements of Statistical Learning](https://statweb.stanford.edu/~tibs/ElemStatLearn/), by Jerome H. Friedman, Robert Tibshirani, and Trevor Hastie
- [Machine Learning: A Probabilistic Perspective](http://www.cs.ubc.ca/~murphyk/MLbook/index.html), by Kevin Patrick Murphy (covers the theory of ML)


### Course #5 - [Deep Learning: Theory and Numerics](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/5-neural-networks-theory.pdf)

**Content:**
- Review of MLP and its variants (CNN, ResNet)
- Theoretical framework of two-layer MLPs
- Gradient and Jacobians in neural networks
- Introduction to backpropagation

**Materials:**
- [Slides on deep learning](https://mathematical-tours.github.io/maths-ia-course-data/material/course-ml-deep-learning.pdf)
- [Slides on automatic differentiation](https://www.dropbox.com/s/4jr93tzz2uxgcwv/course-ml-autodiff.pdf?dl=0)
- My course notes on [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Notebook on deep learning](http://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_6_nn.ipynb)
- [Notebook on texture synthesis with deep networks](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_8_deep_texture_synthesis.ipynb)

**Bibliography:**
- [Christopher Olah's Blog](http://colah.github.io/)


### Course #6 - [Differential Programming](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/6-autodiff.pdf)

**Content:**
- Recap on Gradient and Jacobian
- Forward and reverse mode automatic differentiation
- Introduction to PyTorch
- The adjoint method in computational mathematics

**Materials:**
- [Slides on automatic differentiation](https://mathematical-tours.github.io/maths-ia-course-data/material/course-ml-autodiff.pdf)
- My course notes on [Optimization for Machine Learning](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML.pdf)
- [Code example: Multilayer perceptron and autograd](http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/python/ml_5_mlp.ipynb)

**Bibliography:**
- [Reverse-mode Automatic Differentiation: A Tutorial](https://rufflewind.com/2016-12-30/reverse-mode-automatic-differentiation)

### Course #7 - [Sampling and Diffusion Models](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/7-diffusion-models.pdf)

**Content:**
- Refresher on Stochastic Gradient Descent (SGD)
- Introduction to Langevin dynamics
- Overview of diffusion models

**Materials:**
- [Numerical tour on diffusion models](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/ml_12_diffusion_models.ipynb)
- Course notes on [Diffusion Models](https://mathematical-tours.github.io/book-sources/optim-ml/OptimML-DiffusionModels.pdf)

### Course #8 - [LLM and Generative AI](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/8-llm.pdf)

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
- Online [tokenizer demo by OpenAI](https://platform.openai.com/tokenizer)
- [Rotary Position Embedding paper](https://arxiv.org/abs/2104.09864)
- Codes: [Flash attention](https://github.com/Dao-AILab/flash-attention), [xFormers](https://github.com/facebookresearch/xformers), [Triton](https://openai.com/research/triton)
- [Theory paper on Flash Attention](https://arxiv.org/abs/2402.07443)
- [Mamba paper](https://arxiv.org/abs/2312.00752), Blog [on Mamba](https://jackcook.com/2024/02/23/mamba.html) (SSM), [Parallel Prefix Sum algorithm](https://en.wikipedia.org/wiki/Prefix_sum)

### Course #9 - [Generative Models](https://mathematical-tours.github.io/maths-ia-course-data/transcripts/9-generative.pdf)

**Content:**
- Understanding generative models as density fitting techniques.
- Basics of Maximum Likelihood Estimation and f-divergences.
- Gaussian mixtures and the Expectation-Maximization algorithm.
- Variational Autoencoders (VAE).
- Introduction to Normalizing Flows.
- Generative Adversarial Networks (GANs), Wasserstein GANs (WGANs).
- Diffusion Models.

**Materials:**
- [Course notes on the EM algorithm for mixtures](https://mathematical-tours.github.io/maths-ia-course-data/material/em-algo.pdf)
- [Tutorial code on generative models (EM/VAE/NF/GANs)](https://mathematical-tours.github.io/maths-ia-course-data/material/generative-models.ipynb)

**Bibliography:**
- [Tutorial on Variational Autoencoders](https://arxiv.org/abs/1606.05908)
- [Tutorial on Normalizing Flows](https://pyro.ai/examples/normalizing_flows_i.html)
- [Slides on Normalizing Flows](https://www.shakirm.com/slides/DeepGenModelsTutorial.pdf)

### Course #10 - Optimal Transport

**Content:**
- Introduction to Monge and Kantorovich formulations.
- The Sinkhorn algorithm.
- Training of generative models.
- Duality and Wasserstein GANs.

**Materials:**
- [Slides on Optimal Transport](https://speakerdeck.com/gpeyre/numerical-optimal-transport-1?slide=120)
- [Notebook on Linear Programming for Optimal Transport](http://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/optimaltransp_1_linprog.ipynb)
- [Notebook on the Sinkhorn algorithm](http://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/optimaltransp_5_entropic.ipynb)

**Bibliography:**
- [Computational Optimal Transport](https://arxiv.org/abs/1803.00567), by Gabriel Peyré and Marco Cuturi
- [Optimal Transport for Applied Mathematicians](https://www.math.u-psud.fr/~filippo/OTAM-cvgmt.pdf), by Filippo Santambrogio (advanced)
- [Python POT (Python Optimal Transport) toolbox](https://pythonot.github.io/)
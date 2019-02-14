% Fast and Accurate Bilateral Filtering using 
% Gaussian-Polynomial Approximation
%%%%%%%%%%%%%%%%%%%%%%%%
%  f           : input 8-bit grayscale image
%  g           : output of Gaussian bilateral filter
%  b           : output of box bilateral filter
%  sigmas   : width of spatial Gaussian
%  B           : width of spatial box
%  sigmar   : width of range Gaussian
%  eps        : kernel accuracy

%  Authors:    Kunal Chaudhury, Swapnil Dabhade.
%  Date:       March 25, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all force;
f  =  double( imread('./images/cameraman.tif') );
% filter parameters
sigmar = 40;
eps    = 1e-3;
% compute Gaussian bilateral filter
sigmas = 3;
[g,Ng] = GPA(f, sigmar, sigmas, eps, 'Gauss');     
% compute box bilateral filter
B = 9;
[b,Nb] = GPA(f, sigmar, B, eps, 'box');    
% results
figure('Units','normalized','Position',[0 0.5 1 0.5]);
colormap gray,
subplot(1,3,1), imshow(uint8(f)),
title('Input', 'FontSize', 10), axis('image', 'off');
subplot(1,3,2), imshow(uint8(g)),
title('Gaussian Bilateral Filter','FontSize', 10),
axis('image', 'off');
subplot(1,3,3), imshow(uint8(b)),
title('Box Bilateral Filter','FontSize', 10),
axis('image', 'off');

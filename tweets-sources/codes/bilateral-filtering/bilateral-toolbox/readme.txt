This is the MATLAB implementation of the fast and accurate 
approximation of the Bilateral filtering described in the following paper:

[1] K. N. Chaudhury and S. Dabhade, "Fast and provably accurate 
bilateral filtering", IEEE Transactions on Image Processing, vol. 25, 
no. 6, pp. 2519-2528, 2016.

Authors  : Kunal Chaudhury and Swapnil Dabhade.
GUI       :    Anmol Popli.
Date       : March 2016

To run the software:
============

(1) Gaussian Bilateral filter:
    [g,Ng] = GPA(f, sigmar, sigmas, eps,'Gauss')
%% Input:
    % f              : input image 
    % sigmar     : width of range Gaussian
    % sigmas     : width of spatial Gaussian
    % eps          : kernel approximation accuracy
%% Output:
    % g             : output image
    % Ng          : approximation order

(2) Box bilateral filter:
    [b,Nb] = GPA(f, sigmar, B, eps,'box')
%% Input:
    % f              : input image 
    % sigmar     : width of range Gaussian
    % B             : width of box kernel
    % eps           : kernel approximation accuracy
%% Output:
    % b               : output image
    % Nb            : approximation order

(In either case, the range kernel is a Gaussian with width sigmar).

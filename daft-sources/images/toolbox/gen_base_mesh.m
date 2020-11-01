function [vertex,face] = gen_base_mesh(type);

% gen_base_mesh - generate a simple triangulation.
%
%   [vertex,face] = gen_base_mesh(type);
%
%   'type' can be one 'triangle', 'square', 'square1', 'L', 'L1', 'tetra',
%   'oct', or 'ico'.
%
%   Copyright (c) 2004 Gabriel Peyré

switch  lower(type)
    case 'triangle'
        vertex = [0,0,0; 1,0,0; 0.5,1/sqrt(2),0];
        face = [1,2,3];
    case 'square'
        vertex = [0,0,0; 1,0,0; 1,1,0; 0,1,0; 0.5,0.5,0];
        face = [1,2,5; 2,3,5; 3,4,5; 4,1,5];
    case 'square1'
        vertex = [0,0,0; 1,0,0; 1,1,0; 0,1,0];
        face = [1,2,3; 3,4,1];
    case 'l'
        vertex = [  0,0,0; 1,0,0; 2,0,0;
                    0,1,0; 1,1,0; 2,1,0;
                    0,2,0; 1,2,0;
                    0,3,0; 1,3,0];
        face = [1,2,5; 1,5,4;
                2,3,6; 2,6,5;
                4,5,8; 4,8,7;
                7,8,10; 7,10,9];
    case 'l1'
        vertex = [  0,0,0; 1,0,0; 2,0,0;
                    0,1,0; 1,1,0; 2,1,0;
                    0,2,0; 1,2,0;
                    0.5,0.5,0; 1.5,0.5,0; 0.5,1.5,0];
        face = [1,2,9; 2,5,9; 5,4,9; 4,1,9; 
                2,3,10; 3,6,10; 6,5,10; 5,2,10;
                4,5,11; 5,8,11; 8,7,11; 7,4,11];  
    case 'tetra',
        sqrt_3 = 0.5773502692;
        vertex = [  sqrt_3,  sqrt_3,  sqrt_3 ;
                -sqrt_3, -sqrt_3,  sqrt_3 ;
                -sqrt_3,  sqrt_3, -sqrt_3 ;
                 sqrt_3, -sqrt_3, -sqrt_3 ]; 
        face = [ 1, 2, 3;
                1, 4, 2;
                3, 2, 4;
                4, 1, 3 ]; 
    case 'oct',
        vertex = [  1,  0,  0 ;
              -1,  0,  0 ;
               0,  1,  0 ;
               0, -1,  0 ;
               0,  0,  1 ;
               0,  0, -1 ];
        face = [ 1 5 3 ;
              3 5 2 ;
              2 5 4 ;
              4 5 1 ;
              1 3 6 ;
              3 2 6 ;
              2 4 6 ;
              4 1 6 ];    
    case 'ico',
        tau = 0.8506508084;
        one = 0.5257311121;
        vertex = [  tau,  one,    0;
                -tau,  one,    0
                -tau, -one,    0;
                tau, -one,    0;
                one,   0 ,  tau;
                one,   0 , -tau;
                -one,   0 , -tau;
                -one,   0 ,  tau;
                0 ,  tau,  one;
                0 , -tau,  one;
                0 , -tau, -one;
                0 ,  tau, -one ];
        face = [  5,  9,  8 ;
               5,  8, 10 ;
               6,  7, 12 ;
               6, 11,  7 ;
               1,  5,  4 ;
               1,  4,  6 ;
               3,  8,  2 ;
               3,  2,  7 ;
               9,  1, 12 ;
               9, 12,  2 ;
              10, 11,  4 ;
              10,  3, 11 ;
               9,  5,  1 ;
              12,  1,  6 ;
               5, 10,  4 ;
               6,  4, 11 ;
               8,  9,  2 ;
               7,  2, 12 ;
               8,  3, 10 ;
               7, 11,  3 ];
    
    otherwise
        error('Unknown type of mesh.');
end
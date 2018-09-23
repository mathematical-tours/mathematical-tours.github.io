%%
% Just a 4x4 orthobasis test.


% 2ab + b2 + (a+2b)^2 = 0
% a2 + 5b2 + 6ab = 0
% 1 + 6b + 5b2=0  ==>  Delta=36-20 = 16
% b = (-6 + 4)/2 = -1

a = 1; b = -1;
c = -a-2*b;
U = [1, 1, 1, 1; ...
    a, b, b, c; ...
    b, a, b, c; ...
    b, b, a, c];
 
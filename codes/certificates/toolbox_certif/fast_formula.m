function v = fast_formula(f,X,Y)

% fast_formula - fast evaluation of a symbolic function 
%
%   v = fast_formula(f,X,Y);
%
%   f should be a symbolic function of two variable f(x,y)
%
%   return v=f(X,Y) where xx and yy are usual vectors or matrices.
%   
%   Replaces *,/,^ by vectorized .*,./,.^
%
%   Copyright (c) 2017 Gabriel Peyre

syms xx yy;
A = char(f(xx,yy));
A = strrep(A, '*', '.*');
A = strrep(A, '/', './');
A = strrep(A, '^', '.^');
A = strrep(A, 'xx', 'X');
A = strrep(A, 'yy', 'Y');
eval(['v=' A ';']);

end
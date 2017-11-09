% P: number of nodes to generate randomly
% N: number of derivatives to expand up to
%
% M= M_fac{1}*M_fac{2}*M_fac{3} is a matrix describing the Taylor expansion, up to order 'N'
% of a kernel at each node point in 'Nodes'.
%
% C = C_sep{1}+C_sep{2} says which partial derivative each column of M corresponds to.
% For example, if
% C = [1 y x y^2 xy+t*x^3 x^2];
% then M(:,1) corresponds to coefficient in front of f(0,0)
% M(:,2) corresponds to coefficient in front of f_{y}(0,0)
% M(:,3) corresponds to coefficient in front of f_{x}(0,0)
% M(:,4) corresponds to coefficient in front of f_{yy}(0,0)
% M(:,5) corresponds to coefficient in front of f_{xy}(0,0) + t*f_{xxx}(0,0) ...
%
% The contribution of C_sep{2} is insignificant as t->0.
%
% The rows are ordered in terms of derivatives, then nodes: 
% for example, if Node 1 = N1= (n1,m1) and Node 2 = N2 = (n2,m2),
% then row 1 --> f(N1), row 2 --> f(N2), row 3--> f_x(N1)
% row 4 --> f_x(N2), row 5 --> f_y(N1), row 6--> f_y(N2)

function [M_fac,C_sep, Nodes,max_deg] = TaylorMtx(X,N)

% Nodes is (P,2)

P = size(X,1);

syms t u v

syms x y

% theta = pi/14;
% RotM = @(p) [cos(p) sin(p);-sin(p) cos(p)];
% Nodes = RotM(theta)*Nodes';
% Nodes = Nodes';

% Nodes = sym([randi([-100 100],P,1),randi([-100 100],P,1)]);

Nodes = sym(round(X));
Nodes = t*Nodes; % Nodes(1,:) = sym([0 0]);

% Nodes = [0 0;  1 0;u v; 0 1];
len_n = length(Nodes(:,1));

%Create the Taylor expansion matrix
M=[]; i = 1; colInfo=sym('A', [1 (N+2)*(N+1)/2]);

M0 = []; M1 = []; M2 = [];
for p = 1:len_n
    r0 = [];
    r11 = [];
    r12 = [];
    for n=0:N
        for j=0:n
            k= n-j;
            if p==1
                %                         colInfo{i} = sprintf('C%d%d',j,k); i =i+1;
                colInfo(i) = x^j*y^k; i=i+1;
            end
            r0 = [r0 taylorcoeff(Nodes(p,1), Nodes(p,2), j,k)];
            r11 = [r11 taylorcoeff(Nodes(p,1), Nodes(p,2), j-1,k)];
            r12 = [r12 taylorcoeff(Nodes(p,1), Nodes(p,2), j,k-1)];
            
        end
    end
%     M = [M;r0;r11;r12];
M0 = [M0;r0];
M1 = [M1;r11];
M2 = [M2;r12];
end
M = [M0;M1;M2];
M_taylor = M;

%perform column echelon
for j=1:length(M(1,:))-1
    
    %stop if remaining columns all zeros
    status = char(sum(sum(abs(M(:,j:end)))));
    if length(status)>1 || str2num(status)~=0
        
        %move zero columns to the end
        while str2num(char(sum(abs(M(:,j)))))== 0
            M = movetoend(M, j);
        end
        
        if M(j,j)==0
            l = findnonzero(M,j);
            fprintf('Move col %d to col %d \n', l,j);
            M = movecol(M,l,j);
            
        end
        
        if M(j,j)~= 0
            for N = j+1:length(M(1,:))
                if M(j,N)~=0
                    d = M(j,N)/M(j,j);
                    M(:,N) = M(:,N) - d*M(:,j);
                end
            end
        end
    else
        break
    end
    
end

M = simplify(M);

L = length(M(:,1));

%If assertion fails, then you should have chosen larger N to be larger
assert(str2double(char(sum(sum(abs(M(:,L+1:end))))))==0);


M = M(:,1:L);
C = M\M_taylor*colInfo.';

%factorize M
t=0;
C2 = subs(C); 
C_sep = {C2; C-C2};
syms t
D3 = diag([sym(ones(len_n,1)); sym(ones(2*len_n,1))/t]);

x=rand(1)*t; y=rand(1)*t;
C2 = subs(C2);
t=1;
coeffs = subs(C2); C2 = C2./coeffs;

D2 = diag(C2);
T = D3\M/D2;
M_fac = {D3,T,D2};

t=2;
max_deg = max(double(log2(subs(diag(D2)))));

end

%find the first nonzero entry in row j after index j
function l = findnonzero(M,j)
for l=j+1:length(M(1,:))
    if M(j,l)~=0
        
        break;
    end
end
end

%remove kth column and place it above column j
function Mnew = movecol(M, k,j)
assert(k>j);
r = M(:,k);
Mnew = M;
M(:,k) = [];
Mnew(:,j+1:end) = M(:,j:end);
Mnew(:,j) = r;



end

%make kth column the last column
function M = movetoend(M,k)
c = M(:,k);
M(:,k) = [];
M = [M, c];

end

%get taylor coefficient of order (n,m), centred at (0,0) in direction (u,v)
function coeff = taylorcoeff(u,v, n,m)
if n>=0
    un = u^n;
    fn = factorial(n);
else
    un = 0;  fn = 1;
end

if m>=0
    vm = v^m;  fm = factorial(m);
else
    vm = 0;  fm = 1;
end

coeff = un*vm/fn/fm;
end

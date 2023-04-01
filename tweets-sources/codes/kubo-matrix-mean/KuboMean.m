% geometric mean : A # B = A #_{1/2} B
% geodesic : A #_t B := A^{1/2} ( A^{-1/2} B A^{-1/2} )^t A^{1/2}
% solution of Ricatti equation : X A^{-1} X = B <=> X=A#B

% Unique mean satisfying:
% Congruence invariance : T^* (A # B) T = (T^* A T) # (T^* B T)
% operator monotone : A<=A', B<=B' ==> A#B <= A'#B'

% symmetry: A#B=B#A
% inversion : A^{-1} # B^{-1} = (A#B)^{-1}
% operator concave : ((1-t)A+tA')#((1-t)B+tB') >= (1-t)A#B + tA'#B'
% Arithmetic–geometric–harmonic mean inequality: : 
%% (A^{-1}/2+B^{-1}/2)^{-1} <= A#B <= (A+B)/2
% geodesic for d(A,B) = tr( log(A^{-1}B)^2 )^{1/2} = | A^{-1/2}BA^{-1/2} |_2
% Fumio Kubo, Tsuyoshi Ando , Means of positive linear operators, Math. Ann., 1980
% http://math.bme.hu/~petz/nagymedal.html

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

% generate symmetric matrices
R = @(t)[cos(t), sin(t); -sin(t), cos(t)];
S = @(t,a)R(t)*diag([1,a])*R(-t);

%
A{1} = S( .5 * pi/2,.03);
A{2} = S(-.5 * pi/2,.25);
%
A{1} = S(.7 * pi/2,.03);
A{2} = S(0,1);
q = 50;
xmax = 10;
for it=1:q
    clf; hold on;
    for jt=1:it
        t = (jt-1)/(q-1);
        col = [1-t 0 t];
        % linear mean
        B = (1-t)*A{1} + t*A{2};
        % kubo 
        Kubo = @(X,Y,t) real( X^(1/2) * (X^(-1/2) * Y * X^(-1/2) )^t * X^(1/2) );
        B = Kubo(A{1},A{2},t);
        %
        m = t*xmax;
        DrawEllipse(B,m,col,.2);
    end
    DrawEllipse(B,m,col,1);
    axis([-1 xmax+1, -1 1]);
    axis equal; axis off;
    drawnow;
    mysaveas(it);
end









axis equal;
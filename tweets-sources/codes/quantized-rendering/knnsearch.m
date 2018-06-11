function [idx,D]=knnsearch(varargin)
% KNNSEARCH   Linear k-nearest neighbor (KNN) search
% IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by eahc row of Q (m x d array).
% The results are stored in the (m x K) index array, IDX. 
%
% IDX = knnsearch(Q,R) takes the default value K=1.
%
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.
%
% Rationality
% Linear KNN search is the simplest appraoch of KNN. The search is based on
% calculation of all distances. Therefore, it is normally believed only
% suitable for small data sets. However, other advanced approaches, such as
% kd-tree and delaunary become inefficient when d is large comparing to the
% number of data points. On the other hand, the linear search in MATLAB is
% relatively insensitive to d due to the vectorization. In  this code, the 
% efficiency of linear search is further improved by using the JIT
% aceeleration of MATLAB. Numerical example shows that its performance is
% comparable with kd-tree algorithm in mex.
%
% See also, kdtree, nnsearch, delaunary, dsearch

% By Yi Cao at Cranfield University on 25 March 2008

% Example 1: small data sets
%{
R=randn(100,2);
Q=randn(3,2);
idx=knnsearch(Q,R);
plot(R(:,1),R(:,2),'b.',Q(:,1),Q(:,2),'ro',R(idx,1),R(idx,2),'gx');
%}

% Example 2: ten nearest points to [0 0]
%{
R=rand(100,2);
Q=[0 0];
K=10;
idx=knnsearch(Q,R,10);
r=max(sqrt(sum(R(idx,:).^2,2)));
theta=0:0.01:pi/2;
x=r*cos(theta);
y=r*sin(theta);
plot(R(:,1),R(:,2),'b.',Q(:,1),Q(:,2),'co',R(idx,1),R(idx,2),'gx',x,y,'r-','linewidth',2);
%}

% Example 3: cputime comparion with delaunay+dsearch I, a few to look up
%{
R=randn(10000,4);
Q=randn(500,4);
t0=cputime;
idx=knnsearch(Q,R);
t1=cputime;
T=delaunayn(R);
idx1=dsearchn(R,T,Q);
t2=cputime;
fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
fprintf('CPU time for knnsearch = %g\n',t1-t0);
fprintf('CPU time for delaunay  = %g\n',t2-t1);
%}
% Example 4: cputime comparion with delaunay+dsearch II, lots to look up
%{
Q=randn(10000,4);
R=randn(500,4);
t0=cputime;
idx=knnsearch(Q,R);
t1=cputime;
T=delaunayn(R);
idx1=dsearchn(R,T,Q);
t2=cputime;
fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
fprintf('CPU time for knnsearch = %g\n',t1-t0);
fprintf('CPU time for delaunay  = %g\n',t2-t1);
%}
% Example 5: cputime comparion with kd-tree by Steven Michael (mex file) 
% <a href="http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7030&objectType=file">kd-tree by Steven Michael</a> 
%{
Q=randn(10000,10);
R=randn(500,10);
t0=cputime;
idx=knnsearch(Q,R);
t1=cputime;
tree=kdtree(R);
idx1=kdtree_closestpoint(tree,Q);
t2=cputime;
fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
fprintf('CPU time for knnsearch = %g\n',t1-t0);
fprintf('CPU time for delaunay  = %g\n',t2-t1);
%}

% Check inputs
[Q,R,K,fident] = parseinputs(varargin{:});

% Check outputs
error(nargoutchk(0,2,nargout));

% C2 = sum(C.*C,2)';
[N,M] = size(Q);
L=size(R,1);
idx = zeros(N,K);
D = idx;

if K==1
    % Loop for each query point
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [D(k),idx(k)]=min(d);
    end
else
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [s,t]=sort(d);
        idx(k,:)=t(1:K);
        D(k,:)=s(1:K);
    end
end
if nargout>1
    D=sqrt(D);
end

function [Q,R,K,fident] = parseinputs(varargin)
% Check input and output
error(nargchk(1,3,nargin));

Q=varargin{1};

if nargin<2
    R=Q;
    fident = true;
else
    fident = false;
    R=varargin{2};
end

if isempty(R)
    fident = true;
    R=Q;
end

if ~fident
    fident = isequal(Q,R);
end

if nargin<3
    K=1;
else
    K=varargin{3};
end

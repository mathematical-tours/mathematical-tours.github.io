function GameOfLife
% This is a simple simulation of Conway Game of life GoL
% it is good for understanding Cellular Automata (CA) concept
% GoL Rules:
%   1. Survival: an alive cell live if it has 2 or 3 alive neighbors
%   2. Birth: a dead cell will be alive if it has 3 alive neighbors
%   3. Deaths: 
%        Lonless: alive cell dies if it has 0 or 1 alive neighbors    
%        Overcrowding: alive cell dies if it has 4 or more alive neighbors    
% Any questions related to CA are welcome
% By: Ibraheem Al-Dhamari

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size= 500x500 
% different random initial values
%A= rand(500,500);
%A= ones(500,500);

% periodic configuration
%A= zeros(500,500);
% A(100,200)=1;
% A(100,201)=1;
% A(100,202)=1;

% initial from image
 A=imread('CA02.JPG');

% convert to binary--> % states={0,1}
A=im2bw(A);
% visualize the initial states 
disp('the binary image')
imshow(A);
% pause
% boundary type: 0= reflection 
%                1= doublication
%                2= null, zeros

% this step enlarge A with 4 virtual vectors
A=Bnd(A,0); 
[d1,d2]=size(A);
% disp('the extended binary image')
% whos A
imshow(A);
pause
B=A;
t=0;
stp=false; % to stop when if no new configrations
%B is the CA in time t
%A is the CA in time t+1
%t is the number of generations 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Play ^_^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~stp & (t<10) % repeat for 10 generations
    % for each cell in the CA 
    for i=2:d1-1
        for j=2:d2-1           
            % apply Game of life rule     
            A(i,j)=GOL(A,B,i,j)            ;
        end
    end 
    % visualize what happened
    disp('the CA image')
    imshow(A);    
    drawnow;
%    pause
    % save B 
    if A==B
       stp=true; % no more new states
    end
    B=A;  
    t=t+1  
end  

%==========================================
%   Game of Life Rules 
%==========================================
function s=GOL (A,B,i,j)
% game of life rule
sm=0;
% count number of alive neighbors
sm=sm+ B(i-1,j-1)+B(i-1,j)+B(i-1,j+1);
sm=sm+ B(i,j-1)+           B(i,j+1);
sm=sm+ B(i+1,j-1)+B(i+1,j)+B(i+1,j+1);

% compute the new state of the current cell
s=B(i,j);
if B(i,j)==1
    if (sm>1)&&(sm<4)
        s=1;
    else
        s=0 ;   
    end
else
    if sm==3
       s=1;
    end
end
    
%==========================================
%   Boundary Type
%==========================================
function bA= Bnd(A,k)
% add new four vectors based on boundary type
[d1, d2]=size(A);
d1=d1+2; d2=d2+2;
X=ones(d1,d2);
X=im2bw(X);
X(2:d1-1,2:d2-1)=A;
imshow(X);
whos A X
if k==0 % Reflection
   X(  1  , 2:d2-1)=A(end , :);
   X(  d1 , 2:d2-1)=A( 1  , :);
   X( 2:d1-1 , 1  )=A(: , end);
   X( 2:d1-1 , d2 )=A(: ,  1 );
   
   X(1,1)    =A(end,end);
   X(1,end)  =A(end,1);
   X(end,1)  =A(1,end);
   X(end,end)=A(1,1);    
elseif k==1 % Double
   X(  1  , 2:d2-1)=A( 1  , :);
   X(  d1 , 2:d2-1)=A(end , :);
   X( 2:d1-1 , 1  )=A(: ,  1 );
   X( 2:d1-1 , d2 )=A(: , end);
   
   X(1,1)    =A(end,1);
   X(1,end)  =A(end,end);
   X(end,1)  =A(1,1);
   X(end,end)=A(1,end);

else % k==2 % zeros
   X(  1  ,:)=0;
   X( end ,:)=0;
   X(: ,  1  )=0;
   X(: , end )=0;
end
bA=X;


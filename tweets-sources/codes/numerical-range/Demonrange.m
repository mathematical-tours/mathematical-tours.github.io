% Demo for the nrange function
% This program uses the nrange function to illustrate
% the numerical ranges of various square matrices.

% Reference: The procedure nrange we wrote is
% repetedly used. See the m-file nrange.m for
% this procedure and the way to use it. This demo
% will not function without it.

% Copyright:
% Dora Matache and Valentin Matache
% University of Nebraska Omaha
% Department of Mathematics
% Year 2001

disp('  ');
disp('The following is a presentation of the');
disp('routine nrange.');
disp('It graphs the numerical range of any square');
disp('matrix with complex entries');
disp('and calculates its numerical radius.')
disp('The standard command is nrange(A) where A');
disp('is a previously introduced matrix.');
disp('To see a first example press enter.');
disp('   ')
pause;
disp('Numerical range of ')
A=[1 1 1; 0 1 1; 0 0 1]
nrange(A);
disp('2x2 matrices have only')
disp('elliptic numerical ranges,');
disp('posibly degenerate, (circular disks, segments, and poits).');
disp('In 3x3 matrices, like this one, the numerical');
disp('range can have various shapes.');
disp('In the following we will illustrate this');
disp('diversity in several examples.');
disp('You can stop at any time by pressing control+c');
disp('Press enter to continue.');
disp(' ')
pause;

disp('Numerical range of')
A=[1+i 0 0; 0 1/2 1; 0 0 -1/2]
nrange(A,500,'m');
disp('Press Enter to continue');
pause;

disp('Numerical range of')
A=[1 0 0; i/3 i/2 0; -1/9 -1/3 -1/4]
nrange(A,500,'g');

disp('Press Enter to continue');
pause;

disp('Numerical range of')
A=[2 4 9; 6 8 1; 7 3 5]
nrange(A,500,'c');
disp('Press Enter to continue');
pause;

disp('Numerical range of')
A=[0 1 0; 0 0 1; 1/2 0 0]
nrange(A,500,'r');
disp('Press Enter to continue');
pause;

disp('Numerical range of')
A=[1 i i+1; 1/2 0 1; i+2 0 0]
nrange(A,500,'g');
disp('Press Enter to continue');
pause;


disp('Numerical range of')
A=[1 0 0;0 i 0; 0 0 -i]
nrange(A,500,'k');
disp('One can graph only the boundary.');
disp('Let us use the matrix above to illustrate this feature.')
disp('Press Enter to continue');
pause;

disp('Boundary of the numerical range of')
A=[1 0 0;0 i 0; 0 0 -i]
nrange(A,500,'k','n');
disp('The larger the dimension, the more diverse');
disp('the numerical range gets.');
disp('Here is a 6x6 matrix.');
disp('Press enter to continue.');
pause;

disp('Numerical range of')
A=[4 0 0 0 0 0; i/3 2+i/2 0 0 0 0; -1/9 -1/3 1.75 0 0 0; 0 0 0 0 1 -1; 0 0 0 -1 0 1; 0 0 0 -1 -1 0 ]
nrange(A);
disp('Here is a beautiful example: another 6x6 matrix.');
disp('Press enter to continue');
pause;

disp('Numerical range of')
A=[0 0 0 1 0 0; 0 0 0 0 0 0; 1 0 0 0 0 0; 0 0 1 0 0 0; 0 1 0 0 0 0; 0 0 0 0 1 0]
nrange(A,500,'m'); 

disp('Ready for a little motion?');
disp('Here is a "movie" where the numerical');
disp('ranges of 3x3 matrices randomly generated');
disp('with the command rand(3)+i*rand(3)');
disp('are graphed consecutively.');
disp('Press enter to play it.');
pause;

for j=1:10
    A=rand(3)+i*rand(3);
    nrange(A,500,'m');
    fmat(:,j)=getframe;
    A=rand(3)+i*rand(3);
    nrange(A,500,'c');
    fmat(:,j)=getframe;
    A=rand(3)+i*rand(3);
    nrange(A,500,'g');
    fmat(:,j)=getframe;
    A=rand(3)+i*rand(3);
    nrange(A,500,'r');
    fmat(:,j)=getframe;A=rand(3)+i*rand(3);
    nrange(A,500,'k');
    fmat(:,j)=getframe;
    A=rand(3)+i*rand(3);
    nrange(A,500,'r');
    fmat(:,j)=getframe;
    A=rand(3)+i*rand(3);
    nrange(A,500,'y');
    fmat(:,j)=getframe;
    A=rand(3)+i*rand(3);
    nrange(A,500,'b');
    fmat(:,j)=getframe;
end

movie(fmat,1);
    
disp('THE END')
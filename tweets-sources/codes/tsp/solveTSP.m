function varargout = solveTSP( cities, display, maxIteration, order)
% cities = solveTSP( cities, maxItt, display)
%
% cities - An Nx2 matrix containing cartesian coordinates of the "cities"
% beeing visited. The initial trail is assumed from the first city to the
% scond and so on...
%
% display - bolean flag decide if to display the progress of the program (slows the running time).
%                       default = false;
%
% maxIteration - maximum iterations for the program
%                                   default = 10,000
%
%
% [cities ind] = solveTSP( cities, display) returns the aranged cities and
% an index vector of the visiting order
%
% [cities ind totalDist] = solveTSP( cities, display)
% totalDist is the route total distance
%
% demo1:
% cities = solveTSP( rand(100,2), true );
%
% demo2:
% t = (0:999)' /1000;
% cities = [ t.^2.*cos( t*30 ) t.^2.*sin( t*30 ) ];
% [ans ind] = sort( rand(1000,1) );
% [cities ind] = solveTSP( cities(ind,:), true );

if nargin < 2
    display = false;
end
if nargin<3
%maxIteration = 1e5;
maxIteration = 1000;
end

siz = size(cities);
if siz(2) ~= 2
    error( 'The program is expecting cities to be an Nx2 matix of cartesian coordinates' );
end
N = siz(1);

if nargin<4
    order = (1:N)';     % initial cities visit order
end

if display
    hFig = gcf; % figure;
    hAx = gca;
    updateRate = ceil( N/50 );
end


itt = 1;
maxItt = min(20*N,maxIteration);
noChange = 0;

while itt < maxItt  && noChange < N
    
    dist = calcDistVec( cities(order,:),1 );    % travel distance between the cities
    
    %% ----------- Displaying current route -----------------------
    if display && ~mod(itt,updateRate) && ishandle( hFig )
        hold(hAx,'off');
        plot( hAx, cities( order,1),cities( order,2),'r.' );
        hold( hAx,'on');
        plot( hAx, cities( order,1),cities( order,2) );
        str = {[ 'iteration: ' num2str( itt ) ] ;
            [ 'total route: ' num2str( sum( dist) ) ] };
        title(  hAx,str );
        pause(0.02)
    end
    
    flip = mod( itt-1, N-3 )+2 ;
    
    untie = dist(1:end-flip) + dist(flip+1:end);    % the distance saved by untying a loop
    shufledDist = calcDistVec( cities( order,:),flip );
    connect =  shufledDist(1:end-1) + shufledDist( 2:end);  % the distance payed by connecting the loop (after flip)
    benifit = connect - untie;  % "what's the distance benifit from this loop fliping
    
    %% --------------- Finding the optimal flips (most benficial) ----------------
    localMin = imerode(benifit,ones(2*flip+1,1) );
    minimasInd = find( localMin == benifit);
    reqFlips = minimasInd( benifit(minimasInd) < -eps );
    
    %% -------- fliping all loops found worth fliping --------------------
    prevOrd = order;
    for n=1:numel( reqFlips )
        order( reqFlips(n) : reqFlips(n)+flip-1 ) = order( reqFlips(n) +flip-1: -1 :reqFlips(n) );
    end
    
    %% -------  counting how many iterations there was no improvement
    if isequal( order,prevOrd )
        noChange = noChange + 1;
    else
        noChange = 0;
    end
    
    itt = itt+1;
    
end     % while itt < maxItt  && noChange < N

output = {cities( order,:), order, sum( dist)};
varargout = output(1:nargout);

function dist = calcDistVec( cord,offset )
% dist = calcDistVec( cord,offset )
% offset is the number of cities to calculate the distence between
% the distance for the first city is allway 0

dist = zeros( size(cord,1)-offset+1,1 );
temp = cord( 1:end-offset,:) - cord( offset+1:end,:);
dist(2:end) = sqrt( sum(temp.^2,2) );
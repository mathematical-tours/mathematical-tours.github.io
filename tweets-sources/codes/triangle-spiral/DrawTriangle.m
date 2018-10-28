function DrawTriangle(T, delta, col, lw)

if size(T,1)>1
    hold on; axis equal;
    for i=1:size(T,1)
        Tr = T(i,:);
        DrawTriangle(Tr, delta, col, lw);
    end
    axis off;
    return;
end

if iscell(T)
    hold on; axis equal;
    for i=1:length(T)
        DrawTriangle(T{i}, delta, col, lw);
    end
    axis off;
    return;
end

hold on;
while true
    plot(T([1:end 1]), 'Color', col, 'LineWidth', lw);
    % decide next move
    L = mean( abs( T-T([2:end 1]) ) );
    if L<delta
        break;
    end
    t = min(delta/L,1);%
%    T = (1-t)*T + t*T([2:end 1]);
    T1 = T;
    for i=1:length(T1)
        j = mod(i,length(T1))+1;
        l = abs(T1(i)-T1(j)); t = min(delta/l,.95);
        T(i) = (1-t)*T1(i) + t*T1(j);
    end
end

end

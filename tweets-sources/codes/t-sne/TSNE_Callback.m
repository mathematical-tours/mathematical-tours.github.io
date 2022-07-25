function stop = TSNE_Callback(optimValues,state,species)

global param_t
global disp_list
global cnt

stop = false;

switch state
    case 'init'
    case 'iter'
        if optimValues.iteration==disp_list(cnt)
            Y = optimValues.Y;
            clf;
            scatter(Y(:,1),Y(:,2), 20, param_t(:,1), 'filled'); % 
            axis equal;
            axis off;
            a = mean( sqrt( sum(Y.^2,2) ) );
            axis([-1 1 -1 1]*a*2);
            drawnow;
            global rep
            saveas(gcf, [rep 'anim-' znum2str(cnt,3) '.png']);
            cnt = cnt+1;
        end
    case 'done'
end

end
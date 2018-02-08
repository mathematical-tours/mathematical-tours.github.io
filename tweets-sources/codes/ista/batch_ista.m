%% 
% batch run ista

eta_list = [0 .25 .5 1 2];

for ie=1:length(eta_list)
    eta = eta_list(ie);
    PlotIsta;
end
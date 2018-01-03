%%
% Plot Daubechies Wavelets

rep = '../results/daubechies/';
[~,~] = mkdir(rep);

addpath('../toolbox/');

N = 1024;
j = 3;

vm_list = linspace(1,8,50);

for i=1:length(vm_list)
    vm = vm_list(i);
    h1 = compute_wavelet_filter('Daubechies',floor(vm)*2); h1(end+1:end+2) = 0;
    h2 = compute_wavelet_filter('Daubechies',(floor(vm)+1)*2);
    u = vm - floor(vm);
    options.h = h1*(1-u)+h2*u;
    fw = zeros(N,1);
    fw(1 + 2^j) = 1;
    f = perform_wavortho_transf(fw,1,-1,options);
    f = circshift(f,N/2);
    % display
    t = linspace(0,1,N);
    a = [-.15 .2];
    s = (i-1)/(length(vm_list)-1);
    clf;
    plot(t, f/max(abs(f)), 'Color', [s 0 1-s], 'LineWidth', 2);
    axis([.1 .9 -1.02 1.02]); axis off;
    drawnow;
    saveas(gcf, [rep 'daubechies-' num2str(i), '.png']);
end
%%
% Display Cepstrum

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20,'FontName','Helvetica');

% y, and a sample rate for that data, Fs
name = 'cello';
[y,Fs] = audioread([name '.wav']);
t = linspace(0,1,length(y));
% sound(y,Fs);


tb = [.23 .27];
tb = [.22 .28];
I = find(t>tb(1) & t<tb(2));
f = y(I);

n = length(f);
% duration: 
T = n/Fs;

clf;
plot(linspace(0,T,n),f, 'k');
axis([0,T,-max(abs(f)), max(abs(f))]);
SetAR(1/2);
saveas(gcf, [rep name '-sound.png']);


dB = @(u)20*log10(u/max(u));
Lf = dB(abs(fft(f)));
plot(linspace(0,2*pi/T*n/2,n/2) / 1e3, Lf(1:end/2), 'b');
axis([0,2*pi/T*n/2/1e3,-120,0]);
SetAR(1/2);
saveas(gcf, [rep name '-frequency.png']);

% quefrency 
Sf = real( ifft( log(abs(fft(f))) ) );
plot(abs(Sf(1:end/2)), 'r');
SetAR(1/2);
axis([0 5000 0 0.04]);
saveas(gcf, [rep name '-cepstrum.png']);

SetAR(1/2);
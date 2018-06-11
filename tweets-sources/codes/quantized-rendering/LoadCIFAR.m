function Z = LoadCIFAR()

% change your user name if you are a mac user
username = 'gpeyre';

CIFAR = 100;
% change the path CIFAR database here
proot = ['/Users/' username '/Documents/MATLAB/data/'];
if CIFAR==10
    a = [];
    for k=1:5
        load([proot 'cifar-10/data_batch_' num2str(k)]);
        a = [a; data];
    end
elseif CIFAR==100
    load([proot 'cifar-100-matlab/train']);
    a = data;
    load([proot 'cifar-100-matlab/test']);
    a = [a; data];
end
Z = reshape(a', [32 32 3 size(a,1)]);
Z = double(Z);
clear data; clear a;

end
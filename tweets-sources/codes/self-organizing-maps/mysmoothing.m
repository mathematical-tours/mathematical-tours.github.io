function H = mysmoothing(Z, niter)

% smooth out
H = Z;
for k=1:niter
    H = H + H([1 1:end-1],:) + H([2:end end],:);
    H = H + H(:,[1 1:end-1]) + H(:,[2:end end]);
    H = H/9;
end

end
function f = mybil(f,sigmar, sigmas)

eps    = 1e-3;
for k=1:size(f,3)
    f(:,:,k) = GPA(f(:,:,k)*255, sigmar, sigmas, eps, 'Gauss')/255;
end

end
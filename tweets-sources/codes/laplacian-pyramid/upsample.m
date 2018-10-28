function fjU = upsample(fj)

% up-sampling by linear interpolation
fjU = zeros(2*size(fj));
fjU(1:2:end,1:2:end) = fj;
%
fjU(2:2:end,1:2:end) = (fj + fj([2:end 1],:))/2;
fjU(1:2:end,2:2:end) = (fj + fj(:,[2:end 1]))/2;
fjU(2:2:end,2:2:end) = (fj + fj([2:end 1],:) + fj(:,[2:end 1]) + + fj([2:end 1],[2:end 1]))/4;
% details


end
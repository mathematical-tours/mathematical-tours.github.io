function b = tensor_mult(K1,a)

% tensor_mult - multiply tensor by vector

b = zeros(size(K1,1),size(K1,3));
for i=1:size(K1,3)
    b(:,i) = K1(:,:,i)*a;
end

end
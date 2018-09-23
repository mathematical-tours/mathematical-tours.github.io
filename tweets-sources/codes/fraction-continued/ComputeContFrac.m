function a = ComputeContFrac(x,n)

a = [];
for i=1:n
    a(i) = floor(x);
    x = 1/( x-a(i) );
end

end
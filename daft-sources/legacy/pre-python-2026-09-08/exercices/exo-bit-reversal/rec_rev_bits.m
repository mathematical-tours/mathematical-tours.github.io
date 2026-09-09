function y = rec_rev_bits(x)
n = length(x);
if n>1
   y = [rec_rev_bits(x(1:2:n)) ; rec_rev_bits(x(2:2:n))]; 
else
    y = x;
end
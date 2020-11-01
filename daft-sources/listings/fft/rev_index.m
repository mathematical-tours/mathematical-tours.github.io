% Calcule l'inversion des bits d'un entier.
function res = rev_index(t,index)
res = 0;
tmp = index;
for i=0:t-1
    bit = mod(tmp,2);
    tmp = floor(tmp/2);
    res = res*2 + bit;
end
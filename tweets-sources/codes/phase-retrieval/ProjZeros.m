function g = ProjZeros(g,sel)

if 0
g1 = g;
g = g*0;
g(sel,sel) = g1(sel,sel);
end
g(sel) = 0;

end
function x = ProjDisk(x,c,r)

if abs(x-c)>r
    x = c + r * (x-c)/abs(x-c);
end

end

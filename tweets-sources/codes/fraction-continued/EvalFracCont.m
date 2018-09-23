function [x] = EvalFracCont(a)

if length(a)==1
    x = a;
else
    x = a(1) + 1/EvalFracCont(a(2:end));
end

end
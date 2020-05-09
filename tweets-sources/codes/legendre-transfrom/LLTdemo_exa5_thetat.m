function thetat=LLTdemo_exa5_thetat(x)
%Compute abs(x)<=t, thetat=t*(1-sqrt(1-(x/t)^2)); else thetat=inf;
% x is a real number. t is a global parameter and a real number.
global t;
if abs(x)<=t,
	thetat=t*(1-sqrt(1-(x/t)^2));
else
	thetat=inf;
end




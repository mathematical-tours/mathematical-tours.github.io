function y = transformee_sinus_2d(x)
y = zeros(size(x)); n = length(x);
for(i=1:n) y(i,:) = transformee_sinus(x(i,:)')'; end;
for(j=1:n) y(:,j) = transformee_sinus(y(:,j)); end;
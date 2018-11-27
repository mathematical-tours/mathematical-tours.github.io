function P = randsubsample(P,n)

P = P(:);
I = randperm(length(P)); 
I = I(1:n);
P = P(I);

end
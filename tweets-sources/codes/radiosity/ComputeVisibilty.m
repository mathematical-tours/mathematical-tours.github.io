function V = ComputeVisibilty(oc,r, P, R)

% Oc: occluders centers
% r: occluders radius


dotp = @(u,v)sum(u.*v,3);

p = size(P,1);

Pz = repmat(P,[1 p]);
V = ones(p,p);  % visiblity 
for k=1:length(oc)
    % distance btween point Oc and  line passing through Pz and directed by
    % R
    Oc = repmat(reshape(oc{k}, [1 1 3]), [p p]);
    dProj = dotp(Oc-Pz,R);
    d1 = sum( (Oc-Pz).^2, 3);
    d2 = d1-dProj.^2;
    %
    V = V .* (d2>=r(k)^2);    
end

end
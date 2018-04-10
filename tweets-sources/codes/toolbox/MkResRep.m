function rep = MkResRep(name)

if nargin==0
    name = [];
end

s = pwd(); a = strfind(s, '/'); s = s(a(end)+1:end);
rep = ['../results/' s '/'];
if not(isempty(name))
    rep = [rep name '/'];
end
[~,~] = mkdir(rep);

end
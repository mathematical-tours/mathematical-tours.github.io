function v = getoptions(options, name, v, mandatory)

% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0, mendatory);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    error('Not enough arguments.');
end
if nargin<4
    mandatory = 0;
end

if isfield(options, name)
    v = options.(name);
elseif mandatory
    error(['You have to provide options.' name '.']);
end 
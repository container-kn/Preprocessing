function [X,Zf] = moving_average(N,X,Zi)
% Run a moving-average filter along the second dimension of the data.
% [X,Zf] = moving_average(N,X,Zi)
%
% In:
%   N : filter length in samples
%   X : data matrix [#Channels x #Samples]
%   Zi : initial filter conditions (default: [])
%
% Out:
%   X : the filtered data
%   Zf : final filter conditions
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2012-01-10

if nargin <= 2 || isempty(Zi)
    Zi = zeros(size(X,1),N); end

% pre-pend initial state & get dimensions
Y = [Zi X]; M = size(Y,2);
% get alternating index vector (for additions & subtractions)
I = [1:M-N; 1+N:M];
% get sign vector (also alternating, and includes the scaling)
S = [-ones(1,M-N); ones(1,M-N)]/N;
% run moving average
X = cumsum(bsxfun(@times,Y(:,I(:)),S(:)'),2);
% read out result
X = X(:,2:2:end);

if nargout > 1
    Zf = [-(X(:,end)*N-Y(:,end-N+1)) Y(:,end-N+2:end)]; end
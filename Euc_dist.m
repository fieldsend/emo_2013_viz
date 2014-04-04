function d = Euc_dist(a,b)

% function d = Euc_dist(a,b)
%
%
% INPUTS
%
% a = matrix of points - n1 by m
% b = matrix of points - n2 by m
%
%
% returns the Euclidean distance between all n1 vectors in a and all n2 
% vectors in b, returns in d


[n1, m1] = size(a);
[n2, m2] = size(b);
if ( m1 ~= m2 )
    error('Cannot calculate distances between matrices with col numbers');
end

d = sqrt((ones(n2,1) * sum(a.^2,2)')' + ones(n1,1) * sum(b.^2,2)' - 2.*(a*(b')));
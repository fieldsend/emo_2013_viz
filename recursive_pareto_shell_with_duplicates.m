function [shell] = recursive_pareto_shell_with_duplicates(Y,index)

% function [shell] = recursive_pareto_shell(Y,index)
%
% function recusively computes Pareto shell membership
%
% Y = A n by m matrix of objectives, where m is the number of objectives
% and n is the number of points
% index = the integer value you wish to attribute to the estimated Pareto
% optimal shell (typically 0 or 1)
%
% shell = n by 1 array of corresponding shell membership values
%
% copes with duplicates
% assumes minimisation
%
% Jonathan Fieldsend, 2011, University of Exeter

[n,m] = size(Y);
S = zeros(n,1);
shell = zeros(n,1);

for i=1:n
    % get number of points that dominate Y
    S(i) = sum((sum(Y<=repmat(Y(i,:),n,1),2) == m) & (sum(Y<repmat(Y(i,:),n,1),2) > 0));
end
S=S+1;
shell(S==1)=index;

if sum(S)==n % if at last shell, terminate and chain back
    return
else
    shell(S~=1) = recursive_pareto_shell_with_duplicates(Y(S~=1,:),index+1);
end


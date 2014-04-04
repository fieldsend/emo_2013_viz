function dm = get_dom_between_matrices(A,B)

% A is np by nd matrix of fitness values
% B is is a npp by nd matrix of fitness values
%
% return a matrix npp by np which details which members of B
% are dominated by which members of A
% 1 in element ji indicates A member i dominates B member j
% 0 indicates member A i does not dominate member B j
%
%can halve this complexity by just calculating upper or lower
%triangle
%
%
% Copyright (c) Jonathan E. Fieldsend 2013
        
[np, nd] = size(A);
[npp, ndd] = size(B);
if (ndd~=nd)
    error('Number of objectives must be the same in both matrices');
end
dm = zeros(npp,np);

for i=1:np; %for each matrix A member
    dom = zeros(npp,1);
    for j=1:nd;
        dom = dom + ( A(i,j) <= B(:,j) );
    end
    dom = (dom==nd);
    dm(:,i) = dom;
end

end


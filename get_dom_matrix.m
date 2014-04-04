function dm = get_dom_matrix(A)

% A is np by nd matrix of fitness values
% return a matrix np by np
% 1 in element ji indicate member i dominates member j
% 0 indicates member i does not dominate member j
%
% Copyright (c) Jonathan E. Fieldsend 2013
        
%can halve this complexity by just calculating upper or lower
%triangle
    
    
[np, nd] = size(A);
dm = zeros(np,np);

for i=1:np; %for each matrix member
    dom = zeros(np,1);
    for j=1:nd;
        dom = dom + ( A(i,j) <= A(:,j) );
    end
    dom = (dom==nd);
    dm(:,i) = dom;
end

end


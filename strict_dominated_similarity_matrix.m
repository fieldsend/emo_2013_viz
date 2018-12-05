function [ S ] = strict_dominated_similarity_matrix( R )

% function [ S ] = strict_dominated_similarity_matrix( R )
%
% 
% Copyright (c) Jonathan E. Fieldsend 2013
    
[np, nd] = size(R);
S = zeros(np,np);
for i=1:np
    for j=1:np
        di = zeros(np,1);
        dj = zeros(np,1);
        dboth = zeros(np,1);
        for m =1:nd
            dit = (R(i,m)>=R(:,m));
            djt = (R(j,m)>=R(:,m));
            dboth = dboth + (dit & djt);
            di = di + dit;
            dj = dj + djt;
        end
        dboth = (dboth==nd);
        di = (di==nd);
        dj = (dj==nd);
        S(i,j) = 2*sum(dboth)/(sum(di+dj));
    end
end


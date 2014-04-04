function [ S ] = dominance_similarity_matrix( R )

% function [ S ] = dominance_similarity_matrix( R )
%
% 
% Copyright (c) Jonathan E. Fieldsend 2013
    
[np, nd] = size(R);
S = zeros(np,np);
for i=1:np
    for j=1:np
        d = zeros(np,1);
        for m =1:nd
            d = d + ((R(i,m)<R(:,m)) & (R(j,m)<R(:,m)));
            d = d + ((R(i,m)==R(:,m)) & (R(j,m)==R(:,m)));
            d = d + ((R(i,m)>R(:,m)) & (R(j,m)>R(:,m)));
        end
        d = d/nd;
        S(i,j) = mean(d);
    end
end


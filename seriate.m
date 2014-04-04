function [ permutation ] = seriate( A )


% Will probably refactor this at version 2.0 as there are more
% efficient ways to compute this
%
% Copyright (c) Jonathan E. Fieldsend 2013
    

K = cumsum(A);
K = K(end,:);
D = diag(K);
L = D-A;

[Evec, Eval] = eig(L);
% will almost always be second eigen vector, but just in case.
temp = Eval>eps;
temp = sum(temp);
%I=find(temp>0);
%diag(Eval)'
%Fiedler_vector = Evec(:,I(1)); % as eig orders from smallest to largest, and smallest is zero, get smallest non-zero
Fiedler_vector = Evec(:,2);
[temp, permutation] = sort(Fiedler_vector);



end


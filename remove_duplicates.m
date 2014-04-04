function A = remove_duplicates(A)

% function A = remove_duplicates(A)
%
% Service method to remove any duplicate elements in the
% set of points A
%
% Copyright (c) Jonathan E. Fieldsend 2013
    
% Remove any duplicates
to_remove =[];
for i=1:size(A,1)
    eq = sum(sum(A==repmat(A(i,:),size(A,1),1),2)==size(A,2));
    if (eq>1)
        to_remove = [to_remove i];
    end
end
A(to_remove,:)=[];

disp(['removed ' int2str(length(to_remove)) ' duplicates']);
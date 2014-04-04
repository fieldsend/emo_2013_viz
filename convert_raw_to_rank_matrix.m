function [ R ] = convert_raw_to_rank_matrix( solution_set )

%function [ R ] = convert_raw_to_rank_matrix( solution_set )
%
% assumes minimisation when ranking raw data
% lowest rank is best
% will share ranks across multiple inidividuals (mid rank used, so not
% all ranks will necessarily be integers)
%
% Copyright (c) Jonathan E. Fieldsend 2013
    
[np, nd] = size(solution_set);
R=zeros(np,nd);

for i=1:nd
    [temp, I] = sort(solution_set(:,i));
    R(I(1),i)=1; %assign top rank 1
    sharing = 0;
    sharing_list =[];
    for j=2:np
        if (solution_set(I(j),i)~=solution_set(I(j-1),i)) % if no need to share ranks with previous
           R(I(j),i) = j;
           if (sharing) % imediately previous ranks were being shared, so calculate
              val = sharing/length(sharing_list);
              R(sharing_list,i) = val; 
              sharing =0; 
              sharing_list =[];
           end
        else  % need to share and keep track
            if sharing ==0 % first instance to share, so need to keep track
                sharing =j+j-1;
                sharing_list = [I(j-1) I(j)];
            else
               sharing = sharing+j;
               sharing_list = [sharing_list I(j)];
            end
        end
    end
    if sharing ~= 0 % end of column sharing to be fixed
        val = sharing/length(sharing_list);
        R(sharing_list,i) = val;
    end
end


end


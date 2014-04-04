function [F N] = get_fronts_via_dom_matrix(tdm)

% tdm = pre-computed domination matrix
% returns front (shell) membership for each individual
%
% Copyright (c) Jonathan E. Fieldsend 2013
        
shell=1;
n=length(tdm);
F=zeros(n,1);



while (sum(sum(tdm))>0)
    sum(sum(tdm));
    dc = sum(tdm,2); %get those individuals not dominated by anyone but themselves
    I = (dc == 1);
    F(I)=shell;
    %for i=1:length(I);
    tdm(:,I)=0;  %remove from dominating list
    %end
    shell=shell+1;
    
end

N = zeros(max(F),1);
for i=1:max(F)
    N(i) = sum(F==i); % get total numbers in each front
end

end
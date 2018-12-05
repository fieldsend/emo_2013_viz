function [A, S, R, F, N, p, shell0_locs, domination_matrix,p0angles] ...
    = seriate_pareto_front_and_convert_to_2Dshell( A, seriation_matrix_function, local_seriation)

% function [A, S, R, F, N, p, shell0_locs, domination_matrix,p0angles] = 
%       seriate_pareto_front_and_convert_to_2Dshell( A, seriation_type)
%
% Description
%
% General method to seriate a D-dimensional (D>=2) population of
% objective vectors, and map to 2-dimensional representation of the
% 0th shell (front)
%
% INPUTS
%
% A = set of objective vectors (columns objectives, rows elements)
% seriation_matrix_function = name of function to generate similarity matrix (in string)
% local_seriation = optional input, if not zero then local
%   seriation of just the shell zero elements is undertaken, rather
%   than seriation of the complete population. Default is not to have
%   local seriation.
%
% OUTPUTS
%
% A = original data with any duplicates removed
% S = matrix A converted into ranks
% R = Rank matrix equivalent of A
% F = front membership of each element of A
% N = number of elements in each front
% p = permutation of A based on spectral seriation (depending on
%   the third input argment, this may be a complete permutation, or
%   just a permutation of the shell zero elements)
% shell0_locs = mapping of shell 0 (i.e. nondominated) elements of
%   A onto the unit radius curve
% domination_matrix = matrix tracking which elements dominate which
%   others (1 in element ji indicate member i dominates member j
%   0 indicates member i does not dominate member j
% p0angles = angles of of shell 0 (i.e. nondominated) elements of
%   A when mapped to the curve in 2D
%
% Copyright (c) Jonathan E. Fieldsend 2013

if ( exist('local_seriation','var')==0 )
    local_seriation = 0;
end

A = remove_duplicates(A); % removes duplicates from A -- if you know
                         % there will not be duplicates in your
                         % code, can comment out

S = feval(seriation_matrix_function, A); %get similarity matrix
R  = convert_raw_to_rank_matrix( A );
domination_matrix = get_dom_matrix( A ); % keep track of which members dominate which
[F N] = get_fronts_via_dom_matrix( domination_matrix ); % get shell membership for elements of A

if (local_seriation)
    [p, original_D_dimension_distances] = local_approach(S, F, A, N);
else
    [p, original_D_dimension_distances] = global_approach(S, F, A, N);
end

% given order provided by seriation, now map non-dominated members
% to the unit radius curve
total_traversal = sum(original_D_dimension_distances);
norm_dist = cumsum(original_D_dimension_distances)/total_traversal;
p0angles = norm_dist*1/2*pi; 
% now converted locations in D dimensions to angles for mapping shell 0 points in 2D
shell0_locs = [cos(p0angles), sin(p0angles)]; % shell on unit  radius

end

%---------------------------------
function [p, original_D_dimension_distances] = ...
        local_approach(S, F, A, N)

permutation = seriate( S(F==1,F==1) ); % get best permutation of shell0 given similarity    
I = find(F==1);
II=I(permutation);
% pre computed permutation to use
p=1:size(A,1);
p(I)=II;
A0=A(F==1,:); % extract the non-dominated members
original_D_dimension_distances = zeros(N(1),1);
for i=1:N(1)-1
    original_D_dimension_distances(i+1) = Euc_dist(A0(permutation(i),:),A0(permutation(i+1),:));
end

end

%------------------
function [p, original_D_dimension_distances] = ...
        global_approach(S, F, A, N)

p = seriate( S ); % get best permutation of entire data set given similarity
I = find(F(p)==1); % indices of shell 1
     
original_D_dimension_distances = zeros(N(1),1);
for i=1:N(1)-1
    original_D_dimension_distances(i+1) = Euc_dist(A(p(I(i)),:),A(p(I(i+1)),:));
end
    
end

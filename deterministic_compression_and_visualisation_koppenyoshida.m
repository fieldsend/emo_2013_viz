function [ S, R, F, permutation, domination_matrix,est_dm,c] = ...
        deterministic_compression_and_visualisation_koppenyoshida( A, similarity_type, print_type, err_calc,P )

% function [ A, S, R, F, p, domination_matrix,c,dom_err,diff2] =
% deterministic_compression_and_visualisation_koppenyoshida( A,
% similarity_type, print_type, P )
%
% Method provides Pareto visualisation based on the description in
%
% M. Koppen and K. Yoshida. 
% Visualization of Pareto-sets in evolutionary multi-objective 
% optimization. 
% In Proceedings of the 7th International Conference on Hybrid Intelligent 
% Systems, pages 156?161, Washington, DC, USA, 2007
%
% As repoduced in
%
% Fieldsend JE, Everson RM. 
% Visualising high-dimensional Pareto relationships in two-dimensional 
% scatterplots, 
% Evolutionary Multi-criterion Optimisation, EMO 2013.
% LNCS pp 558-572
%    
%
% INPUTS
%
% A = set of objective vectors (columns objectives, rows elements)
% similarity_type = optional string holding name of function to generate 
%   similarity matrix from set A. There are four functions bundled with 
%   this group of files, though of course you may defined your own. The 
%   pre-defined ones are: 
%        (default) 'dominance_similarity_matrix' for dominance similarity 
%        (counts the objectives you have equivalent relationships with each 
%        other point to with respect to all other points in the set).
%        'strict_dominance_similarity_matrix' for strict dominance 
%         similarity (measures the proportion of points in the set that two
%         points both dominate)
%        'strict_dominated_similarity_matrix' for strict dominated 
%         similarity (measures the proportion of points in the set that two
%         points are both dominated by)
%        'combined_dominance_similarity_matrix' (measures the proportion of
%         points in the set that two points are either both dominated by or
%         both dominate)
% print_type = optional integer holding how links between points will be
%   plotted. (default) 1 for no links, 2 for domination links bewteen
%   adjecent shells/fronts and 3 for all domination links.
% err_calc = optional integer to denote calculation and printing of error
%   of fit measures (non-zero will print, default assumption is zero)
% P = a permuation of the non-dominated elements in A. In the original work
%   by Koppen and Yoshida this was accomplished via using NSGA-II (as also
%   used in the publication by Fieldsend and Everson). This is an optional 
%   input, if it is not provided, then spectral seriation is used to fix 
%   the order using the similarity matrix method passed in as the second 
%   argument (as described in Fieldsend and Everson). If P is passed in,
%   then the second argument will not be used (and can be set as e.g. an
%   empty string)
%
% OUTPUTS
%
% 
% S = matrix A converted into ranks
% R = Rank matrix equivalent of A
% F = front membership of each element of A
% p = permutation of A based on spectral seriation (depending on
%   the third input argment, this may be a complete permutation, or
%   just a permutation of the shell zero elements)
% domination_matrix = matrix tracking which elements dominate which
%   others (1 in element ji indicate member i dominates member j
%   0 indicates member i does not dominate member j
% est_dm = domination_matrix as inferred from 2D projection
% c = locations of mapped A points in 2D
%
% First plot produced is for a single pass through the data just using
% domination relationships to previous shells to set the position on the
% next shell, second plot is after 'passes' iterations though the data,
% where the dominating members of preceeding shells and the dominated 
% members of following shells are used to readjust the positions
%
% Colour of a point denotes the number of set members for which the
% location in 2D gives the wrong visual Pareto relationship
%
% Copyright (c) Jonathan E. Fieldsend 2013


if exist('similarity_type', 'var')==0
    similarity_type = 'dominance_similarity_matrix';
end

if exist('print_type', 'var')==0
    print_type =1;
end

if exist('err_calc','var')==0
    err_calc =0;
end
    
if exist('P','var')==0 % if a permuation of the shell0 elements is not provided
    [A, S, R, F, N, p] = seriate_pareto_front_and_convert_to_2Dshell( A, similarity_type);
    Ip = F(p)==1; % sorted order of nondom solutions
    indices_of_nondom = p(Ip); % non-dom indices shuffled by permuation
    [temp,P] = sort(indices_of_nondom); % convert indices down to range of front, but ordered appropriately
else
    [F] = recursive_pareto_shell_with_duplicates(A,1);
    for i=1:max(F);
        N(i) = sum(F==i);
    end
end
% P now contains the permuted order of just the rank 0 elements

I = find(F==1);
II=I(P);

% pre computed permutation to use
permutation =P; % permuation of the non dominated points only
p=1:size(A,1);
p(I)=II;

S=[];
R  = convert_raw_to_rank_matrix( A );
domination_matrix = get_dom_matrix( A ); % keep track of which members dominate which
%[F N] = get_fronts_via_dom_matrix( domination_matrix ); % get shell membership for elements of A
original_D_dimension_distances = zeros(N(1),1);


for i=1:N(1)-1
    original_D_dimension_distances(i+1) = Euc_dist(A(II(i),:),A(II(i+1),:));
end
total_traversal = sum(original_D_dimension_distances);
norm_dist = cumsum(original_D_dimension_distances)/total_traversal;
p0angles = norm_dist*1/2*pi; % now converted locations in D dimensions to angles for mapping shell 0 points in 2D

initial_set = [cos(p0angles), sin(p0angles)]; % shell on unit  radius

np = size(A,1);
c = zeros(np,2);
% now flip between axis
initial_set = initial_set *-1;
III = length(II):-1:1;
initial_set = initial_set(III,:);
c(II,:) = initial_set;



ndI = find(F(p)~=1);
for i=1:length(ndI);
    dI= domination_matrix(p(ndI(i)), II)==1;
    dps = c(II(dI),:); % get which members of shell 0 dominate this point
    c(p(ndI(i)),:) = max(dps,[],1);
end


est_dm = get_dom_matrix( c );
b = est_dm-domination_matrix;


dist_plotting(F,c,domination_matrix,p,[-1.1 0.1 -1.1 0.1],sum(abs(b),2),print_type);

if (err_calc)
    [Fc] = recursive_pareto_shell_with_duplicates(c,1);
    for i=1:max(Fc);
        Nc(i) = sum(Fc==i);
    end
    err=0;
    m = length(N);
    if (length(Nc)<length(N))
        m = length(Nc);
        err = err+N(m+1:end); % fewer shells so have definitely missed all the members of the shells above
    end
    for i=1:m
        a = setdiff(find(F==i),find(Fc==i));
        err = err+length(a);
    end
    
    fprintf('dominance error %d\n non-dom error %d\n shell error %d\n', ...
        sum(sum(abs(b==-1))), sum(sum(abs(b==1))),err);
end

end

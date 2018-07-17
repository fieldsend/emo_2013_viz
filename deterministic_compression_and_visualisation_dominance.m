function [ S, R, F, p, domination_matrix, est_dm, c] = ...
    deterministic_compression_and_visualisation_dominance( A, similarity_type,  print_type, err_calc)

% [ S, R, F, p, domination_matrix, est_dm, c] = 
%    deterministic_compression_and_visualisation_dominance( A, similarity_type,  print_type, err_calc)
%
% Method provides dominance-based visualisation described in:
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
%        'strict_dominance_similarity_matrix' for strict dominated 
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
% passes = optional non-negative integer holding the number of passes to be
%   undertaken in the refinement phase, defaults to 10 
%
% OUTPUTS
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

if (size(A,2) < 2)
   error('Number of objectives must be at least 2'); 
end

if exist('err_calc','var')==0
    err_calc = 0;
end

% first, seriate population and map the non-dominated subset to
% shell/front 0
disp('Seriating');
[A, S, R, F, N, p, shell0_locs, domination_matrix] = ...
   seriate_pareto_front_and_convert_to_2Dshell( A, similarity_type);
disp('Assigning shell 0');
%permutation =p;
np = size(A,1);
c = zeros(np,2);
% get indices and locations of shell 0
I = find(F(p)==1);
c(p(I),:)=shell0_locs;

% now go through each shell in turn, placing each individual in turn in the
% appropriate region

shellm = zeros(max(F),1);
% first set distance of each shell
shellm(1)=1;
for i=2:max(F)
    shellm(i) = sqrt(2)*shellm(i-1);
end

% now process each shell in turn
for i=2:max(F)
    disp(['Assigning shell ' int2str(i-1)]);
    % get the bounds of ecah of the regions delimited be by the
    % preceeding mapped shell
    r_bounds = zeros(2*length(I),2);
    range = 1:length(I);
    r_bounds(range,:) = [c(p(I),1) sqrt(shellm(i)^2-c(p(I),1).^2)];
    r_bounds(range+length(I),:) = [sqrt(shellm(i)^2-c(p(I),2).^2)  c(p(I),2)];
    [s, index] = sort(r_bounds(:,1));
    r_bounds = r_bounds(index,:);
    
    % get mid-point of each region as the putative point to map to
    % (i.e. to represent each region
    mid_points = (r_bounds(1:end-1,:)+r_bounds(2:end,:))/2;
    mid_points = mid_points./repmat(sqrt(sum(mid_points.^2,2)),1, 2)*shellm(i);
    
    % get the domination relationship between these mid-point and
    % all the solutions in the preceeding fronts, stored in a matrix
    dmAB = get_dom_between_matrices(c(p(I),:),mid_points);
    
    % now have matrix of domination relationships between dominant points
    % and the shell to be considered *and* between the dominant point and
    % all the candidate locations
    II = find(F(p)==i);
    locs = zeros(length(II),1);
    for j=1:length(II)
        %dmAB- repmat(domination_matrix(II(j),p(I)),2*length(I)-1,1)
        % element is 0 if both are domed, non-domed or doming
        % element is 1 r is domed but oringinal isn't
        % element is -1 if r is non-domed but original is
        diff = dmAB - repmat(domination_matrix(p(II(j)),p(I)),2*length(I)-1,1);
        % matrix with ones where solutions which should dominate are not
        diff_domerr = diff;
        diff_domerr(diff_domerr==1)=0;
        % matrix with ones where solutions which should not dominate are
        diff_nondomerr = diff;
        diff_nondomerr(diff_nondomerr==-1)=0;
        
        d=sum(abs(diff_domerr),2);
        m = min(d);
        if (m>0)
            error('theory broken -- a little sanity check needed');
        end
        mI = find(d==m);
       
        % OK, now have set of points with minimal domination error, now find
        % those with minimum non-dom error which are members of this subset
        diff_nondomerr = diff_nondomerr(mI,:);
        d_nd = sum(abs(diff_nondomerr),2);
        m_dnd = min(d_nd);
        mII = find(d_nd==m_dnd);
        
        locs(j) = mI(mII(1)); % take rightmost (cleaner visualisation
        c(p(II(j)),:) = mid_points(mI(mII(1)),:);
    end
    % check if any shared locations, if so, spread out
    for j=1:length(mid_points(:,1))
        if sum(locs==j)>1 % if more than one point in region
            Imm = find(locs==j);
            %c(p(II(Imm)),:) % display midpoint being mapped to
            scale = linspace(0,2,length(Imm)+2);
            mid_mid_points = zeros(length(Imm)+2,2);
            for k=2:length(scale)-1
                mid_mid_points(k,:) = (r_bounds(j,:)*(2-(scale(k)))+r_bounds(j+1,:)*scale(k))/2;
            end
            mid_mid_points = mid_mid_points./repmat(sqrt(sum(mid_mid_points.^2,2)),1,2)*shellm(i); %rescale on curve
            % locally seriate the
            loc_perm = seriate( S(p(II(Imm)),p(II(Imm))) );
            Imm =Imm(loc_perm);
            c(p(II(Imm)),:) = mid_mid_points(2:end-1,:);
        end
    end
    
    % update I to include the front just processed
    I = [I; II];
end

% now plot the data out, and also get some error (misinformation) measures
%figure

%plot(c(p(I),1), c(p(I),2), '.');

est_dm = get_dom_matrix( c );
b = est_dm-domination_matrix;
disp('Plotting...');
dist_plotting(F,c,domination_matrix,p,[0 max(shellm) 0 max(shellm)],sum(b,2), print_type);


if (err_calc)   
    fprintf('dominance error %d\n, non-dom error %d\n, shell error %d\n', sum(sum(abs(b==-1))), sum(sum(abs(b==1))),0);
end
    
end




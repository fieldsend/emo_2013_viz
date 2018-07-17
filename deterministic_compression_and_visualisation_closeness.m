function [ A, S, R, F, p, domination_matrix,est_dm,c] = ...
        deterministic_compression_and_visualisation_closeness( A, similarity_type, print_type, err_calc, passes)

% function [ A, S, R, F, p, domination_matrix,est_dm,c] =
% deterministic_compression_and_visualisation_closeness( A,
% similarity_type, print_type )
%
% Method provides distance-based visualisation described in
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
% Please contact me at J.E.Fieldsend@exeter.ac.uk for any queries
% or corrections
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


if exist('passes','var')==0
    passes = 10;
end

if passes < 0
    error('Cannot have a negative number of passes');
end


% first, seriate population and map the non-dominated subset to
% shell/front 0
disp('Seriating');
[A, S, R, F, N, p, shell0_locs, domination_matrix,ang0] = ...
    seriate_pareto_front_and_convert_to_2Dshell( A, similarity_type);
disp('Assigning shell 0');
%permutation =p;
np = size(A,1);
c = zeros(np,2);
I = find(F(p)==1);

% store in c the location of the set in 2D, and copy across the
% shell 0 solutions fixed via seriation 
c(p(I),:)=shell0_locs;
%shared_tracker=zeros(101,1);

% store in angles the angle used to fix each solution on each of
% the curves, and set the angles of the shell 0 ones
angles = zeros(np,1);
angles(p(I)) = ang0;

% store in doming_index initially the indices of shell 0
doming_index = I;
for i=2:max(F)
    % get the indices of the shell to be processed
    I = find(F(p)==i);
    % fix angles based on dominating solutions
    [locs,angs] = get_locations_and_angles(N(i), domination_matrix(p(I), p(doming_index)), angles(p(doming_index)),i);
    c(p(I),:)=locs;
    angles(p(I)) = angs;
    % the indices of the dominating solutions need to be updated
    % with hte shell just processed
    doming_index=[doming_index; I];
end
disp('Initial positions calculated');
disp('Plotting initial positions');

est_dm = get_dom_matrix( c );
b = est_dm-domination_matrix;
% plot the initial mapping
dist_plotting(F,c,domination_matrix,p,[0 max(F) 0 max(F)],sum(abs(b),2),print_type);

% now iteratively smooth the estimate using dominated as well as
% dominating solutions -- 10 passes has been more than enough for problems
% I've encountered thus far

disp('Refining positions using dominated information');
if (passes>-1)
    diff2=zeros(passes,1);
    for j=1:passes
        cc_new=c;
        for i=2:max(F)
            I = find(F(p)==i);
            [locations, angs] = get_loc_and_ang_all( N(i), domination_matrix(p(I),p), domination_matrix(p,p(I)),angles(p), i);
            cc_new(p(F(p)==i),:) = locations;
            angles(p(F(p)==i)) = angs;
        end
        diff2(j) = sum(sum((c-cc_new).^2));
        c=cc_new;
    end
    disp('Sum of squared location changes across points between each pass:');
    disp(diff2);
    disp('Refined positions calculated');
    disp('Plotting refined positions');
    
    est_dm = get_dom_matrix( c );
    b = est_dm-domination_matrix;
    dist_plotting(F,c,domination_matrix,p,[0 max(F) 0 max(F)],sum(abs(b),2),print_type);
end


if (err_calc)
    fprintf('dominance error %d\n, non-dom error %d\n, shell error %d\n', sum(sum(abs(b==-1))), sum(sum(abs(b==1))),0);
end

end

%--------
function [locs,angs] = get_loc_and_ang_all(N, dm, dm2,angles, shell)

    
% dm contains part of domination matrix relating to the members of the
% front below (columns) which dominate this shell (rows)

% dm2 contains part of domination matrix relating to the members of the
% front above (rows) which are dominated this by shell (columns)
    
locs = ones(N,2);
angs = zeros(N,1);
for i=1:N;
    doming_angles = [angles(dm(i,:)==1); angles(dm2(:,i)==1)]; % the angles of all the solutions which dominate the ith
                                                              % and dominated by the ith
    angs(i) = mean(doming_angles);
    locs(i,:) = [cos(angs(i))*shell, sin(angs(i))*shell];
end
end

%-----------
function [locs,angs] = get_locations_and_angles(N, dm, angles, shell)

    
% dm contains part of domination matrix relating to the members of the
% front below (columns) which dominate this shell (rows)
locs = ones(N,2);
angs = zeros(N,1);
for i=1:N;
    doming_angles = angles(dm(i,:)==1); % the angles of all the solutions which dominate the ith
    angs(i) = mean(doming_angles);
    locs(i,:) = [cos(angs(i))*shell, sin(angs(i))*shell];
end
end



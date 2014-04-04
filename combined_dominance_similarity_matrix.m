function [ S ] = combined_dominance_similarity_matrix( R )

% function [ S ] = dominance_similarity_matrix( R )
%
% 
% Copyright (c) Jonathan E. Fieldsend 2013
    
S = (strict_dominance_similarity_matrix( R ) + strict_dominated_similarity_matrix( R ))/2;


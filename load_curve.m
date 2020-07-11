function [ f, v ] = load_curve( file_name )
% [ f, v ] = load_curve( file_name )
%
% Loads a tabulated curve from a text file, one point per row.
%

fid = fopen( file_name, 'rt');
fv = transpose( fscanf( fid, '%e %e %e', [ 3 Inf ] ) );
fclose( fid );

% Frequency and values
f = fv(:,1);
v = fv(:,2) + i*fv(:,3);

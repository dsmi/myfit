function [ f, v ] = load_s2p( file_name )
% [ f, v ] = load_s2p( file_name )
%
% 
%

fid = fopen( file_name, 'rt');

% This reads the options line -- but so far we do nothing with it
optl = fgetl(fid);
    
fv = transpose( fscanf( fid, '%e %e %e %e %e %e %e %e %e', [ 9 Inf ] ) );

fclose( fid );

% Frequency and values
f = fv(:,1);

v = zeros( 2, 2, length(f) );
v(1,1,:) = fv(:,2) + i*fv(:,3);
v(1,2,:) = fv(:,4) + i*fv(:,5);
v(2,1,:) = fv(:,6) + i*fv(:,7);
v(2,2,:) = fv(:,8) + i*fv(:,9);

function [ rmserr, maxerr, err ] = calc_err( f, v, poles, resid, d )
% [ rmserr, maxerr, err ] = calc_err( f, v, poles, resid, d )
%
% Calculate the fitting error.
%

s = 2*pi*i*f;

ns = size( s, 1 ); % number of the samples
nv = size( v, 2 ); % number of the vector elements
np = size( poles, 1 );  % number of the poles

residr = repmat( permute( resid, [ 3 2 1 ] ), ns, 1 );
pr = repmat( permute( poles, [ 3 2 1 ] ), ns, nv );
vt = sum( residr ./ ( repmat( s, [ 1 nv np ] ) - pr ), 3 ) + d;
rmserr = sqrt( sum( sum( abs( ( vt - v ).^2 ), 2 ), 1 ) ) / sqrt( ns*nv );
err = abs( vt - v );
maxerr = max( err );

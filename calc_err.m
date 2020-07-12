function [ rmserr, maxerr ] = calc_err( f, v, poles, resid, d )
% [ f, v ] = calc_err( f, v, poles, resid, d )
%
% Calculate the fitting error.
%

s = 2*pi*i*f;

ns = size( s, 1 );
np = size( poles, 2 );

residr = repmat( resid.', ns, 1 );
pr = repmat( poles, ns, 1 );
vt = sum( residr ./ ( repmat( s, 1, np ) - pr ), 2 ) + d;
rmserr = sqrt( sum( abs( ( vt - v ).^2 ) ) ) / sqrt( ns );
maxerr = max( abs( vt - v ) );

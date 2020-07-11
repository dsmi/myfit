function [ poles, resid, rmserr ] = simplest_fit( f, v, npoles, niter )
% [ poles, resid, rmserr ] = simplest_fit( f, v, npoles, niter )
%
% Simplest possible implementation of the vector fitting algorithm.
%

% Add negative frequencies by mirroring! The first point is not
% mirrored if it is zero.
if 0.0 == f( 1 )
    frq = [ -flip( f(2:end), 1 ) ; f ];
    val = [ flip( real( v(2:end) ) + -i*imag( v(2:end) ), 1 ) ; v ];
else
    frq = [ -flip( f, 1 ) ; f ];
    val = [ flip( real( v ) + -i*imag( v ), 1 ) ; v ];
end

% Laplace variable, n-by-1
s = 2*pi*i*frq;

% Initial poles
startp = max( f(1), 10 );
endp   = f(end);
pval   = linspace( startp, endp, floor( npoles/2 ) );

% Place conjugate poles next to each other
poles = 0*[ pval pval ]; 
poles( 1:2:end ) = pval*1.0e-2 + i*pval;
poles( 2:2:end ) = pval*1.0e-2 - i*pval;

% Another real one if the number is odd
poles = [ poles endp*ones( 1, rem( npoles, 2 ) ) ];

ns = size( s, 1 );
np = size( poles, 2 );

% Iterations
for iter = 1:niter

    A1 = 1 ./ ( repmat( s, 1, np ) - repmat( poles, ns, 1 ) );
    A2 = [ ]; % [ ones( ns, 1 ) s ];
    A3 = -repmat( val, 1, np ) .* A1;

    A = [ A1 A2 A3 ];
    b = val;

    x = A\b;

    cfs = x( 1 : np );          % Residues of val*sigma approximation
    cs  = x( end-np+1 : end );  % Residues of sigma approximation

    H = diag( poles ) - ones( np, 1 ) * transpose( cs );

                                % Zeros of sigma
    zs = eig( H );

    poles = zs.';
end

% Find residues of approximated val
A1 = 1 ./ ( repmat( s, 1, np ) - repmat( poles, ns, 1 ) );
resid = A1\val;

% Calculate error
residr = repmat( resid.', ns, 1 );
pr = repmat( poles, ns, 1 );
vt = sum( residr ./ ( repmat( s, 1, np ) - pr ), 2 );
rmserr = sqrt( sum( abs( ( vt - val ).^2 ) ) ) / sqrt( ns );

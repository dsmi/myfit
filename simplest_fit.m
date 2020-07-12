function [ poles, resid, d, rmserr ] = simplest_fit( f, v, npoles, niter )
% [ poles, resid, d, rmserr ] = simplest_fit( f, v, npoles, niter )
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

% Matrix A, used for the pole identification, consists of three row blocks.
% Here is how a row of the A matrix looks like, below the block boundaries
% are shown:
%   Ak = [  1/(sk-a1) ... 1/(sk-aN) 1 sk  -f(sk)/(sk-a1) ... -f(sk)/(sk-aN) ]
%          |        A1             | A2  |                   A3            |

% Mid-block of the A matrix, which includes constant and linear terms
% and does not change.
A2 = [ ones( ns, 1 ) ];
%% A2 = [ ones( ns, 1 ) s ];

% Right hand vector is always the same.
b = val;

% Pole relocation iterations, plus another one to calculate A1 for
% the final residue calculation
for iter = 1:(niter+1)

    % First column block of the A matrix as shown above
    A1 = 1 ./ ( repmat( s, 1, np ) - repmat( poles, ns, 1 ) );

    % Is it the final n+1 iteration, where we only need to calculate A1?
    if iter > niter
        break
    end

    A3 = -repmat( val, 1, np ) .* A1;
    
    A = [ A1 A2 A3 ];
    
    x = A\b;

    % Residues of val*sigma approximation, and constant term (don't need here)
    % cfs = x( 1 : np );          
    % d = x( np + 1 );

    % Residues of sigma approximation
    cs  = x( end-np+1 : end );  

    H = diag( poles ) - ones( np, 1 ) * transpose( cs );

    % Zeros of sigma
    zs = eig( H );

    poles = zs.';
    
end

% Doing final residue calculation
A = [ A1 A2 ];

x = A\b;

% Residues and constant term
resid = x( 1 : np );          
d = x( np + 1 );

% Calculate error
residr = repmat( resid.', ns, 1 );
pr = repmat( poles, ns, 1 );
vt = sum( residr ./ ( repmat( s, 1, np ) - pr ), 2 ) + d;
rmserr = sqrt( sum( abs( ( vt - val ).^2 ) ) ) / sqrt( ns );

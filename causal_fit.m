function [ poles, resid, d, rmserr ] = causal_fit( f, v, weight, np, nit, asymp )
% [ poles, resid, d, rmserr ] = causal_fit( f, v, weight, np, nit, asymp )
%
% Simple implementation of the vector fitting, does not use negative
% frequencies, ensures perfect causal pole pairs, supports weight.
%

frq = f;
val = v;

% Laplace variable, n-by-1
s = 2*pi*i*frq;

% Initial poles
startp = max( f(1), 10 )*2*pi;
endp   = f(end)*2*pi;
pval   = linspace( startp, endp, floor( np/2 ) );

% Place complex poles
poles = [ ( -pval*1.0e-2 + i*pval )  ( -pval*1.0e-2 - i*pval ) ];

% Another real one if the number is odd
poles = [ poles -endp*ones( 1, rem( np, 2 ) ) ];

% Sort into complex conjugate pairs. Positive imaginary first.
poles = flipdim( cplxpair( poles ), 2 );

ns = size( s, 1 );
np = size( poles, 2 );

% Matrix A, used for the pole identification, consists of three row blocks.
% Here is how a row of the A matrix looks like, below the block boundaries
% are shown:
%   Ak = [  1/(sk-a1) ... 1/(sk-aN) 1 sk  -f(sk)/(sk-a1) ... -f(sk)/(sk-aN) ]
%          |        A1             | A2  |                   A3            |

% Mid-block of the A matrix, which includes constant and linear terms
% and does not change.
if asymp == 1
    A2 = [ ones( ns, 1 ) ];
else
    A2 = [ ];
end
%% A2 = [ ones( ns, 1 ) s ];

% Pole relocation iterations, plus another one to calculate A1 for
% the final residue calculation
for iter = 1:(nit+1)

    A1 = 1 ./ ( repmat( s, 1, np ) - repmat( poles, ns, 1 ) );

    % Indices of the complex poles
    cpi = find( imag( poles ) != 0 );

    % First and second elements of the conjugate pole pair (pairs are sorted)
    cpi1 = cpi( 1:2:end );
    cpi2 = cpi( 2:2:end );

    % Modify the corresponding elements of A
    ai = poles(cpi1);
    A1( :,cpi1 ) = 1 ./ ( s - ai ) + 1 ./ ( s - conj( ai ) );
    A1( :,cpi2 ) = j ./ ( s - ai ) - j ./ ( s - conj( ai ) );

    % Is it the final n+1 iteration, where we only need to calculate A1?
    if iter > nit
        break
    end

    A3 = -repmat( val, 1, np ) .* A1;

    % A and b both scaled by weight
    A = diag(weight)*[ A1 A2 A3 ];
    b = diag(weight)*val;

    Ar = [ real(A) ; imag(A) ]; % to preserve conjugacy
    br = [ real(b) ; imag(b) ];

    x = Ar \ br;

    %cfs = x( 1 : np );          % Residues of val*sigma approximation
    %cs  = x( end-np+1 : end );  % Residues of sigma approximation

    % Residues of sigma approximation, with conjugate pairs
    cs = x( end-np+1 : end );  % real
    cs(cpi1) = real( x( end-np+cpi1 ) ) + j*real( x( end-np+cpi2 ) );
    cs(cpi2) = real( x( end-np+cpi1 ) ) - j*real( x( end-np+cpi2 ) );

    % Calculate zeroes of sigma
    A = diag( poles );
    A( sub2ind( size(A), cpi1, cpi1 ) ) =  real( poles(cpi1) );
    A( sub2ind( size(A), cpi2, cpi2 ) ) =  real( poles(cpi1) );
    A( sub2ind( size(A), cpi1, cpi2 ) ) =  imag( poles(cpi1) );
    A( sub2ind( size(A), cpi2, cpi1 ) ) = -imag( poles(cpi1) );

    b = ones( np, 1 );
    b( cpi1 ) = 2;
    b( cpi2 ) = 0;

    cst = transpose( cs );
    cst( cpi1 ) = real( cs( cpi1 ) );
    cst( cpi2 ) = imag( cs( cpi1 ) );
    
    H = A - b * cst;

    % Zeros of sigma
    zs = eig( H );

    % Flip unstable poles
    unstables = real(poles)>0;  
    poles(unstables) = poles(unstables) - 2*real(poles(unstables));

    % Sort into complex conjugate pairs. Positive imaginary first.
    poles = flipdim( cplxpair( zs.' ), 2 );

end

% Doing final residue calculation
% A and b both scaled by weight
A = diag(weight)*[ A1 A2 ];
b = diag(weight)*val;

% to preserve conjugacy
Ar = [ real(A) ; imag(A) ]; 
br = [ real(b) ; imag(b) ];

x = Ar \ br;

% Residues, with conjugate pairs
resid = x( 1 : np ); % real
resid(cpi1) = real( x( cpi1 ) ) + j*real( x( cpi2 ) );
resid(cpi2) = real( x( cpi1 ) ) - j*real( x( cpi2 ) );

% Constant term
if asymp == 1
    d = x( np + 1 );
else
    d = 0;
end


% Calculate error
residr = repmat( resid.', ns, 1 );
pr = repmat( poles, ns, 1 );
vt = sum( residr ./ ( repmat( s, 1, np ) - pr ), 2 ) + d;
rmserr = sqrt( sum( abs( ( vt - val ).^2 ) ) ) / sqrt( ns );

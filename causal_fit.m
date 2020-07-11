function [ poles, resid, rmserr ] = causal_fit( f, v, npoles, niter )
% [ poles, resid, rmserr ] = causal_fit( f, v, npoles, niter )
%
% Simple implementation of the vector fitting, does not use negative
% frequencies, ensures perfect causal pole pairs.
%

frq = f;
val = v;   
    
% Laplace variable, n-by-1
s = 2*pi*i*frq;

% Initial poles
startp = max( f(1), 10 );
endp   = f(end);
pval   = linspace( startp, endp, floor( npoles/2 ) );

% Place complex poles
poles = [ ( pval*1.0e-2 + i*pval )  ( pval*1.0e-2 - i*pval ) ];
poles = cplxpair( poles ); % sort into complex conjugate pairs

% Another real one if the number is odd
poles = [ poles endp*ones( 1, rem( npoles, 2 ) ) ];

ns = size( s, 1 );
np = size( poles, 2 );

% Iterations
for iter = 1:niter

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

    %% % Modify the corresponding elements of A
    %% ai = poles(1);
    %% A1( :,1 ) = 1 ./ ( s - ai ) + 1 ./ ( s - conj( ai ) );
    %% A1( :,2 ) = j ./ ( s - ai ) - j ./ ( s - conj( ai ) );
    
    A2 = [ ]; % [ ones( ns, 1 ) s ];
    A3 = -repmat( val, 1, np ) .* A1;

    A = [ A1 A2 A3 ];
    b = val;

    x = [ real(A) ; imag(A) ] \ [ real(b) ; imag(b) ]; % to preserve conjugacy

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

    poles = cplxpair( zs.' ); % sort into complex conjugate pairs

end

% Find residues of approximated val
A1 = 1 ./ ( repmat( s, 1, np ) - repmat( poles, ns, 1 ) );
resid = A1\val;

% Calculate error
residr = repmat( resid.', ns, 1 );
pr = repmat( poles, ns, 1 );
vt = sum( residr ./ ( repmat( s, 1, np ) - pr ), 2 );
rmserr = sqrt( sum( abs( ( vt - val ).^2 ) ) ) / sqrt( ns );

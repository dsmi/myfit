function [ poles, resid, d, rmserr ] = vector_fit( f, v, weight, np, nit, asymp )
% [ poles, resid, d, rmserr ] = vector_fit( f, v, weight, np, nit, asymp )
%
% Almost complete implementation of the vector fitting, unity non-triviality
% constant.
%

% Laplace variable, n-by-1
s = 2*pi*i*f;

% number of the vector elements
nv = size( v, 2 );

% number of the frequency samples
ns = size( s, 1 );

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

    
    A = [ ];
    b = [ ];
    
    for iv = 1:nv
        vi = v(:,iv);            % samples of this element
        wi = diag(weight(:,iv)); % weights of this element

        A12i = [ real( wi*[ A1 A2 ] )     ; imag( wi*[ A1 A2 ] ) ];
        A3i  = [ real( -diag( vi ) * A1 ) ; imag( -diag( vi ) * A1 ) ];
        
        %% 'slow' VF
        %% A12 = blkdiag( A12, A12i );
        %% A3  = [ A3 ; A3i ];
        %% bb  = [ bb ; real(wi*vi) ; imag(wi*vi) ];

        %% fast VF -- diagonalize A matrix, keep the last column block only
        [ Q, R ] = qr( [ A12i A3i ], 0 );
        R22 = R(end-np+1:end,end-np+1:end);
        
        A = [ A ; R22 ];
        b = [ b ; Q(:,end-np+1:end)'*[ real(wi*vi) ; imag(wi*vi) ] ];
        
    end

    %% 'slow' VF
    %% A = [ A12 A3 ];
    %% x = A\b;

    x = A \ b;

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
resid = zeros( np, nv );
d = zeros( 1, nv );

for iv = 1:nv

    vi = v(:,iv);            % samples of this element
    wi = diag(weight(:,iv)); % weights of this element
    
    % A and b both scaled by weight
    A = wi*[ A1 A2 ];
    b = wi*vi;

    % to preserve conjugacy
    A = [ real(A) ; imag(A) ]; 
    b = [ real(b) ; imag(b) ];

    x = A \ b;

    % Residues, with conjugate pairs
    resid(:,iv) = x( 1 : np ); % real
    resid(cpi1,iv) = real( x( cpi1 ) ) + j*real( x( cpi2 ) );
    resid(cpi2,iv) = real( x( cpi1 ) ) - j*real( x( cpi2 ) );

    % Constant term
    if asymp == 1
        d(iv) = x( np + 1 );
    end

end

% To have poles and residues arranged the same way, as the column vectors
poles = transpose( poles );

% Calculate error
residr = repmat( permute( resid, [ 3 2 1 ] ), ns, 1 );
pr = repmat( permute( poles, [ 3 2 1 ] ), ns, nv );
vt = sum( residr ./ ( repmat( s, [ 1 nv np ] ) - pr ), 3 ) + d;
rmserr = sqrt( sum( sum( abs( ( vt - v ).^2 ), 2 ), 1 ) ) / sqrt( ns*nv );

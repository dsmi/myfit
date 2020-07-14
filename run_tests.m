
% Load the tabulated curve to fit
[ f, v2 ] = load_s2p( "ladder3.s2p" );

% Order of approximation, number of poles
npoles = 5;

% Number of the pole relocating iterations
niter = 3;

% First one to fit
v = squeeze( v2(1,1,:) );

[ poles, resid, d ] = simplest_fit( f, v, npoles, niter );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid, d );
assert( rmserr < 1e-6 )
assert( maxerr < 1e-6 )

[ poles, resid, d ] = causal_fit( f, v, f*0+1, npoles, niter, 1 );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid, d );
assert( rmserr < 1e-6 )
assert( maxerr < 1e-6 )

% second one to fit
v = squeeze( v2(2,1,:) );

[ poles, resid, d ] = simplest_fit( f, v, npoles, niter );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid, d );
assert( rmserr < 1e-6 )
assert( maxerr < 1e-6 )

[ poles, resid, d ] = causal_fit( f, v, f*0+1, npoles, niter, 1 );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid, d );
assert( rmserr < 1e-7 )
assert( maxerr < 1e-7 )

%
% Weighted fit test
%

% Load the tabulated curve to fit
[ fhz, v2 ] = load_s2p( "test_pdn.s2p" );

% curve to fit
v = squeeze( v2(2,1,:) );

% Order of approximation, number of poles
npoles = 200;

% Number of the pole relocating iterations
niter = 2;

% Uniform weights
weight = 1+v*0;
[ poles, resid, d ] = causal_fit( fhz, v, weight, npoles, niter, 0 );
[ rmserr, maxerr, err ] = calc_err( fhz, v, poles, resid, d );
assert( rmserr < 1e-3 )
assert( err(1) < 1e-4 )
assert( err(1) > 1e-6 )

% Inverted sample magnitude weights
weight = 1./max(sqrt(abs(v)),1e-5);
[ poles, resid, d ] = causal_fit( fhz, v, weight, npoles, niter, 0 );
[ rmserr, maxerr, err ] = calc_err( fhz, v, poles, resid, d );
assert( rmserr < 1e-3 )
assert( err(1) < 1e-12 )





%% % number of ports
%% Np = size( v2, 1 );

%% % vector to be fitted, Nc-by-Ns
%% v = reshape( v2, [], Np*Np );
%% size(v)

%% [ poles, resid, d ] = vector_fit( f, v, npoles, niter );

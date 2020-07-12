
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

[ poles, resid, d ] = causal_fit( f, v, npoles, niter );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid, d );
assert( rmserr < 1e-11 )
assert( maxerr < 1e-11 )

% second one to fit
v = squeeze( v2(2,1,:) );

[ poles, resid, d ] = simplest_fit( f, v, npoles, niter );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid, d );
assert( rmserr < 1e-6 )
assert( maxerr < 1e-6 )

[ poles, resid, d ] = causal_fit( f, v, npoles, niter );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid, d );
assert( rmserr < 1e-7 )
assert( maxerr < 1e-7 )

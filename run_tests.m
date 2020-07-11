
%% Load the tabulated curve to fit
[ f, v ] = load_curve( 'ladder3_s21.txt' );
    
% Order of approximation, number of poles
npoles = 6;

% Number of the pole relocating iterations
niter = 3;

[ poles, resid ] = simplest_fit( f, v, npoles, niter );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid )

[ poles, resid ] = causal_fit( f, v, npoles, niter );
[ rmserr, maxerr ] = calc_err( f, v, poles, resid )

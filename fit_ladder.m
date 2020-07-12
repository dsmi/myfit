
% Load the tabulated curve to fit
[ f, v2 ] = load_s2p( "ladder3.s2p" );
v = squeeze( v2(2,1,:) );
%% [ f, v ] = load_curve( 'ladder3_s21.txt' );
    
% Order of approximation, number of poles
npoles = 5;

% Number of the pole relocating iterations
niter = 5;

[ poles, resid, d, rmserr ] = causal_fit( f, v, npoles, niter );

poles.'
resid
rmserr

% poles
do_plots = 1;

if do_plots

    % Mirror frequencies for the plots
    frq = [ -flip( f(2:end), 1 ) ; f ];
    val = [ flip( real( v(2:end) ) + -i*imag( v(2:end) ), 1 ) ; v ];

    % Laplace frequencies
    s = 2*pi*i*frq;

    ns = size( s, 1 );
    np = size( poles, 2 );

    % Plot source and approximated function
    figure(1)
    residr = repmat( resid.', ns, 1 );
    pr = repmat( poles, ns, 1 );
    vt = sum( residr ./ ( repmat( s, 1, np ) - pr ), 2 ) + d;
    maxerr = max( abs( vt - val ) );
    plot( s/i, abs(val), '*r', s/i, abs(vt), '-b' );
    legend( 'source', 'fit' )

    %% % Plot poles
    %% figure(2)
    %% plot( real(poles), imag(poles), '*' )
    %% xlim( [ -1e10 1e10 ] )
    %% ylim( [ -6e10 6e10 ] )
    %% xlabel('real(s)')
    %% ylabel('imag(s)')

    %% % 3d surface plot
    %% figure(3)

    %% maxr = max( abs( real( poles ) ) );
    %% maxi = max( abs( imag( poles ) ) );

    %% nr = 20*2;
    %% ni = 120*2;

    %% mr = linspace( -2e10, 0, nr );
    %% mi = linspace( -6e10, 6e10, ni );

    %% [ meshr, meshi ] = meshgrid( mr, mi );

    %% s = reshape( meshr + i*meshi, [ ], 1 );

    %% ns = size( s, 1 );
    %% residr = repmat( resid.', ns, 1 );
    %% pr = repmat( poles, ns, 1 );
    %% ft = sum( residr ./ ( repmat( s, 1, np ) - pr ), 2 );

    %% ftmesh = reshape( ft, ni, nr );

    %% surf( meshr, meshi, min( 2.5, abs(ftmesh) ) )
    %% xlim( [ -6e10 6e10 ] )
    %% ylim( [ -6e10 6e10 ] )
    %% zlim( [ 0 2.6 ] )
    %% caxis( [ 0 2.6 ] )
    %% xlabel('real(s)')
    %% ylabel('imag(s)')
    %% colormap(jet)
    
end

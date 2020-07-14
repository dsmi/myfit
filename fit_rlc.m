% Fitting of RLC circuit

R=10;
L=5e-9;
C=5e-12;

% frequency samples for the Y plot
fhz = linspace(0,1e10,30000);

% Laplace variable
s=i*2*pi*fhz;

% Laplace-domain admittance
Y=s./(s.*s.*L+s*R+1/C);

%% plot(fhz, real(Y), '-r', fhz, imag(Y), '-b');

% Laplace-domain admittance, partial fractions
D=R*R-4*L/C;
s1=(-R-sqrt(D))/2/L;
s2=(-R+sqrt(D))/2/L;
c1=s1./(L*(s1-s2));
a1=s1;
c2=s2./(L*(s2-s1));
a2=s2;
Y2=c1./(s-a1)+c2./(s-a2);

% Order of approximation, number of poles
npoles = 2;

% Number of the pole relocating iterations
niter = 1;

[ poles, resid, d, rmserr ] = causal_fit( fhz.', Y.', (Y.')*0+1, npoles, niter, 0 );
rmserr

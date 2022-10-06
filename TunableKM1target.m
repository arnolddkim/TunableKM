%% TunableKM1target.m
%
% Compute the KM and modified KM images for a single point target. 
%
% This code was used to produce the results shown in Fig. 2 of the 
% manuscript, "Tunable high-resolution synthetic aperture radar imaging" 
% by A. D. Kim and C. Tsogka.
%
% Written by A. D. Kim on 3/2/2022

clear;

%% FIGURE PARAMETERS

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 
  
%% PHYSICAL PARAMETERS

% system bandwidth

f0 = 9.6e9; % central frequency
B  = 622e6; % bandwidth
M  = 16;
f  = linspace( f0 - B / 2, f0 + B / 2, 2 * M - 1 );
c0 = 299705e3;
k0 = 2 * pi * f0 / c0;
k  = 2 * pi * f / c0;

% add negative wavenumbers

%k = [ -fliplr(k) k ];
%

lambda0 = c0 / f0;

% synthetic aperture

H = 7300;
R = 7100;
L = sqrt( H^2 + R^2 );
N = 124;
a = 130;

xa = linspace( -a / 2, a / 2, N );
ya = R;
za = H;

% location of point target

x1 = 1.0;
y1 = 1.0;
z1 = 0.0;

% complex reflectivity of the point target

rho = 3.4 * 1j;

%% MEASUREMENT MATRIX

% meshgrid of wavenumbers and distances

[ K, R ] = ndgrid( k, sqrt( ( xa - x1 ).^2 + ( ya - y1 )^2 + ( za - z1 )^2 ) );

% Green's functions

G = exp( 1j * K .* R ) ./ ( 4 * pi * R );

% measurement matrix

BB = rho * G .* G;

%% KIRCHHOFF MIGRATION

% mesh for the imaging window

Ns  = 51;

xscale = 100 * lambda0;
yscale = 50 * lambda0;

xIW = linspace( -xscale / 2, xscale / 2, Ns );
yIW = linspace( -yscale / 2, yscale / 2, Ns );

[ X, Y ] = meshgrid( xIW, yIW );

% allocate memory for the image

KM = 0 * X;

% compute the KM image

for n = 1 : N
    
    % isolate column of data matrix

    dhat  = BB(:,n);

    for i = 1 : Ns * Ns
        
        % illumination vector
        
        rn = sqrt( ( xa(n) - x1 - X(i) ).^2 + ( ya - y1 - Y(i) ).^2 + za.^2 );
        an = exp( 1j * 2 * k' * rn );
        an = an / norm( an );

        % update imaging function

        KM(i) = KM(i) + dhat' * an;
        
    end
    
end

KM = abs( KM ).^2;
KM = KM / max( KM(:) );

%% compute tunable-KM image

epsilon = 1e-2;

IKM = epsilon ./ ( 1 - ( 1 - epsilon ) * KM );

%% IMAGE

figure(1)
pcolor( xIW / lambda0, yIW / lambda0, 10 * log10( IKM ) );
colorbar;
shading flat;
xlabel( '$( x - x_{1} )/\lambda_{0}$', 'Interpreter', 'LaTeX' );
ylabel( '$( y - y_{1} )/\lambda_{0}$', 'Interpreter', 'LaTeX' );

figure(2)
pcolor( xIW / lambda0, yIW / lambda0, 10 * log10( KM ) );
colorbar;
shading flat;
xlabel( '$( x - x_{1} )/\lambda_{0}$', 'Interpreter', 'LaTeX' );
ylabel( '$( y - y_{1} )/\lambda_{0}$', 'Interpreter', 'LaTeX' );

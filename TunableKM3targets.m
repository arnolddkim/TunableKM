%% TunableKM3targets.m
%
% Compute the KM and modified KM images for three point targets. 
%
% This code was used to produce the results shown in Figs. 6 and 7 of the 
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

% % add negative wavenumbers
% 
% k = [ -fliplr(k) k ];

% synthetic aperture

H = 7300;
R = 7100;
L = sqrt( H^2 + R^2 );
N = 124;
a = 129.15;

xa = linspace( -a / 2, a / 2, N );
ya = R;
za = H;

% location of point target

x1 = -1.4;
y1 = -0.5;
z1 =  0.0;

x2 = -0.6;
y2 = -1.2;
z2 = 0.0;

x3 = 1.2;
y3 = 1.1;
z3 = 0.0;

% complex reflectivity of the point target

rho1 = 3.4 * 1j;
rho2 = 3.4 * 1j;
rho3 = 3.4 * 1j;

%% MEASUREMENT MATRIX

% meshgrid of wavenumbers and distances

[ K, R1 ] = ndgrid( k, sqrt( ( xa - x1 ).^2 + ( ya - y1 )^2 + ( za - z1 )^2 ) );
[ ~, R2 ] = ndgrid( k, sqrt( ( xa - x2 ).^2 + ( ya - y2 )^2 + ( za - z2 )^2 ) );
[ ~, R3 ] = ndgrid( k, sqrt( ( xa - x3 ).^2 + ( ya - y3 )^2 + ( za - z3 )^2 ) );

% Green's functions

G1 = exp( 1j * K .* R1 ) ./ ( 4 * pi * R1 );
G2 = exp( 1j * K .* R2 ) ./ ( 4 * pi * R2 );
G3 = exp( 1j * K .* R3 ) ./ ( 4 * pi * R3 );

% measurement matrix

BB = rho1 * G1 .* G1 + rho2 * G2 .* G2 + rho3 * G3 .* G3;

%% add measurement noise

%rng(1);

sn     = 1e-0;
B0     = BB;
enr    = sqrt(sum(abs(BB(:)).^2)/length(BB(:)));
Bnoise = enr*sn*(randn(size(BB))+1i*randn(size(BB)))/sqrt(2);
BB     = B0 + Bnoise;
 
SNR_est  = -10 * log10(sn);
SNR_real = 10 * log10(norm(B0)/norm(Bnoise));

disp( ' ' );
disp( '---' );
disp( 'SNR' );
disp( '---' );
disp( ' ' );
disp( ['  Realization SNR = ', num2str(SNR_real) ] );
disp( ['  Average SNR     = ', num2str(SNR_est)  ] );
disp( ' ' );

%% KIRCHHOFF MIGRATION

% coarse mesh for the imaging window

Ns = 65;

xscale = 4.8;
yscale = 4.8;

xIW = linspace( -xscale/2, xscale/2, Ns );
yIW = linspace( -yscale/2, yscale/2, Ns );

[ X, Y ] = meshgrid( xIW, yIW );

% fine mesh for the imaging window

xscale0 = 5 * c0 / f0;
yscale0 = 5 * c0 / f0;

xIW0 = linspace( -xscale0/2, xscale0/2, Ns );
yIW0 = linspace( -yscale0/2, yscale0/2, Ns );

[ X0, Y0 ] = meshgrid( xIW0, yIW0 );

% allocate memory for the images

KM = 0 * X;

KM1 = 0 * X0;
KM2 = 0 * X0;
KM3 = 0 * X0;

% compute the KM images

for n = 1 : N
    
    % isolate column of data matrix

    dhat  = BB(:,n) / norm( BB(:,n) );

    for i = 1 : Ns * Ns
        
        % illumination vectors

        r = sqrt( ( xa(n) - X(i) ).^2 + ( ya - Y(i) ).^2 + za.^2 );
        a = exp( 1j * 2 * k' * r );
        a = a / norm( a );
        
        r1 = sqrt( ( xa(n) - x1 - X0(i) ).^2 + ( ya - y1 - Y0(i) ).^2 + za.^2 );
        a1 = exp( 1j * 2 * k' * r1 );
        a1 = a1 / norm( a1 );

        r2 = sqrt( ( xa(n) - x2 - X0(i) ).^2 + ( ya - y2 - Y0(i) ).^2 + za.^2 );
        a2 = exp( 1j * 2 * k' * r2 );
        a2 = a2 / norm( a2 );

        r3 = sqrt( ( xa(n) - x3 - X0(i) ).^2 + ( ya - y3 - Y0(i) ).^2 + za.^2 );
        a3 = exp( 1j * 2 * k' * r3 );
        a3 = a3 / norm( a3 );

        % update imaging functions

        KM(i) = KM(i) + dhat' * a;
        
        KM1(i) = KM1(i) + dhat' * a1;
        KM2(i) = KM2(i) + dhat' * a2;
        KM3(i) = KM3(i) + dhat' * a3;
        
    end
    
end

KM = abs( KM ).^2;
KM = KM / max( KM(:) );

KM1 = abs( KM1 ).^2;
KM2 = abs( KM2 ).^2;
KM3 = abs( KM3 ).^2;

KM1 = KM1 / max( KM1(:) );
KM2 = KM2 / max( KM2(:) );
KM3 = KM3 / max( KM3(:) );

%% tunable KM

epsilon = 1e-4;

IKM = epsilon ./ ( 1 - ( 1 - epsilon ) * KM );

IKM1 = epsilon ./ ( 1 - ( 1 - epsilon ) * KM1 );
IKM2 = epsilon ./ ( 1 - ( 1 - epsilon ) * KM2 );
IKM3 = epsilon ./ ( 1 - ( 1 - epsilon ) * KM3 );

%% IMAGES

lambda0 = c0 / f0;

figure(10)

pcolor( xIW / lambda0, yIW / lambda0, 10*log10(KM) );
colorbar;
shading flat;
hold on;
plot( x1 / lambda0, y1 / lambda0, 'r+' );
plot( x2 / lambda0, y2 / lambda0, 'r+' );
plot( x3 / lambda0, y3 / lambda0, 'r+' );
hold off;
xlabel( '$x/\lambda_{0}$', 'Interpreter', 'LaTeX' );
ylabel( '$y/\lambda_{0}$', 'Interpreter', 'LaTeX' );

figure(11)

pcolor( xIW / lambda0, yIW / lambda0, 10*log10(IKM) );
colorbar;
shading flat;
hold on;
plot( x1 / lambda0, y1 / lambda0, 'r+' );
plot( x2 / lambda0, y2 / lambda0, 'r+' );
plot( x3 / lambda0, y3 / lambda0, 'r+' );
hold off;
xlabel( '$x/\lambda_{0}$', 'Interpreter', 'LaTeX' );
ylabel( '$y/\lambda_{0}$', 'Interpreter', 'LaTeX' );

figure(1)

pcolor( xIW0 / lambda0, yIW0 / lambda0, 10*log10(IKM1) );
colorbar;
shading flat;
xlabel( '$(x - x_{1})/\lambda_{0}$', 'Interpreter', 'LaTeX' );
ylabel( '$(y - y_{1})/\lambda_{0}$', 'Interpreter', 'LaTeX' );

figure(2)

pcolor( xIW0 / lambda0, yIW0 / lambda0, 10*log10(IKM2) );
colorbar;
shading flat;
xlabel( '$(x - x_{2})/\lambda_{0}$', 'Interpreter', 'LaTeX' );
ylabel( '$(y - y_{2})/\lambda_{0}$', 'Interpreter', 'LaTeX' );

figure(3)

pcolor( xIW0 / lambda0, yIW0 / lambda0, 10*log10(IKM3) );
colorbar;
shading flat;
xlabel( '$(x - x_{3})/\lambda_{0}$', 'Interpreter', 'LaTeX' );
ylabel( '$(y - y_{3})/\lambda_{0}$', 'Interpreter', 'LaTeX' );
function SNorm = CalcVonKarmanIEC( n, U, z, normalise )
%CALCVONKARMANIEC Calculate normalised Von Karman model
%   Based on equations from:
%       T. Burton, N. Jenkins, D. Sharpe, and E. Bossanyi, "Wind Energy
%       Handbook (2nd Edition)," John Wiley & Sons, 2011, pp. 9-21.
%   Inputs:
%       - n: frequencies to compute the spectrum at (Hz)
%       - U: mean wind speed (m/s)
%       - z: height above the ground (m)
%       - nomalise: whether to normalise the PSD by the frequency n
%   Outputs:
%       - SNorm: power spectral density, nomalised by mean wind speed U and
%                (optionally) by frequency n
%   Written by:  Z.J. Chen, 2017
%   Modified by: J.X.J. Bannwarth, 2020/11/17

% Default parameters for countryside with trees and hedges
z0          = 0.1;            % (m) surface roughness length 
zi          = 1000*z0^0.18;   % (m) min height for isotropic turbulence
angVelEarth = 7.2921159e-5;   % (rad/s) angular rotation speed of earth
latAkl      = 36.8485*pi/180; % (rad) latitude of Auckland
kappa       = 0.4;            % (-) von Karman constant

% Turbulence intensity
fCor = 2 * angVelEarth * sin( abs(latAkl) ); % (rad/s) Coriolis parameter
uFric = (kappa*U - 34.5*fCor*z) / log(z/z0); % (m/s) friction velocity
eta = 1 - 6*fCor*z/uFric;
p = eta^16;

h = uFric / (6*fCor);
% Standard deviation of longitudinal component
stdU = ( 7.5*eta*(0.538 + 0.09*log(z/z0))^p*uFric ) / ...
    ( 1 + 0.156*log(uFric/fCor*z0) );

% Standard deviation of the other two components
stdV = stdU * ( 1 - 0.22*cos(pi*z/(2*h))^4 );
stdW = stdU * ( 1 - 0.45*cos(pi*z/(2*h))^4 );

% Return normalised von Karman spectrum - Change normalisation?
SNorm = zeros( length(n), 3 );
L = zeros( 3, 1 );
for i = 1:3
    if (i == 1)
        L(i) = 280*(z/zi)^0.35;
        SNorm(:,i) = (4*n*L(i)/U)./(1+70.8*(n*L(i)/U).^2).^(5/6);
        stdWind = stdU;
    elseif (i == 2)
        L(i) = 140*(z/zi)^0.48;
        SNorm(:,i) = (4*(n*L(i)/U).*(1+755.2*(n*L(i)/U).^2))./(1+283.2*(n*L(i)/U).^2).^(11/6);
        stdWind = stdV;
    elseif (i == 3)
        L(i) = 0.35*z;
        SNorm(:,i) = (4*(n*L(i)/U).*(1+755.2*(n*L(i)/U).^2))./(1+283.2*(n*L(i)/U).^2).^(11/6);
        stdWind = stdW;
    end
    SNorm(:,i) = SNorm(:,i)*stdWind^2 ./ (U^2);
    if ~normalise
        SNorm(:,i) = SNorm(:,i) ./ n';
    end
end
return
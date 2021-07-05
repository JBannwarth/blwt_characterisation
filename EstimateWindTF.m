%ESTIMATEWINDTF Estimate wind transfer function using measured PSD
%   By 'wind transfer function', we mean a transfer function that shapes
%   white noise to obtain a similar PSD to that of the wind.
%   Written by: J. X. J. Bannwarth, 2020/11/18
%% Clean-up and load data
clearvars; clc; close all
folderIn = 'data';

% Get files
files = dir( fullfile( folderIn, 'coarsegrid_center_*_2.thA' ) );
fileNames = {files.name}';
namesNoExt = erase( fileNames, '.thA' );

% Load and get frequency response
UMean = zeros( length(namesNoExt), 3 );
UStd = UMean;
UTi = UMean;
f = logspace(-2,2,501)'; % Hz
for ii = 1:length(namesNoExt)
    % Load data from file
    [ u, v, w, Ps, Settings ] = ReadTHFile( fullfile( folderIn, fileNames{ii}) );
    
    % Fix misalignments
    [ u, v, w ] = ConditionWindSpeed( u, v, w, Settings.DataRate );
    
    % Write to struct
    windData(ii).U = u;
    windData(ii).V = v;
    windData(ii).W = w;
    windData(ii).Ps = Ps;
    
    % Statistics
    U = [u v w];
    UMean(ii,:) = mean(U);
    UStd(ii,:) = std(U);
    UTi(ii,:) = 100.* UStd(ii,:) ./ UMean(ii,1);
    
    % Get spectrum
    [windData(ii).pU, windData(ii).f] = pwelch( ...
        windData(ii).U-UMean(ii,1), hann(1024, 'periodic'), ...
        [], f, Settings.DataRate );
    windData(ii).pU = windData(ii).pU';
    windData(ii).f = windData(ii).f';
    
    [windData(ii).pV, ~] = pwelch( ...
        windData(ii).V-UMean(ii,2), hann(1024, 'periodic'), ...
        [], f, Settings.DataRate );
    windData(ii).pV = windData(ii).pV';
    
    [windData(ii).pW, ~] = pwelch( ...
        windData(ii).W-UMean(ii,3), hann(1024, 'periodic'), ...
        [], f, Settings.DataRate );
    windData(ii).pW = windData(ii).pW';
    
    % Convert to decibels and rad/s
    windData(ii).MU = 20*log10( windData(ii).pU );
    windData(ii).MV = 20*log10( windData(ii).pV );
    windData(ii).MW = 20*log10( windData(ii).pW );
    windData(ii).w = 2*pi*windData(ii).f;
end

% Export
figSize = [15, 10];
fontSize = 11;

%% Plot
% Select magnitude to plot
w = windData(3).w;
M = windData(3).MU;
p = windData(3).pU;
f = windData(3).f;

% Magnitude response
figure( 'name', 'Magnitude response' )
semilogx( f, M )
hold on; grid on; box on
xlabel( 'Frequency (Hz)' )
ylabel( 'Magnitude (dB)' )
SetFigProp( figSize, fontSize )

sys = frd( p, w, 1/1250 );

% Convert to TF
windTF = tfest( sys, 2 );

% Manually fit - way number 1
dcGain = mean( M(w < 1) );
MHighF = M( w > 50 );
fHighF = f( w > 50 );
c = polyfit( log10(fHighF), MHighF, 1 );
fCornerLog = ( dcGain - c(2) ) / c(1);
fCorner = 10^( fCornerLog );
wCorner = 2*pi*fCorner;

plot( [f(1) fCorner f(end)], dcGain*[1 1 1] )
plot( [fCorner fHighF(end)], c(2) + c(1).*log10([fCorner fHighF(end)]) )
plot( fCorner, dcGain, 'o' )

% Manual fit - way number 2
ind = find( M <= (dcGain-3), 1 );
fCorner2 = f( ind );
wCorner2 = 2*pi*fCorner2;

plot( fCorner2, dcGain, '*' )

legend( { 'Measured';
            sprintf( 'DC gain = %.1f dB', dcGain );
            sprintf( 'High $f$ asymptote, -%.1f dB/dec', -c(1) );
            sprintf( '$\\omega_c$ = %.1f rad/s', wCorner );
            sprintf( '$\\omega_{c,2}$ = %.1f rad/s', wCorner2 ) }, ...
        'location', 'southwest' )

% Create filters
[num, den] = butter( 2, wCorner2, 's');
windTFManual1 = 10^(dcGain/20)*tf( wCorner2, [1 wCorner2]);
windTFManual2 = 10^(dcGain/20)*tf( num, den );

% Compare
figure( 'Name', 'Bode plots' )
bode( windTF )
hold on; grid on; box on;
bode( windTFManual1 )
bode( windTFManual2 )
bode( sys )
legend( {'Estimated', 'Estimated (manual 1st order)', ...
    'Estimated (manual 2nd order)', 'Experiment'} )
xlim( [w(1) w(end)] )
SetFigProp( figSize, fontSize )

% Based on all of this, a 1st order transfer function seems to be the most
% appropriate
% Alternative - Use the theoretical PSD and just try to match the high freq
% asymptote? Since outdoors we would expect the corner freq to be much
% lower. But we could just say we tuned to the wind tunnel.

%% Compare different axes
% Get corner frequencies - assume the high frequency asymptote is the same
dcGainV = max( windData(3).MV(w < 1) );
dcGainW = max( windData(3).MW(w < 1) );

indV = find( windData(3).MV <= (dcGainV-3), 1 );
indW = find( windData(3).MW <= (dcGainW-3), 1 );
fCornerV = f( indV );
fCornerW = f( indW );

% fCornerV = 10^( log10(fCorner2) + -(dcGainV - dcGain)/20);
% fCornerW = 10^( log10(fCorner2) + -(dcGainW - dcGain)/20);
wCornerV = 2*pi*fCornerV;
wCornerW = 2*pi*fCornerW;

% Transfer functions
sysV = frd( windData(3).pV, w, 1/1250 );
sysW = frd( windData(3).pW, w, 1/1250 );

windTFManualV = 10^(dcGainV/20)*tf( wCornerV, [1 wCornerV]);
windTFManualW = 10^(dcGainW/20)*tf( wCornerW, [1 wCornerW]);

% Plot
figure( 'name', 'Magnitude response - All axes' )
semilogx( w, windData(3).MU )
hold on; grid on; box on
semilogx( w, windData(3).MV )
semilogx( w, windData(3).MW )
plot( [wCorner2 windData(3).w(end)], dcGain -20.*( -log10(wCorner2) + log10([wCorner2 windData(3).w(end)])) )

plot( wCorner2, dcGain, 'x' )
plot( wCornerV, dcGainV, '+' )
plot( wCornerW, dcGainW, 'o' )
xlabel( 'Frequency (rad/s)' )
ylabel( 'Magnitude (dB)' )
legend( {'$U$', '$V$', '$W$', '$\omega_{c,U}$', '$\omega_{c,V}$', ...
    '$\omega_{c,W}$'} )
SetFigProp( figSize, fontSize )

% Bode plot
figure( 'name', 'Bode plot - All axes' )
bode( sys, sysV, sysW )
hold on; box on; grid on;
bode( windTFManual1 )
bode( windTFManualV )
bode( windTFManualW )
legend( {'$e_U$', '$e_V$', '$e_W$', '$TF_U$', '$TF_V$', '$TF_W$'} )

SetFigProp( figSize, fontSize )

%% Export
% Get magnitude response
windTFOut = tf(  tf( db2mag(dcGain)*wCorner, [1 wCorner ] ) );
[magExp, phaseExp, w] = bode( sys );
[magFit, phaseFit] = bode( windTFOut, w );
magExp = squeeze( 20*log10(magExp) );
phaseExp = squeeze( phaseExp );
magFit = squeeze( 20*log10(magFit) );
phaseFit = squeeze( phaseFit );

% Plot to verify
semilogx( w, magExp, w, magFit )
axis tight

% Export
fileOut = 'bode_weight_wind.csv';
fid = fopen( fileOut, 'w' );
fprintf( fid, 'w magE phaseE magF phaseF' );
fclose( fid );
writematrix( [w magExp phaseExp magFit phaseFit], fileOut, 'WriteMode', 'append', 'Delimiter', ' ' )
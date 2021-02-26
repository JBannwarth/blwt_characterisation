%GENERATETURBSIMSETTINGS Generate Turbsim settings file from given data.
%
%   Generate a Turbsim settings file that can be used to generate wind
%   profiles with equivalent characteristics to provided wind speed
%   measurement data.
%
%   Written: 2021/02/25, J.X.J. Bannwarth

%% Load data and compute necessary values
clearvars; clc
folderIn  = 'data';
folderOut = 'turbsim';
templateFile = 'turbsim_input_blank.inp';
seedFile = 'turbsim_seeds.mat';
batchFile = 'batch_run.bat';
regenSeed = false;

% Get files
files = dir( fullfile( folderIn, 'coarsegrid_center_*_2.thA' ) );
fileNames = {files.name}';
namesNoExt = erase( fileNames, '.thA' );
inputSettings = split( namesNoExt, '_' );
inputSettings = cellfun( @(x)(str2double(x)), inputSettings(:,3) );

UMean = zeros( length(namesNoExt), 3 );
UStd  = UMean;
UTi   = UMean;
Li    = zeros( length(namesNoExt), 1 );
turbScaleParam = zeros( length(namesNoExt), 1 );
hubHeight      = zeros( length(namesNoExt), 1 );
for ii = 1:length(namesNoExt)
    % Load data from file
    [ u, v, w, Ps, Settings ] = ReadTHFile( fullfile( folderIn, fileNames{ii}) );
    
    % Fix misalignments
    [ u, v, w ] = ConditionWindSpeed( u, v, w, Settings.DataRate );
    
    % Statistics
    U = [u v w];
    UMean(ii,:) = mean(U);
    UStd(ii,:) = std(U);
    UTi(ii,:) = 100.* UStd(ii,:) ./ UMean(ii,1);
    UZeroed = U - UMean(ii,:);
    
    % Sampling parameters
    fs = Settings.DataRate;
    ts = 1 / fs;
    
    % Calculate double-sided spectrum using FFT
    L = size( UZeroed, 1);
    Y = fft( UZeroed, L, 1);
    f = ( 0:(fs/L):(fs/2-fs/L) )';
    S = ( ts / L ) * abs(Y).^2; % [(m/s)^2/Hz]
    
    % Convert to single-sided spectrum (DC and Nyquist frequency occur once
    % only)
    S(1,:) = 0; % DC frequency
    S(2:L/2,:) = 2 * S(2:L/2,:);
    S = S(1:L/2,:);
    
    SNorm = zeros( size( S ) );
    for jj = 1:3
        SNorm(:,jj) = S(:,jj) .* f / UStd(ii,jj)^2;
    end
    
    % Remove DC point
    f = f(2:end);
    SNorm = SNorm(2:end,:);
    
    % Fit theoretical spectrum to obtain length scales
    hCostFun = @(x)WindSpectrumCostFun( x, f, UMean(ii,1), 1.25, ...
        SNorm(:,jj), jj );
    if ii > 1
        L0 = Li(ii-1);
    else
        L0 = 1;
    end
    
    % Fit Li
    optOptions         = optimoptions('lsqnonlin');
    optOptions.Display = 'iter-detailed';
    optOptions.TolFun  = 1e-9;
    optOptions.TolX    = 1e-9;
    
    Li(ii) = lsqnonlin(hCostFun, L0, 0, 100, optOptions);
    
    % Calculate turbulence scale parameter
    turbScaleParam(ii) = Li(ii)/3.5;
    
    % Calculate hub height for TurbSim based on Edition 3
    hubHeight(ii) = turbScaleParam(ii)/0.7;
end

UStdRatio = UStd ./ UStd(:,1);

nWsp      = length( UMean );
nTurbSeed = 6;
UStdNorm = zeros( size(UStdRatio,1)*nTurbSeed, 3 );
for ii = 1:nWsp
    UStdNorm((1+(ii-1)*nTurbSeed):(ii*nTurbSeed), :) = repmat( ...
        UStdRatio(ii,:), nTurbSeed, 1);
end

% Other parameters
usableTime = 700;
timeStep = 0.0008;

% Save seeds for reproducibility          
if ~isfolder( folderOut )
    mkdir( folderOut )
end

if ~isfile( fullfile( folderOut, seedFile ) ) || regenSeed
    seedNb = randi( [-2147483648, 2147483647], nWsp, nTurbSeed );
    save( fullfile( folderOut, seedFile ), 'seedNb')
else
    load( fullfile( folderOut, seedFile ), 'seedNb' )
end

%% Generate input files
% Read in template
templateStr = fileread( templateFile );

% Parameter list
paramName = { '[Seed]';
              '[TimeStep]';
              '[HubHeight]';
              '[TurbInt]';
              '[MeanWind]';
              '[UsableTime]' };

% Format input files
fileStr = cell( nWsp, nTurbSeed );
fileName = cell( nWsp, nTurbSeed );
for ii = 1:nWsp
    for jj = 1:nTurbSeed
        % File name
        fileName{ii,jj} = sprintf( 'turbsim_%02d_%02d.inp', ...
            inputSettings(ii), jj );
        
        % Modify template
        paramStr = { ...
            sprintf( '%d'  , seedNb(ii,jj) );
            sprintf( '%.4f', timeStep      );
            sprintf( '%.8f', hubHeight(ii) );
            sprintf( '%.2f', UTi(ii)       );
            sprintf( '%.2f', UMean(ii)     );
            sprintf( '%d'  , usableTime    ) };
        
        fileStr{ii,jj} = replace( templateStr, paramName, paramStr );
    end
end

% Clear output folder
delete( fullfile( folderOut, '*.inp' ) )

% Write to input files
for ii = 1:nWsp
    for jj = 1:nTurbSeed
        fid = fopen( fullfile( folderOut, fileName{ii,jj}), 'w' );
        fprintf( fid, '%s', fileStr{ii,jj} );
        fclose( fid );
    end
end

%% Generate batch run file
filePaths = dir( fullfile( folderOut, '*.inp' ) );
filePaths = fullfile( {filePaths.folder}', {filePaths.name}' );
batchStr = strjoin( strcat( 'Turbsim', {' '}, filePaths ), '\n' );
fid = fopen( fullfile( folderOut, batchFile ), 'w' );
fprintf( fid, '%s', batchStr );
fclose( fid );
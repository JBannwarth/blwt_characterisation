%GENERATETURBSIMSETTINGS Generate Turbsim settings file from given data.
%
%   Generate a Turbsim settings file that can be used to generate wind
%   profiles with equivalent characteristics to provided wind speed
%   measurement data.
%
%   Written: 2021/02/25, J.X.J. Bannwarth

%% Load data and compute necessary values
clearvars; clc
folderIn = 'data';

% Get files
files = dir( fullfile( folderIn, 'coarsegrid_center_*_2.thA' ) );
fileNames = {files.name}';
fileNames = fileNames( endsWith( fileNames, '.thA' ) );
namesNoExt = erase( fileNames, '.thA' );

UMean = zeros( length(namesNoExt), 3 );
UStd = UMean;
UTi = UMean;
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
    
    % Compute hub height
    fs = Settings.DataRate;
    ts = 1 / fs;
    UZeroed = U - UMean(ii,:);
    
    %
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
for ii = 1:nTurbSeed
    UStdNorm = [ UStdNorm; repmat(UStdRatio, nTurbSeed, 1) ];
end

% Other parameters
usableTime = 700;
timeStep = 0.0008;
% Parameter list
paramName = { '[Seed]';
              '[TimeStep]';
              '[HubHeight]';
              '[TurbInt]';
              '[MeanWind]';
              '[UsableTime]' };

if ~exist( fullfile( tmpFolder, 'turbSimSeeds.mat' ), 'dir' ) || regenSeed
    seedNo = randi([-2147483648,2147483647], nWsp, nTurbSeed);
    save( fullfile( tmpFolder, 'turbSimSeeds.mat' ), 'seedNo')
else
    load( fullfile( tmpFolder, 'turbSimSeeds.mat' ) )
end

%% GENERATE INPUT FILES
% Read in template
fid = fopen( fullfile( turbSimFolder, outputFolder, 'TurbSimInputFileBlank.inp' ), 'r' );
n = 1;
while ~feof(fid)
    try
    data{n, 1} = fgetl(fid);
    catch
    end
    n = n + 1;
end
fclose(fid);

% Write TurbSim input files
n = 1;
fileName{nWsp*nTurbSeed,1} = '';
for ii=1:nWsp
    for j=1:nTurbSeed
        fileName{n} = sprintf( 'TurbSim_%02d_%02d.inp', tfwt(ii), j );
        fid = fopen( fullfile( turbSimFolder, outputFolder, fileName{n} ), 'w' );
        paramString = { sprintf('%d',seedNo(ii,j) );
            sprintf( '%.4f', timeStep );
            sprintf( '%.8f', hubHeight(ii) );
            sprintf( '%.2f', windTI(ii) );
            sprintf( '%.2f', windAvg(ii) );
            sprintf( '%d',   usableTime) };
        for k = 1:length(data)
            tempData = data{k};
            for m = 1:length(paramName)
                tempData = strrep( tempData, paramName{m}, paramString{m} );
            end
            fprintf( fid, '%s\n', tempData );
        end
        fclose( fid );
        n = n + 1;
    end
end
save( fullfile( tmpFolder, 'turbSimFileNames.mat' ), 'fileName', 'windStdNorm' )

%% GENERATE BATCH RUN FILE
fid = fopen( fullfile( turbSimFolder, outputFolder, 'batchRun.bat' ), 'w' );
for ii=1:length(fileName)
    fprintf( fid, '.\\TurbSim .\\%s\\%s\n', outputFolder, fileName{ii} );
end
fclose( fid );
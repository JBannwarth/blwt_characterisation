%LOADDATA2019_08_01
%   Created:       2019/08/01, by J.X.J. Bannwarth
%   Last modified: 2019/08/01, by J.X.J. Bannwarth

%% Load data and compute stats
clearvars;
clc
folderIn = 'data';
files = dir( folderIn );

fileNames = {files.name}';
fileNames = fileNames( endsWith( fileNames, '.thA' ) );
namesNoExt = erase( fileNames, '.thA' );

T = table( 'Size', [0 15], ...
    'VariableTypes', { 'string', 'categorical', 'categorical', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string', 'string', 'double', 'double' },  ...
    'VariableNames', { 'filename', 'grid', 'location', 'inputPct', 'height', 'Umean', 'Uti', 'Vti', 'Wti', 'pressure', 'temperature', 'date', 'time', 'dataRate', 'deviceID'} );
T.Properties.VariableUnits = { '', '', '', 'Pct', '', 'm/s', 'Pct', 'Pct', 'Pct', 'kPa', 'degC', '', '', 'Hz', ''  };

for ii = 1:length(namesNoExt)
    % Load data from file
    [u, v, w, Ps, Settings] = ReadTHFile( fullfile( folderIn, fileNames{ii}) );
    
    [ u, v, w ] = ConditionWindSpeed( u, v, w, Settings.DataRate );
    windData(ii).U = u; windData(ii).V = v; windData(ii).W = w;
    windData(ii).Ps = Ps;
    
    U = [u v w];
    Umean = mean(U);
    Ustd = std(U);
    Uti = 100.* Ustd ./ Umean(1);
    
    C = textscan( namesNoExt{ii}, '%s %s %2d %1d', 'Delimiter', '_' ); % grid location inputPct height
    C{1} = erase( C{1}, 'grid' );
    T = [T; [fileNames{ii}, C, Umean(1), Uti(1), Uti(2), Uti(3), Settings.Pbaro/100, Settings.Tmean, Settings.FirstSampleDate, Settings.FirstSampleTime, Settings.DataRate, Settings.DeviceID]];
end

% PSDs
L1u = 5.67*1.5;
f = logspace(-2,2,501);

% Calculate PSDs
for ii = 1:length(windData)
    [windData(ii).p, windData(ii).f] = pwelch( windData(ii).U-mean(windData(ii).U), ...
        hann(1024, 'periodic'), [], f, Settings.DataRate );
end

%% Set-up Plots
set( groot, 'defaulttextinterpreter', 'latex' ); % necessary to avoid warnings
close all

% Categories
locations = categories( T.location );
grids = categories( T.grid );
grids = grids(end:-1:1, :);
inputPcts = unique( T.inputPct );
heights = unique( T.height );

% Formatting
cMap = lines();
decFactor = 5;
figSize = [10, 15];
fontSize = 12;
lineStyles = { '-', '--', '-.', ':' };
legendStr = arrayfun( @(f)( sprintf( 'in = %d\\%%', f ) ), inputPcts, ...
    'UniformOutput', false );

%% Plot PSDs
% Create figures
hFig = zeros( size(locations) );
hAx  = zeros( length(locations), length(inputPcts) );
for ii = 1:length(locations)
    hFig(ii) = figure( 'Name', sprintf( 'PSD - %s', locations{ii} ) );
    for jj = 1:length(heights)
        hAx(ii,jj) = subplot( 2, 2, jj );
        hold on; grid on; box on;
        xlabel( 'Frequency, $f$ (Hz)' )
        ylabel( '$fS_U/U_\mathrm{mean}^2$ (-)' )
        title( sprintf( '$h$ = %d', heights(jj) ) )
        set( gca, 'XScale', 'log', 'YScale', 'log' )
    end
end

hold on; grid on; box on

% Fill figures
for ii = 1:length(locations)
    for jj = 1:length(heights)
        for kk = 1:length(inputPcts)
            for ll = 1:length(grids)
                ind = (T.location == locations{ii}) & ...
                    (T.height == heights(jj)) & (T.inputPct == inputPcts(kk)) & ...
                    (T.grid == grids{ll});
                ind = find( ind, 1 );
                
                % Plot experimental data
                plot( hAx(ii, jj), windData(ind).f(1:decFactor:end), ...
                    windData(ind).p(1:decFactor:end).*windData(ind).f(1:decFactor:end)/T.Umean(ind)^2, ...
                    lineStyles{kk}, 'Color', cMap(ll,:) )
            end
            
            % Plot theoretical data
            inds = (T.location == locations{ii}) & ...
                    (T.height == heights(jj)) & (T.inputPct == inputPcts(kk));
            Umean = mean( T.Umean(inds) );
            SNorm = CalcVonKarmanIEC( f, Umean );
            plot( hAx(ii, jj), f, SNorm(:,1), lineStyles{kk}, ...
                'color', 'k' )
        end
    end
end

% Finishing touches
for ii = 1:length(hFig)
    figure( hFig(ii) )
    legend( hAx(ii,end), legendStr, 'location', 'best' )
    SetFigProp( figSize, fontSize )
end

%% Plot velocities
figure;
a1 = axes; hold on; grid on; box on;
colours = get(gca,'colororder');
linestyles = { '-', '-.', '--', ':' };
for j = 1:length(grids)
    for ii = 1:length(locations)
        toPlot = (T.grid == grids{j}) & (T.location == locations{ii}) & (T.height == 2);
        plot( T.inputPct(toPlot), T.Uti(toPlot), "color", colours(ii,:), "LineStyle", linestyles{j} );
    end
end
xlabel( 'Wind tunnel input (\%)')
ylabel( 'Turbulence intensity (\%)')

a2 = axes( 'position', get(gca,'position'), 'visible', 'off'); hold on;
h(1) = plot(NaN,NaN,'k-');
h(2) = plot(NaN,NaN,'k-.');
h(3) = plot(NaN,NaN,'k--');
legend( a2, h, { 'No grid', 'Fine grid', 'Coarse grid' }, ...
    'location', "eastoutside" )
axes(a1)
legend( locations, 'Location', "northeastoutside" )
SetFigProp( figSize, fontSize )
a2.Position = a1.Position;
MatlabToLatexEps( 'TIvsInput' )

% Wind speed
locations = categories( T.location );
inputPcts = unique(T.inputPct);
grids = categories( T.grid );
grids = grids(end:-1:1, :);
figure;
a1 = axes; hold on; grid on; box on;
colours = get(gca,'colororder');
linestyles = { '-', '-.', '--', ':' };
for j = 1:length(grids)
    for ii = 1:length(inputPcts)
        toPlot = (T.grid == grids{j}) & (T.inputPct == inputPcts(ii)) & (T.height == 2);
        plot( T.location(toPlot), T.Umean(toPlot), "color", colours(ii,:), "LineStyle", linestyles{j} );
    end
end
xlabel( 'Location')
ylabel( 'Mean wind speed (m/s)')

a2 = axes( 'position', get(gca,'position'), 'visible', 'off'); hold on;
h(1) = plot(NaN,NaN,'k-');
h(2) = plot(NaN,NaN,'k-.');
h(3) = plot(NaN,NaN,'k--');
legend( a2, h, { 'No grid', 'Fine grid', 'Coarse grid' }, ...
    'location', "eastoutside" )
axes(a1)
legendStr = strsplit( sprintf( 'Input = %d\\%%\n', inputPcts ), '\n' );
legendStr = legendStr(~cellfun(@isempty, legendStr));
legend( legendStr, 'Location', "northeastoutside" )
SetFigProp( figSize, fontSize )
a2.Position = a1.Position;
MatlabToLatexEps( 'UmeanvsInput' )
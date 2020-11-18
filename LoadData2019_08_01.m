%LOADDATA2019_08_01 Load and process experimental data
%   Created:       2019/08/01, by J.X.J. Bannwarth
%   Last modified: 2020/11/18, by J.X.J. Bannwarth

%% Load data and compute stats
clearvars; clc
folderIn = 'data';

% Get files
files = dir( folderIn );
fileNames = {files.name}';
fileNames = fileNames( endsWith( fileNames, '.thA' ) );
namesNoExt = erase( fileNames, '.thA' );

% Allocate results table
T = table( 'Size', [0 15], ...
    'VariableTypes', { 'string', 'categorical', 'categorical', ...
        'double', 'double', 'double', 'double', 'double', 'double', ...
        'double', 'double', 'string', 'string', 'double', 'double' }, ...
    'VariableNames', { 'filename', 'grid', 'location', 'inputPct', ...
        'height', 'Umean', 'Uti', 'Vti', 'Wti', 'pressure', ...
        'temperature', 'date', 'time', 'dataRate', 'deviceID' } );
T.Properties.VariableUnits = { '', '', '', '%', '', 'm/s', '%', ...
    '%', '%', 'kPa', 'degC', '', '', 'Hz', '' };

% Probe heights
h = [0.1, 0.75, 1.25, 1.5]; % (m)

% Load data into table
for ii = 1:length(namesNoExt)
    % Load data from file
    [ u, v, w, Ps, Settings ] = ReadTHFile( fullfile( folderIn, fileNames{ii}) );
    
    % Fix misalignments
    [ u, v, w ] = ConditionWindSpeed( u, v, w, Settings.DataRate );
    windData(ii).U = u; windData(ii).V = v; windData(ii).W = w;
    windData(ii).Ps = Ps;
    
    % Compute statistics
    U = [u v w];
    Umean = mean(U);
    Ustd = std(U);
    Uti = 100.* Ustd ./ Umean(1);
    
    % Filename format: <grid>_<location>_<inputPct>_<height>
    C = textscan( namesNoExt{ii}, '%s %s %2d %1d', 'Delimiter', '_' );
    C{1} = erase( C{1}, 'grid' );
    C{1} = replace( C{1}, 'center', 'centre' ); % Change to UK spelling
    C{end} = h(C{end}+1); % Get corresponding height
    
    % Fill table
    T = [ T; [fileNames{ii}, C, Umean(1), Uti(1), Uti(2), Uti(3), ...
        Settings.Pbaro/100, Settings.Tmean, Settings.FirstSampleDate, ...
        Settings.FirstSampleTime, Settings.DataRate, Settings.DeviceID] ];
end

% Power spectral density (PSD) setup
normalisePSD = true;
f = logspace(-2,2,501);
if normalisePSD
    post = '';
else
    post = '_no_norm';
end

% Calculate PSDs
for ii = 1:length(windData)
    [windData(ii).p, windData(ii).f] = pwelch( ...
        windData(ii).U-mean(windData(ii).U), hann(1024, 'periodic'), ...
        [], f, Settings.DataRate );
end

%% Set-up Plots
% Set default interpreter to avoid symbol warnings
set( groot, 'defaultTextInterpreter', 'latex' );
set( groot, 'defaultLegendInterpreter', 'latex' );
close all

% Categories
locations = categories( T.location );
grids = categories( T.grid );
grids = grids(end:-1:1, :);
inputPcts = unique( T.inputPct );
heights = unique( T.height );

% Formatting
colours = lines();
decFactor = 5;
lineStyles = { '-', '--', '-.', ':' };
gridLegend = cellfun( @(f)( sprintf( '%s grid', f ) ), grids, ...
    'UniformOutput', false );

% Export
figSize = [15, 10];
figSizeLarge = [15, 15];
fontSize = 11;
folderOut = 'figures';

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
        if normalisePSD
            ylabel( '$fS_U/U_\mathrm{mean}^2$ (-)' )
        else
            ylabel( '$S_U/U_\mathrm{mean}^2$ (-)' )
        end
        xlim( [1E-2 1E2] )
        ylim( [1E-8 1E2] )
        title( sprintf( '$h$ = %.2f m', heights(jj) ) )
        set( gca, 'XScale', 'log', 'YScale', 'log' )
    end
end

hold on; grid on; box on

% Fill figures
for ii = 1:length(locations)
    for jj = 1:length(heights)
        for kk = 1:length(inputPcts)
            for ll = 1:length(grids)
                % Get corresponding index from table
                ind = (T.location == locations{ii}) & ...
                    (T.height == heights(jj)) & ...
                    (T.inputPct == inputPcts(kk)) & (T.grid == grids{ll});
                ind = find( ind, 1 );
                
                % Plot experimental data
                if normalisePSD
                    SNormExp = windData(ind).p(1:decFactor:end) .* ...
                        windData(ind).f(1:decFactor:end) ./ T.Umean(ind)^2;
                else
                    SNormExp = windData(ind).p(1:decFactor:end) ./ ...
                        T.Umean(ind)^2;
                end
                
                plot( hAx(ii, jj), windData(ind).f(1:decFactor:end), ...
                    SNormExp, lineStyles{kk}, 'Color', colours(ll,:) )
            end
            
            % Plot theoretical data
            inds = (T.location == locations{ii}) & ...
                    (T.height == heights(jj)) & (T.inputPct == inputPcts(kk));
            Umean = mean( T.Umean(inds) );
            SNorm = CalcVonKarmanIEC( f, Umean, 1.25, normalisePSD );
            plot( hAx(ii, jj), f, SNorm(:,1), lineStyles{kk}, ...
                'color', 'k' )
        end
    end
end

% Finishing touches
for ii = 1:length(hFig)
    figure( hFig(ii) )
    legend( hAx(ii,end), [gridLegend; 'von K{\''a}rm{\''a}n'], ...
        'location', 'northwest' )
    SetFigProp( figSizeLarge, fontSize )
    MatlabToLatexEps( fullfile( folderOut, ...
        sprintf('psd_%s%s', locations{ii}, post) ), 300 )
end

%% Plot velocities at mid altitude
% Turbulence internsity vs wind tunnel dial setting
figure( 'Name', 'TI vs input' );
hAx1(1) = axes; hold on; grid on; box on;
linestyles = { '-', '-.', '--', ':' };

% Plot data
for j = 1:length(grids)
    for ii = 1:length(locations)
        toPlot = (T.grid == grids{j}) & ...
            (T.location == locations{ii}) & (T.height == 1.25);
        plot(hAx1(1), T.inputPct(toPlot), T.Uti(toPlot), "color", ...
            colours(ii,:), "LineStyle", linestyles{j} );
    end
end
xlabel( 'Wind tunnel fan setting (\%)')
ylabel( 'Turbulence intensity (\%)')

% Add second legend for line styles
hAx1(2) = axes( 'position', get(hAx1(1),'position'), 'visible', 'off');
hold on;
plot( hAx1(2), NaN, NaN, 'k-' );
plot( hAx1(2), NaN, NaN, 'k-.' );
plot( hAx1(2), NaN, NaN, 'k--' );
legend( hAx1(2), { 'No grid', 'Fine grid', 'Coarse grid' }, ...
    'location', "eastoutside" )

% Main legend
axes(hAx1(1))
legend( replace( locations, 'center', 'centre' ), 'Location', ...
    'northeastoutside' )

% Export
SetFigProp( figSize, fontSize )
hAx1(2).Position = hAx1(1).Position;
MatlabToLatexEps( fullfile( folderOut, 'ti_vs_input' ), 300 )

% Mean wind speed vs location
figure( 'name', 'U vs location' );
hAx2(1) = axes; hold on; grid on; box on;
for j = 1:length(grids)
    for ii = 1:length(inputPcts)
        toPlot = (T.grid == grids{j}) & ...
            (T.inputPct == inputPcts(ii)) & (T.height == 1.25);
        plot( T.location(toPlot), T.Umean(toPlot), "color", ...
            colours(ii,:), "LineStyle", linestyles{j} );
    end
end
xlabel( 'Location')
ylabel( 'Mean wind speed (m/s)')
ylim( [0 15] )

% Add second legend for line styles
hAx2(2) = axes( 'position', get(hAx2(1),'position'), 'visible', 'off');
hold on;
plot( hAx2(2), NaN, NaN, 'k-' );
plot( hAx2(2), NaN, NaN, 'k-.' );
plot( hAx2(2), NaN, NaN, 'k--' );
legend( hAx2(2), { 'No grid', 'Fine grid', 'Coarse grid' }, ...
    'location', "eastoutside" )
axes(hAx2(1))
legendStr = strsplit( sprintf( 'Input = %d\\%%\n', inputPcts ), '\n' );
legendStr = legendStr(~cellfun(@isempty, legendStr));
legend( legendStr, 'Location', "northeastoutside" )

% Export
SetFigProp( figSize, fontSize )
hAx2(2).Position = hAx2(1).Position;
MatlabToLatexEps( fullfile( folderOut, 'windspeed_vs_location' ), 300 )
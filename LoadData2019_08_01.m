%LOADDATA2019_08_01
%   Created:       2019/08/01, by J.X.J. Bannwarth
%   Last modified: 2019/08/01, by J.X.J. Bannwarth

%% Load data and compute stats
folderIn = 'data';
files = dir( folderIn );

fileNames = {files.name}';
fileNames = fileNames( endsWith( fileNames, '.thA' ) );
namesNoExt = erase( fileNames, '.thA' );

T = table( 'Size', [0 15], ...
    'VariableTypes', { 'string', 'categorical', 'categorical', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string', 'string', 'double', 'double' },  ...
    'VariableNames', { 'filename', 'grid', 'location', 'inputPct', 'height', 'Umean', 'Uti', 'Vti', 'Wti', 'pressure', 'temperature', 'date', 'time', 'dataRate', 'deviceID'} );
T.Properties.VariableUnits = { '', '', '', 'Pct', '', 'm/s', 'Pct', 'Pct', 'Pct', 'kPa', 'degC', '', '', 'Hz', ''  };

for i = 1:length(namesNoExt)
    % Load data from file
    [u, v, w, Ps, Settings] = ReadTHFile( fullfile( folderIn, fileNames{i}) );
    
    [ u, v, w ] = ConditionWindSpeed( u, v, w, Settings.DataRate );
    windData(i).U = u; windData(i).V = v; windData(i).W = w;
    windData(i).Ps = Ps;
    
    U = [u v w];
    Umean = mean(U);
    Ustd = std(U);
    Uti = 100.* Ustd ./ Umean(1);
    
    C = textscan( namesNoExt{i}, '%s %s %2d %1d', 'Delimiter', '_' ); % grid location inputPct height
    C{1} = erase( C{1}, 'grid' );
    T = [T; [fileNames{i}, C, Umean(1), Uti(1), Uti(2), Uti(3), Settings.Pbaro/100, Settings.Tmean, Settings.FirstSampleDate, Settings.FirstSampleTime, Settings.DataRate, Settings.DeviceID]];
end

%% Plot results
set(groot,'defaulttextinterpreter','latex'); % necessary to avoid warnings

close all
locations = categories( T.location );
grids = categories( T.grid );
grids = grids(end:-1:1, :);
figure;
a1 = axes; hold on; grid on; box on;
colours = get(gca,'colororder');
linestyles = { '-', '-.', '--', ':' };
for j = 1:length(grids)
    for i = 1:length(locations)
        toPlot = (T.grid == grids{j}) & (T.location == locations{i}) & (T.height == 2);
        plot( T.inputPct(toPlot), T.Uti(toPlot), "color", colours(i,:), "LineStyle", linestyles{j} );
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
SetFigProp( [15,10], 11 )
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
    for i = 1:length(inputPcts)
        toPlot = (T.grid == grids{j}) & (T.inputPct == inputPcts(i)) & (T.height == 2);
        plot( T.location(toPlot), T.Umean(toPlot), "color", colours(i,:), "LineStyle", linestyles{j} );
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
SetFigProp( [15,10], 11 )
a2.Position = a1.Position;
MatlabToLatexEps( 'UmeanvsInput' )
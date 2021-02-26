function cost = WindSpectrumCostFun( L, n, U, z, dataExp, ii )
%WINDSPECTRUMCOSTFUN Return cost for difference between exp and sim spectra
%
%   Inputs:
%       - L:       Length scale to compute theoretical spectrum for (m)
%       - n:       Frequencies to compute the spectra at (Hz)
%       - U:       Mean wind speed to compute theoretical spectrum for (m/s)
%       - z:       Height above the ground (m)
%       - dataExp: Experimental wind spectrum data
%       - ii:      Index of axis to compute cost for (1=u, 2=v, 3=w)
%   Output:
%       - cost:    Cost (-)
%
%   See also GENERATETURBSIMSETTINGS.
%
%   Written: 2017, Z.J. Chen

    % Calculate theoretical spectrum
    dataTheoretical = CalcVonKarmanIEC( n, U, z, true, [L;L;L] );
    dataTheoretical = dataTheoretical(:, ii);
    
    % Calculate cumulative spectra
    dataExp = cumtrapz( n, dataExp./n );
    dataTheoretical = cumtrapz( n, dataTheoretical./n );
    
    % Find frequency over which to weight spectrum
    [ ~, startInd ] = min( abs( dataExp - 0.1 ) );
    [ ~, endInd   ] = min( abs( dataExp - 0.9 ) );
    
    % Calculate cost function weights
    weighting = zeros( length(n), 1 );
    weighting(startInd:endInd) = 1;
    weighting(startInd:endInd) = weighting(startInd:endInd) ./ ...
        dataExp(startInd:endInd);
    
    % Calculate cost function
    cost = (dataExp - dataTheoretical) .* weighting;
end
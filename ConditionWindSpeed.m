function [ Uout, Vout, Wout ] = ConditionWindSpeed( U, V, W, dataRate )
%CONDITIONWINDSPEED Remove bad data and correct for orientation
%   Created:       2019/08/01, by J.X.J. Bannwarth
%   Last modified: 2019/08/01, by J.X.J. Bannwarth

    ns = length(U);
    ts = 1/dataRate;
    t = (0:ts:(ns-1)*ts)';
    goodDataInd = U>0;
    goodDataPct = sum(goodDataInd)/ns*100;
    tGoodData = t(goodDataInd);
    goodData = U(goodDataInd);
    U = interp1(tGoodData,goodData,t,'linear','extrap');
    goodData = V(goodDataInd);
    V = interp1(tGoodData,goodData,t,'linear','extrap');
    goodData = W(goodDataInd);
    W = interp1(tGoodData,goodData,t,'linear','extrap');
    
    meanUvec = [ mean(U); mean(V); mean(W) ];
    
    v = meanUvec./norm(meanUvec); % Use z axis as ref
    
    az = atan2( v(2), v(1));
    R1 = [ cos(-az), -sin(-az), 0;
           sin(-az), cos(-az),   0;
           0,       0,        1];
    v2 = R1 * v;
    el = atan2( v2(3), v2(1) );
    R2 = [ cos(el), 0, sin(el);
           0,       1, 0;
          -sin(el), 0, cos(el)];
    v3 = R2 *v2;
    Rcor = R2*R1;
    
    dataTmp = Rcor*[ U'; V'; W' ];
    Uout = dataTmp(1,:)';
    Vout = dataTmp(2,:)';
    Wout = dataTmp(3,:)';
end


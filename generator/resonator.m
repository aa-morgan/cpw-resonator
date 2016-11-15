% Provided parameters
  % Resonator 
resRatio    = [9.065/5.0, 9.065/5.0];
resGap      = [5.0, 5.0];
resCond     = resGap.*resRatio;
innerRadius = [40, 40];
steps       = [54, 54];
    % ---
couplerLen  = [120, 120];
intLen      = [5000, 5000];
intDist     = [1300, 1300];
totalLen    = [7624, 7624]; % Only applicable if midCurves=true
    % ---
midCurves   = [true, true];
vertRatio   = [0.3, 0.3];
  % Feedline
feedRatio       = 9.065/5.0;
feedGap         = 10.0;                       
feedCond        = feedGap*feedRatio;
feedDist        = [20, 20];
feedClearanceX  = [100, 100];
feedLen         = 500; % Only applicable if doubleRes=false
  % Multi resonator
doubleRes   = true;
resDist     = 5000;
  % CleWin
cleWin          = false;
cleWinVerbose   = false;

% Start figure
if not(cleWin) figure('Position', [0, 0, 1800, 600]); end;

% Store polygons
polygons = {};

% Start Object generation
if doubleRes numRes=2; else numRes=1; end;
for k = 1:numRes
    % Polygons index;
    m = 1;
    
    % Derived parameters
    innerRadiusInnerGap=innerRadius(k);
    outerRadiusInnerGap=innerRadiusInnerGap+resGap(k);
    innerRadiusOuterGap=outerRadiusInnerGap+resCond(k);
    outerRadiusOuterGap=innerRadiusOuterGap+resGap(k);
    tile = outerRadiusOuterGap;
    tileShort = tile-innerRadiusInnerGap;
    tileLong  = tile-tileShort;
    
    if midCurves(k) numCurved = 6; else numCurved = 2;end;

    vert1Len=(vertRatio(k))*(intDist(k)-(feedDist(k)+(tile+tileLong)+(midCurves(k)*2*(tile+tileLong))));
    vert2Len=(1-vertRatio(k))*(intDist(k)-(feedDist(k)+(tile+tileLong)+(midCurves(k)*2*(tile+tileLong))));
    vertLen = [vert1Len,vert2Len];

    curvedPathLen=0.5*pi*((outerRadiusOuterGap+innerRadiusInnerGap)/2);

    if midCurves(k)
        remainLen = totalLen(k) - (numCurved*curvedPathLen+vert1Len+vert2Len+couplerLen(k)+intLen(k));
        horzLen=remainLen/2;
        mTotalLen(k) = totalLen(k);

        if remainLen <= 0
            error('Error: Sum of lengths are not long enough to satisfy the total length criteria!')
        end
    else
        horzLen = 0;   
        mTotalLen(k) = numCurved*curvedPathLen+vert1Len+vert2Len+couplerLen(k)+intLen(k);
    end
    
    % Overwrite feedline lengths if doubleRed is on
    if doubleRes
        feedLen = ((resDist/2)+feedClearanceX(k)+tile);
    end

    % Generate curved components
    for i = 0:3
        start=i*0.5*pi;
        angles = linspace(0,0.25*2*pi,steps(k));
        pt =[cos(angles-start);sin(angles-start)];

        rInner1=outerRadiusInnerGap*pt;
        rInner2=innerRadiusInnerGap*pt;
        rInner=cat(2,rInner1,flip(rInner2,2));

        rOuter1=outerRadiusOuterGap*pt;
        rOuter2=innerRadiusOuterGap*pt;
        rOuter=cat(2,rOuter1,flip(rOuter2,2));

        if start == 0.0*pi % UR
            quarterURinnerGap = rInner;
            quarterURouterGap = rOuter;
        elseif start == 0.5*pi % LR
            quarterLRinnerGap = [rInner(1,:);rInner(2,:)+outerRadiusOuterGap];
            quarterLRouterGap = [rOuter(1,:);rOuter(2,:)+outerRadiusOuterGap];   
        elseif start == 1.0*pi % LL
            quarterLLinnerGap = [rInner(1,:)+outerRadiusOuterGap;rInner(2,:)+outerRadiusOuterGap];
            quarterLLouterGap = [rOuter(1,:)+outerRadiusOuterGap;rOuter(2,:)+outerRadiusOuterGap];
        elseif start == 1.5*pi % UL
            quarterULinnerGap = [rInner(1,:)+outerRadiusOuterGap;rInner(2,:)];
            quarterULouterGap = [rOuter(1,:)+outerRadiusOuterGap;rOuter(2,:)];
        end

    end
      % Combine
    quarterInner = [quarterURinnerGap;quarterLRinnerGap;quarterLLinnerGap;quarterULinnerGap];
    quarterOuter = [quarterURouterGap;quarterLRouterGap;quarterLLouterGap;quarterULouterGap];

    % Curved sections
    order = [2,1,3,4,2,1];
    xOff = [0,0,-horzLen-tile,-horzLen-tile,0,0];
    yOff = [0,vert1Len+tile,vert1Len+tile+tileLong, ...
        vert1Len+2*tile+tileLong,vert1Len+2*(tile+tileLong), ...
        vert1Len+vert2Len+midCurves(k)*2*(tile+tileLong)+tile];
    for i = 1:length(order)
        if midCurves(k) || i == 1 || i == 6
            indices = [2*(order(i)-1)+1:2*(order(i)-1)+2];
            polygons{k,m} = [quarterInner(indices(1),:,:)+xOff(i);quarterInner(indices(2),:,:)+yOff(i)]'; m=m+1;
            polygons{k,m} = [quarterOuter(indices(1),:,:)+xOff(i);quarterOuter(indices(2),:,:)+yOff(i)]'; m=m+1;
        end
    end

    % Curve connectors
      % Vert's
    for i = 1:2  
        xOff=tileLong;
        yOff=tile+((i-1)*(vertLen(1)+midCurves(k)*2*(tile+tileLong)));
        vertLower = [[xOff+0,yOff+0]; ...
                      [xOff+0,yOff+vertLen(i)]; ...
                      [xOff+resGap(k),yOff+vertLen(i)]; ...
                      [xOff+resGap(k),yOff+0]]';
        vertUpper = [vertLower(1,:)+resGap(k)+resCond(k);vertLower(2,:)];
        polygons{k,m} = [vertLower(1,:);vertLower(2,:)]'; m=m+1;
        polygons{k,m} = [vertUpper(1,:);vertUpper(2,:)]'; m=m+1;
    end

      % Horz's
    if midCurves(k)
        for i = 0:1 
            xOff=0;
            yOff=tile+vert1Len+tileLong+(i*(tile+tileLong));
            horzLower = [[xOff+0,yOff+0]; ...
                          [xOff-horzLen,yOff+0]; ...
                          [xOff-horzLen,yOff+resGap(k)]; ...
                          [xOff+0,yOff+resGap(k)]]';
            horzUpper = [horzLower(1,:);horzLower(2,:)+resGap(k)+resCond(k)];
            polygons{k,m} = [horzLower(1,:);horzLower(2,:)]'; m=m+1;
            polygons{k,m} = [horzUpper(1,:);horzUpper(2,:)]'; m=m+1;
        end
    end

    % Coupler
    xOff=0;
    yOff=0;
    couplerLower = [[xOff+0,yOff+0]; ...
                    [xOff-couplerLen(k),yOff+0]; ...
                    [xOff-couplerLen(k),yOff+resGap(k)]; ...
                    [xOff+0,yOff+resGap(k)]]';
    couplerUpper = [couplerLower(1,:);couplerLower(2,:)+resGap(k)+resCond(k)];
    xOff=-couplerLen(k);
    yOff=0;
    couplerEnd = [[xOff+0,yOff+0]; ...
                  [xOff-resGap(k),yOff+0]; ...
                  [xOff-resGap(k),yOff+resGap(k)+resCond(k)+resGap(k)]; ...
                  [xOff+0,yOff+resGap(k)+resCond(k)+resGap(k)]]';       
    polygons{k,m} = [couplerLower(1,:);couplerLower(2,:)]'; m=m+1;
    polygons{k,m} = [couplerUpper(1,:);couplerUpper(2,:)]'; m=m+1;
    polygons{k,m} = [couplerEnd(1,:);couplerEnd(2,:)]'; m=m+1;

    % Interaction region
    xOff=0;
    yOff=vert1Len+vert2Len+(tileLong+tile)+midCurves(k)*2*(tileLong+tile);
    horz2Lower = [[xOff+0,yOff+0]; ...
                  [xOff-intLen(k),yOff+0]; ...
                  [xOff-intLen(k),yOff+resGap(k)]; ...
                  [xOff+0,yOff+resGap(k)]]';
    horz2Upper = [horz2Lower(1,:);horz2Lower(2,:)+resGap(k)+resCond(k)];
    polygons{k,m} = [horz2Lower(1,:);horz2Lower(2,:)]'; m=m+1;
    polygons{k,m} = [horz2Upper(1,:);horz2Upper(2,:)]'; m=m+1;

    % Feedline
    xOff=tile+feedClearanceX(k);
    yOff=-(feedDist(k)+feedGap+feedCond+feedGap);
    feedLower = [[xOff+0,yOff+0]; ...
                 [xOff-feedLen,yOff+0]; ...
                 [xOff-feedLen,yOff+feedGap]; ...
                 [xOff+0,yOff+feedGap]]';
    feedUpper = [feedLower(1,:);feedLower(2,:)+feedGap+feedCond];
    polygons{k,m} = [feedLower(1,:);feedLower(2,:)]'; m=m+1;
    polygons{k,m} = [feedUpper(1,:);feedUpper(2,:)]'; m=m+1;

end

% Add all polygons to figure/CleWin
for obj = 1:length(polygons(:,1))
    % Remove empty elements
    R=polygons(obj,:);
    tPolygons = R(~cellfun('isempty',R));
    for i = 1:length(tPolygons)
        data = tPolygons{i};
        x = data(:,1)';
        y = data(:,2)';
        
        % Rotate if required
        if doubleRes && obj == 2
            % create a matrix of these points, which will be useful in future calculations
            v = [x;y];
            % choose a point which will be the center of rotation
            x_center = 0;
            y_center = 0;
            % create a matrix which will be used later in calculations
            center = repmat([x_center; y_center], 1, length(x));
            % define a 60 degree counter-clockwise rotation matrix
            theta = pi;       % pi/3 radians = 60 degrees
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            % do the rotation...
            s = v - center;     % shift points in the plane so that the center of rotation is at the origin
            so = R*s;           % apply the rotation about the origin
            vo = so + center;   % shift again so the origin goes back to the desired center of rotation
            % this can be done in one line as:
            % vo = R*(v - center) + center
            % pick out the vectors of rotated x- and y-data
            x = vo(1,:);
            y = vo(2,:);
            % Apply X,Y shifts
            x = x - resDist;
            y = y - (sum(feedDist)+2*feedGap+feedCond);
        end
        
        if not(cleWin)
            hold on;
            if length(data) == 4
                plot(x, y, 'x-');  
            else
                plot(x, y, '.');  
            end
        else
            polygon([x;y]')
        end
    end
end

% Display Info
message = sprintf(strcat(...
        'Info, (units: micrometers)', '\n', ...
        '\tTotal length:\t\t', num2str(mTotalLen), '\n', ...
        '\tCoupler length:\t\t', num2str(couplerLen), '\n', ...
        '\tInteraction length:\t\t', num2str(intLen), '\n', ...
        '\tFeedline-coupler distance:\t', num2str(feedDist), '\n', ...
        '\tFeedline-interaction distance:\t', num2str(intDist), '\n', ...
        '\tLength extension section:\t', num2str(midCurves), '\n' ...
        ));
if not(cleWin)
    disp(message)
end
if cleWin && cleWinVerbose
    error(message)
end
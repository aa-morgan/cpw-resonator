% ------------------------------------------------------------------
% ----------- Parameter set START ----------------------------------
% ------------------------------------------------------------------
pMap= containers.Map();
  % Resonator 
pMap('Resonator conductor:gap ratio')           = [9.065/5.0, 9.065/5.0];
pMap('Resonator gap width')                     = [10.0, 10.0];
pMap('Resonator conductor width')               = pMap('Resonator gap width').*pMap('Resonator conductor:gap ratio');
pMap('Resonator curve inner radius')            = [40, 40];
pMap('Resonator curve resolution')              = [54, 54];
    % ---
pMap('Resonator coupler length')                = [120, 120];
pMap('Resonator interaction length')            = [5000,5000];%[[100:500:4500];[100:500:4500]]';
pMap('Resonator feedline-interaction distance') = [2000, 2000];
pMap('Resonator total length')                  = [7624, 7624]; % Only applicable if using extendable sections
    % ---
pMap('Use extenable section')                   = [false, false];
pMap('Position of extendable section')          = [0.3, 0.3];
  % Feedline
pMap('Feedline conductor:gap ratio:')           = 9.065/5.0;
pMap('Feedline gap width')                      = 10.0;                       
pMap('Feedline conductor width')                = pMap('Feedline gap width')*pMap('Feedline conductor:gap ratio:');
pMap('Feedline-coupler distance')               = [20, 20];
pMap('Feedline-coupler X clearance')            = [300, 300];
pMap('Feedline length')                         = 800;  % Only applicable if not using two resonators
  % End connector
pMap('Use end connectors')                      = true;
pMap('End connector ratio factor')              = 1;%1.1031;
pMap('End connector pad conductor width')       = 300;
pMap('End connector pad conductor length')      = 400;
pMap('End connector pad gap length')            = 150;
pMap('End connector expansion length')          = 610;
  % Multi resonator
pMap('Use two resonators')                      = true;
pMap('Resonator X separation')                  = 5000;
  % Labels
pMap('Use labels')                              = true;
pMap('Labels')                                  = {'A1','A2','A3','A4'};
pMap('Labels offset ([x,y])')                   = [0,0];
pMap('Labels size')                             = 0.1;
  % Tiling
pMap('Number of tiles up')                      = 2;
pMap('Number of tiles right')                   = 2;
pMap('Tile size')                               = [9000,5700];
  % Alignment marks
pMap('Use alignment marks')                     = false;
pMap('Alignment mark size [length, width]')     = [100,10];
pMap('Alignment offset [x, y]')                 = [0,0];
  % Dicing boundaries
pMap('Use dicing boundaries')                   = true;
pMap('Wafer size (inch)')                       = 3;
pMap('Wafer clearance (mm)')                    = (12.7-10);
pMap('Dicing boundary width')                   = 500;
pMap('Dicing boundary offset [x, y]')           = [0,0];
pMap('Dicing tiling offset [x, y]')             = [3,6];
pMap('Use wafer template')                      = false;
  % Substrate
pMap('Use global substrate')                    = false;
pMap('Substrate padding ([b,t,l,r])')           = [1,1,1,1]*0;
  % CleWin
pMap('Using CleWin')                            = false;
  % General
pMap('Print information')                       = false;

% ------------------------------------------------------------------
% ----------- Run once ---------------------------------------------
% ------------------------------------------------------------------

% Single set parameters
  % Tiling
numTileX = pMap('Number of tiles right');
numTileY = pMap('Number of tiles up');
tileSize = pMap('Tile size');
  % Alignment marks
alignMarks       = pMap('Use alignment marks');
alignMarksSize   = pMap('Alignment mark size [length, width]');
alignMarksOffset = pMap('Alignment offset [x, y]');
  % Dicing boundaries
dicingBound        = pMap('Use dicing boundaries');
waferSize          = pMap('Wafer size (inch)');
waferClearance     = pMap('Wafer clearance (mm)');
dicingBoundWidth   = pMap('Dicing boundary width');
dicingBoundOffset  = pMap('Dicing boundary offset [x, y]');
dicingTilingOffset = pMap('Dicing tiling offset [x, y]');
showWafer          = pMap('Use wafer template');
  % Substrate
globalSub  = pMap('Use global substrate');
subPadding = pMap('Substrate padding ([b,t,l,r])');
  % CleWin
cleWin  = pMap('Using CleWin');
  % General
verbose = pMap('Print information');

% Start figure
if not(cleWin) figure('Position', [0, 0, 1800, 600]); end;

% Run once (regardless of the number of tiles)
  % Global Substrate
if globalSub
    xOff=-7000;
    yOff=-2000;
    globalSubCorners = [[0-subPadding(3),0-subPadding(1)]; ...
                       [(tileSize(1)*numTileX)+subPadding(4),0-subPadding(1)]; ...
                       [(tileSize(1)*numTileX)+subPadding(4),(tileSize(2)*numTileY)+subPadding(2)]; ...
                       [0-subPadding(1),(tileSize(2)*numTileY)+subPadding(2)]];
    data = [globalSubCorners(:,1)+xOff, globalSubCorners(:,2)+yOff];
   
    if not(cleWin)
        hold on;
        plot(data(:,1), data(:,2), '.-');  
    else
        polygon([data(:,1),data(:,2)])
    end
end

  % Generate alignment mark
if alignMarks
        % Points
    point1  = [0,0];
    point2  = [-alignMarksSize(1),0];
    point3  = [-alignMarksSize(1),alignMarksSize(2)];
    point4  = [0,alignMarksSize(2)];
    point5  = [0,sum(alignMarksSize)];
    point6  = [alignMarksSize(2),sum(alignMarksSize)];
    point7  = [alignMarksSize(2),alignMarksSize(2)];
    point8  = [sum(alignMarksSize),alignMarksSize(2)];
    point9  = [sum(alignMarksSize),0];
    point10 = [alignMarksSize(2),0];
    point11 = [alignMarksSize(2),-alignMarksSize(1)];
    point12 = [0,-alignMarksSize(1)];
        % Combine points
    alignMarkPoints = [point1;point2;point3;point4;point5;point6; ...
                       point7;point8;point9;point10;point11;point12];
        % Shift to center
    alignMarkPoints =  [alignMarkPoints(:,1)-alignMarksSize(2)/2, ...
                        alignMarkPoints(:,2)-alignMarksSize(2)/2];  
end           
                
  % Dicing boundaries
if dicingBound
    dicingGridSize = ((waferSize*25.4)+(2*waferClearance))*1000;
    numDiceBoundX = floor(dicingGridSize/tileSize(1))+1;
    numDiceBoundY = floor(dicingGridSize/tileSize(2))+1;
    
    xOff=alignMarksOffset(1)-7000;
    yOff=alignMarksOffset(2)-tileSize(2)/2;

    % Store polygons.
    polygons = {};
    
    for i = 1:numDiceBoundX
        point1 = [dicingBoundWidth+((i-1)*tileSize(1)),0];
        point2 = [dicingBoundWidth+((i-1)*tileSize(1)),dicingGridSize];
        point3 = [0+((i-1)*tileSize(1)),dicingGridSize];
        point4 = [0+((i-1)*tileSize(1)),0];
        polygons{end+1} = [point1;point2;point3;point4];
    end

    for i = 1:numDiceBoundY
        point1 = [0,0+((i-1)*tileSize(2))];
        point2 = [dicingGridSize,0+((i-1)*tileSize(2))];
        point3 = [dicingGridSize,dicingBoundWidth+((i-1)*tileSize(2))];
        point4 = [0,dicingBoundWidth+((i-1)*tileSize(2))];
        polygons{end+1} = [point1;point2;point3;point4];
    end

    % Add dicing grid polygons to figure/CleWin
    for i = 1:length(polygons)
        data = polygons{i};
        x = data(:,1)';
        y = data(:,2)';
        
        x = x - (dicingBoundWidth/2) + xOff - (tileSize(1)*dicingTilingOffset(1));
        y = y - (dicingBoundWidth/2) + yOff - (tileSize(2)*dicingTilingOffset(2));

        if not(cleWin)
            hold on;
            plot(x, y, '.-');  
        else
            polygon([x;y]')
        end
    end
    
    if showWafer
        totXoff = -(dicingBoundWidth/2) + xOff - (tileSize(1)*dicingTilingOffset(1));
        totYoff = -(dicingBoundWidth/2) + yOff - (tileSize(2)*dicingTilingOffset(2));
        origin=[(dicingGridSize/2)+totXoff,(dicingGridSize/2)+totYoff];
        r=waferSize*25.4*1000/2;
        if cleWin
            circle(origin(1), origin(2), r);
        else
            ang=[0:0.01:2*pi]; 
            xp=r*cos(ang);
            yp=r*sin(ang);
            hold on;
            plot(origin(1)+xp,origin(2)+yp, '--');
        end
    end
end
                  
% ------------------------------------------------------------------
% ----------- Per tile ---------------------------------------------
% ------------------------------------------------------------------

% Loop per tile
for tileY = 1:numTileY
    for tileX = 1:numTileX 
        % Tile index
        tileIndex=((tileY-1)*numTileY)+tileX;
        
        tileOffsetX = (tileX-1)*tileSize(1);
        tileOffsetY = (tileY-1)*tileSize(2);
        
% ------------------------------------------------------------------
% ----------- Per tile - Parameter set -----------------------------
% ------------------------------------------------------------------
            % Resonator 
        resRatio    = pMap('Resonator conductor:gap ratio'); 
        if not(isequal(size(resRatio),[1,2])) resRatio=resRatio(tileIndex,:); end;
        resGap      = pMap('Resonator gap width');
        if not(isequal(size(resGap),[1,2])) resGap=resGap(tileIndex,:); end;
        resCond     = pMap('Resonator conductor width');
        if not(isequal(size(resCond),[1,2])) resCond=resCond(tileIndex,:); end;
        innerRadius = pMap('Resonator curve inner radius');
        if not(isequal(size(innerRadius),[1,2])) innerRadius=innerRadius(tileIndex,:); end;
        steps       = pMap('Resonator curve resolution');
        if not(isequal(size(steps),[1,2])) steps=steps(tileIndex,:); end;
            % ---
        couplerLen  = pMap('Resonator coupler length');
        if not(isequal(size(couplerLen),[1,2])) couplerLen=couplerLen(tileIndex,:); end;
        intLen      = pMap('Resonator interaction length');
        if not(isequal(size(intLen),[1,2])) intLen=intLen(tileIndex,:); end;
        intDist     = pMap('Resonator feedline-interaction distance');
        if not(isequal(size(intDist),[1,2])) intDist=intDist(tileIndex,:); end;
        totalLen    = pMap('Resonator total length'); % Only applicable if midCurves=true
        if not(isequal(size(totalLen),[1,2])) totalLen=totalLen(tileIndex,:); end;
            % ---
        midCurves   = pMap('Use extenable section');
        if not(isequal(size(midCurves),[1,2])) midCurves=midCurves(tileIndex,:); end;
        vertRatio   = pMap('Position of extendable section');
        if not(isequal(size(vertRatio),[1,2])) vertRatio=vertRatio(tileIndex,:); end;
          % Feedline
        feedRatio       = pMap('Feedline conductor:gap ratio:');
        if not(isequal(size(feedRatio),[1,1])) feedRatio=feedRatio(tileIndex); end;
        feedGap         = pMap('Feedline gap width');
        if not(isequal(size(feedGap),[1,1])) feedGap=feedGap(tileIndex); end;
        feedCond        = pMap('Feedline conductor width');
        if not(isequal(size(feedCond),[1,1])) feedCond=feedCond(tileIndex); end;
        feedDist        = pMap('Feedline-coupler distance');
        if not(isequal(size(feedDist),[1,2])) feedDist=feedDist(tileIndex,:); end;
        feedClearanceX  = pMap('Feedline-coupler X clearance');
        if not(isequal(size(feedClearanceX),[1,2])) feedClearanceX=feedClearanceX(tileIndex,:); end;
        feedLen         = pMap('Feedline length'); % Only applicable if doubleRes=false
        if not(isequal(size(feedLen),[1,1])) feedLen=feedLen(tileIndex); end;
          % End connector
        endConnect              = pMap('Use end connectors');
        if not(isequal(size(endConnect),[1,1])) endConnect=endConnect(tileIndex); end;
        endConnectRatioFac      = pMap('End connector ratio factor');
        if not(isequal(size(endConnectRatioFac),[1,1])) endConnectRatioFac=endConnectRatioFac(tileIndex); end;
        endConnectPadCondWidth  = pMap('End connector pad conductor width');
        if not(isequal(size(endConnectPadCondWidth),[1,1])) endConnectPadCondWidth=endConnectPadCondWidth(tileIndex); end;
        endConnectPadCondLen    = pMap('End connector pad conductor length');
        if not(isequal(size(endConnectPadCondLen),[1,1])) endConnectPadCondLen=endConnectPadCondLen(tileIndex); end;
        endConnectPadGapLen     = pMap('End connector pad gap length');
        if not(isequal(size(endConnectPadGapLen),[1,1])) endConnectPadGapLen=endConnectPadGapLen(tileIndex); end;
        endConnectExpanLen      = pMap('End connector expansion length');
        if not(isequal(size(endConnectExpanLen),[1,1])) endConnectExpanLen=endConnectExpanLen(tileIndex); end;
          % Multi resonator
        doubleRes   = pMap('Use two resonators');
        if not(isequal(size(doubleRes),[1,1])) doubleRes=doubleRes(tileIndex); end;
        resDist     = pMap('Resonator X separation');
        if not(isequal(size(resDist),[1,1])) resDist=resDist(tileIndex); end;
          % Labels
        useLabels       = pMap('Use labels');
        if not(isequal(size(useLabels),[1,1])) useLabels=useLabels(tileIndex); end;
        labels          = pMap('Labels');
        if not(isequal(size(labels),[1,1])) labels=labels(tileIndex); end;
        labelsOffset    = pMap('Labels offset ([x,y])');
        if not(isequal(size(labelsOffset),[1,2])) labelsOffset=labelsOffset(tileIndex,:); end;
        labelsSize      = pMap('Labels size');
        if not(isequal(size(labelsSize),[1,1])) labelsSize=labelsSize(tileIndex); end;

% ------------------------------------------------------------------
% ----------- Per tile - Add components to polygons store ----------
% ------------------------------------------------------------------ 

        % Store polygons. Empty for each new tile.
        polygons = {};
        
        % Add alignment marks
        if alignMarks
              % Polygon store indicies, m, k (k=1: first resonator, k=2: second resonator)
            m = 1;
            k = 3;
              % Add required marks to tile
            xOff=alignMarksOffset(1)-7000;
            yOff=alignMarksOffset(2)-tileSize(2)/2;
                % Always add top right
            polygons{k,m} = [alignMarkPoints(:,1)+tileSize(1)+xOff,alignMarkPoints(:,2)+tileSize(2)+yOff]; m=m+1;
            if tileX == 1 % Add top left
                polygons{k,m} = [alignMarkPoints(:,1)+xOff,alignMarkPoints(:,2)+tileSize(2)+yOff]; m=m+1;
            end
            if tileY == 1 % Add bottom right
                polygons{k,m} = [alignMarkPoints(:,1)+tileSize(1)+xOff,alignMarkPoints(:,2)+yOff]; m=m+1;
            end
            if tileX == 1 && tileY == 1 % Add bottom left
                polygons{k,m} = [alignMarkPoints(:,1)+xOff,alignMarkPoints(:,2)+yOff]; m=m+1;
            end
        end
        
        % Labels
        if cleWin && useLabels
            % Define a transformation matrix:
            xOff=-2500;
            yOff=500;
            x = tileOffsetX+labelsOffset(1)+xOff;
            y = tileOffsetY+labelsOffset(2)+yOff;
            m = [labelsSize 0 x; 0 labelsSize y; 0 0 1];
            % Add labels:
            text(labels{:},m);
        end

% ------------------------------------------------------------------
% ----------- Per tile - Per resonator -----------------------------
% ------------------------------------------------------------------  

        % Start Object generation
        if doubleRes numRes=2; else numRes=1; end;
        for k = 1:numRes
            
% ------------------------------------------------------------------
% ----------- Per tile - Per resonator - Prepare common objects ----
% ------------------------------------------------------------------           

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
            
            
            endConnectRatioFac      = pMap('End connector ratio factor');
            endConnectPadCondWidth  = pMap('End connector pad conductor width');
            endConnectPadCondLen    = pMap('End connector pad conductor length');
            endConnectPadGapLen     = pMap('End connector pad gap length');
            endConnectExpanLen      = pMap('End connector expansion length');
            
            % Generate end connector components
              % Useful params
            yMidOff = (feedGap+feedCond+feedGap)/2; 
            endConnectPadGapWidth = (endConnectPadCondWidth)/(feedRatio*endConnectRatioFac);
              % Points
            point1  = [0,yMidOff-(feedCond/2)];
            point2  = [endConnectExpanLen,(yMidOff-(endConnectPadCondWidth/2))];
            point3  = [endConnectExpanLen+endConnectPadCondLen,(yMidOff-(endConnectPadCondWidth/2))];
            point4  = [endConnectExpanLen+endConnectPadCondLen,(yMidOff+(endConnectPadCondWidth/2))];
            point5  = [endConnectExpanLen,(yMidOff+(endConnectPadCondWidth/2))];
            point6  = [0,yMidOff+(feedCond/2)];
            point7  = [0,yMidOff+(feedCond/2)+feedGap];
            point8  = [endConnectExpanLen,(yMidOff+(endConnectPadCondWidth/2)+endConnectPadGapWidth)];
            point9  = [endConnectExpanLen+endConnectPadCondLen+endConnectPadGapLen,(yMidOff+(endConnectPadCondWidth/2)+endConnectPadGapWidth)];
            point10 = [endConnectExpanLen+endConnectPadCondLen+endConnectPadGapLen,(yMidOff-(endConnectPadCondWidth/2)-endConnectPadGapWidth)];
            point11 = [endConnectExpanLen,(yMidOff-(endConnectPadCondWidth/2)-endConnectPadGapWidth)];
            point12 = [0,0];
              % Combine all
            endConnectorPoints = [point1;point2;point3;point4;point5;point6; ...
                                  point7;point8;point9;point10;point11;point12];
                      
% --------------------------------------------------------------------------
% ----------- Per tile - Per resonator - Add components to polygons store --
% --------------------------------------------------------------------------

            % Polygons index;
            m = 1;
            
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
            
            % End Connectors
            if endConnect
                xOff=tile+feedClearanceX(k);
                yOff=-(feedDist(k)+feedGap+feedCond+feedGap);
                polygons{k,m} = [endConnectorPoints(:,1)+xOff,endConnectorPoints(:,2)+yOff]; m=m+1;
                if not(doubleRes)
                    polygons{k,m} = [(-endConnectorPoints(:,1))+xOff-feedLen,endConnectorPoints(:,2)+yOff]; m=m+1;
                end
            end
           
        end
        
% ------------------------------------------------------------------
% ----------- Per tile - Copy polygons to corresponding canvas -----
% ------------------------------------------------------------------
        
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
                    
                    % Apply X,Y shifts for second resonator
                    x = x - resDist;
                    y = y - (sum(feedDist)+2*feedGap+feedCond);
                end
                
                % Apply X,Y shifts for specific tile
                x = x + tileOffsetX;
                y = y + tileOffsetY;

                if not(cleWin)
                    hold on;
                    plot(x, y, '.-');  
                else
                    polygon([x;y]')
                end
            end
        end
        
% ------------------------------------------------------------------
% ----------- Per tile - Display information -----------------------
% ------------------------------------------------------------------

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
        if not(cleWin) && verbose
            disp(message)
        end
        if cleWin && verbose
            error(message)
        end

    end
end
% ------------------------------------------------------------------
% ----------- Parameter set START ----------------------------------
% ------------------------------------------------------------------
pMap= containers.Map();
  % Resonator 
pMap('CPW conductor:gap ratio')               = 10.0/5.0;
pMap('CPW gap width')                         = 10.0;
pMap('CPW conductor width')                   = pMap('CPW gap width').*pMap('CPW conductor:gap ratio');
pMap('CPW length')                            = 1000;
pMap('CPW curve inner radius')                = 100;
pMap('CPW curve resolution')                  = 30;  
  % End connector
pMap('End connector ratio factor')            = 1;
pMap('End connector pad conductor width')     = 300;
pMap('End connector pad conductor length')    = 400;
pMap('End connector pad gap length')          = 150;
pMap('End connector expansion length')        = 610;
  % CleWin
pMap('Using CleWin')                          = false;

% ------------------------------------------------------------------
% ----------- Run once ---------------------------------------------
% ------------------------------------------------------------------

% Single set parameters
  % CleWin
cleWin  = pMap('Using CleWin');

% Start figure
if not(cleWin) figure('Position', [0, 0, 1800, 600]); end;

% ------------------------------------------------------------------
% ----------- Per tile ---------------------------------------------
% ------------------------------------------------------------------
% Resonator 
cpwRatio        = pMap('CPW conductor:gap ratio'); 
cpwGap          = pMap('CPW gap width');
cpwCond         = pMap('CPW conductor width');
cpwLen          = pMap('CPW length');
cpwInnerRadius  = pMap('CPW curve inner radius');
cpwCurveRes     = pMap('CPW curve resolution');

% End connector
endConnectRatioFac      = pMap('End connector ratio factor');
endConnectPadCondWidth  = pMap('End connector pad conductor width');
endConnectPadCondLen    = pMap('End connector pad conductor length');
endConnectPadGapLen     = pMap('End connector pad gap length');
endConnectExpanLen      = pMap('End connector expansion length');
        
% ------------------------------------------------------------------
% ----------- Per tile - Add components to polygons store ----------
% ------------------------------------------------------------------ 

% Store polygons. Empty for each new tile.
polygons = {};

% Derived parameters
% Resonator
innerRadiusInnerGap=cpwInnerRadius;
outerRadiusInnerGap=innerRadiusInnerGap+cpwGap;
innerRadiusOuterGap=outerRadiusInnerGap+cpwCond;
outerRadiusOuterGap=innerRadiusOuterGap+cpwGap;
tile = outerRadiusOuterGap;
tileShort = tile-innerRadiusInnerGap;
tileLong  = tile-tileShort;

% Generate cpw curved components
for i = 0:3
    start=i*0.5*pi;
    angles = linspace(0,0.25*2*pi,cpwCurveRes);
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


% Generate end connector components
% Useful params
yMidOff = (cpwGap+cpwCond+cpwGap)/2; 
endConnectPadGapWidth = (endConnectPadCondWidth)/(cpwRatio*endConnectRatioFac);
% Points
point1  = [0,yMidOff-(cpwCond/2)];
point2  = [endConnectExpanLen,(yMidOff-(endConnectPadCondWidth/2))];
point3  = [endConnectExpanLen+endConnectPadCondLen,(yMidOff-(endConnectPadCondWidth/2))];
point4  = [endConnectExpanLen+endConnectPadCondLen,(yMidOff+(endConnectPadCondWidth/2))];
point5  = [endConnectExpanLen,(yMidOff+(endConnectPadCondWidth/2))];
point6  = [0,yMidOff+(cpwCond/2)];
point7  = [0,yMidOff+(cpwCond/2)+cpwGap];
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

% CPW Curved sections
order = [1,2,3,4];
xOff = [0,200,400,600];
yOff = [0,0,0,0]+600;
for i = 1:length(order)
    indices = [2*(order(i)-1)+1:2*(order(i)-1)+2];
    polygons{end+1} = [quarterInner(indices(1),:,:)+xOff(i);quarterInner(indices(2),:,:)+yOff(i)]';
    polygons{end+1} = [quarterOuter(indices(1),:,:)+xOff(i);quarterOuter(indices(2),:,:)+yOff(i)]';
end

% CPW
xOff=0;
yOff=0;
horzLower = [[xOff+0,yOff+0]; ...
[xOff-cpwLen,yOff+0]; ...
[xOff-cpwLen,yOff+cpwGap]; ...
[xOff+0,yOff+cpwGap]]';
horzUpper = [horzLower(1,:);horzLower(2,:)+cpwGap+cpwCond];
polygons{end+1} = [horzLower(1,:);horzLower(2,:)]';
polygons{end+1} = [horzUpper(1,:);horzUpper(2,:)]';

% End Connectors  
polygons{end+1} = [endConnectorPoints(:,1),endConnectorPoints(:,2)];
        
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

        if not(cleWin)
            hold on;
            plot(x, y, '.-');  
        else
            polygon([x;y]')
        end
    end
end
        
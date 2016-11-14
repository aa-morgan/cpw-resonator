% Provided parameters
  % Resonator 
ratio = 9.065/5.0;
gap  = 5.0;
cond = gap * ratio;
innerRadius=40;
steps=54;

couplerLen=120;
intLen=5000;
intDist=1300;
totalLen=7624; % Fixed only if midCurves=true

midCurves=true;
vertRatio=0.3;

  % Feedline
feedRatio = 9.065/5.0;
feedGap = 10.0;                       
feedCond = feedGap * feedRatio;
feedDist = 20;
feedClearanceXright = 100;
feedLen = 500;
  % General
cleWin = false;
cleWinVerbose = false;

% Derived parameters
innerRadiusInnerGap=innerRadius;
outerRadiusInnerGap=innerRadiusInnerGap+gap;
innerRadiusOuterGap=outerRadiusInnerGap+cond;
outerRadiusOuterGap=innerRadiusOuterGap+gap;
tile = outerRadiusOuterGap;
tileShort = tile-innerRadiusInnerGap;
tileLong  = tile-tileShort;

if midCurves numCurved = 6; else numCurved = 2;end;

vert1Len=(vertRatio)*(intDist-(feedDist+(tile+tileLong)+(midCurves*2*(tile+tileLong))));
vert2Len=(1-vertRatio)*(intDist-(feedDist+(tile+tileLong)+(midCurves*2*(tile+tileLong))));
vertLen = [vert1Len,vert2Len];

curvedPathLen=0.5*pi*((outerRadiusOuterGap+innerRadiusInnerGap)/2);

if midCurves
    remainLen = totalLen - (numCurved*curvedPathLen+vert1Len+vert2Len+couplerLen+intLen);
    horzLen=remainLen/2;

    if remainLen <= 0
        error('Error: Sum of lengths are not long enough to satisfy the total length criteria!')
    end
else
    horzLen = 0;   
    totalLen = numCurved*curvedPathLen+vert1Len+vert2Len+couplerLen+intLen;
end

% Generate curved components
for i = 0:3
    start=i*0.5*pi;
    angles = linspace(0,0.25*2*pi,steps);
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
        quarterLLinnerGap = [rInner(1,:)+outerRadiusOuterGap;rInner(2,:)+outerRadiusOuterGap];;
        quarterLLouterGap = [rOuter(1,:)+outerRadiusOuterGap;rOuter(2,:)+outerRadiusOuterGap];;
    elseif start == 1.5*pi % UL
        quarterULinnerGap = [rInner(1,:)+outerRadiusOuterGap;rInner(2,:)];
        quarterULouterGap = [rOuter(1,:)+outerRadiusOuterGap;rOuter(2,:)];
    end

end
  % Combine
quarterInner = [quarterURinnerGap;quarterLRinnerGap;quarterLLinnerGap;quarterULinnerGap];
quarterOuter = [quarterURouterGap;quarterLRouterGap;quarterLLouterGap;quarterULouterGap];

% Start figure
if not(cleWin) figure('Position', [0, 0, 1800, 600]); end;

% Curved sections
order = [2,1,3,4,2,1];
xOff = [0,0,-horzLen-tile,-horzLen-tile,0,0];
yOff = [0,vert1Len+tile,vert1Len+tile+tileLong, ...
    vert1Len+2*tile+tileLong,vert1Len+2*(tile+tileLong), ...
    vert1Len+vert2Len+midCurves*2*(tile+tileLong)+tile];
for i = 1:length(order)
    if midCurves || i == 1 || i == 6
        indices = [2*(order(i)-1)+1:2*(order(i)-1)+2];
        if not(cleWin)
            hold on;
            plot(quarterInner(indices(1),:,:)+xOff(i),quarterInner(indices(2),:,:)+yOff(i), '.')
            hold on;
            plot(quarterOuter(indices(1),:,:)+xOff(i),quarterOuter(indices(2),:,:)+yOff(i), '.')  
        else
            polygon([quarterInner(indices(1),:,:)+xOff(i);quarterInner(indices(2),:,:)+yOff(i)]')
            polygon([quarterOuter(indices(1),:,:)+xOff(i);quarterOuter(indices(2),:,:)+yOff(i)]')
        end
    end
end

% Curve connectors
  % Vert's
for i = 1:2  
    xOff=tileLong;
    yOff=tile+((i-1)*(vertLen(1)+midCurves*2*(tile+tileLong)));
    vertLower = [[xOff+0,yOff+0]; ...
                  [xOff+0,yOff+vertLen(i)]; ...
                  [xOff+gap,yOff+vertLen(i)]; ...
                  [xOff+gap,yOff+0]]';
    vertUpper = [vertLower(1,:)+gap+cond;vertLower(2,:)];

    if not(cleWin)
        hold on;
        plot(vertLower(1,:),vertLower(2,:), 'x-')
        hold on;
        plot(vertUpper(1,:),vertUpper(2,:), 'x-');
    else
        polygon([vertLower(1,:);vertLower(2,:)]')
        polygon([vertUpper(1,:);vertUpper(2,:)]')
    end
end

  % Horz's
if midCurves
    for i = 0:1 
        xOff=0;
        yOff=tile+vert1Len+tileLong+(i*(tile+tileLong));
        horzLower = [[xOff+0,yOff+0]; ...
                      [xOff-horzLen,yOff+0]; ...
                      [xOff-horzLen,yOff+gap]; ...
                      [xOff+0,yOff+gap]]';
        horzUpper = [horzLower(1,:);horzLower(2,:)+gap+cond];

        if not(cleWin) 
            hold on;
            plot(horzLower(1,:),horzLower(2,:), 'x-')
            hold on;
            plot(horzUpper(1,:),horzUpper(2,:), 'x-')
        else
            polygon([horzLower(1,:);horzLower(2,:)]')
            polygon([horzUpper(1,:);horzUpper(2,:)]') 
        end
    end
end

% Coupler
xOff=0;
yOff=0;
couplerLower = [[xOff+0,yOff+0]; ...
                [xOff-couplerLen,yOff+0]; ...
                [xOff-couplerLen,yOff+gap]; ...
                [xOff+0,yOff+gap]]';
couplerUpper = [couplerLower(1,:);couplerLower(2,:)+gap+cond];
xOff=-couplerLen;
yOff=0;
couplerEnd = [[xOff+0,yOff+0]; ...
              [xOff-gap,yOff+0]; ...
              [xOff-gap,yOff+gap+cond+gap]; ...
              [xOff+0,yOff+gap+cond+gap]]';       
if not(cleWin) 
    hold on;
    plot(couplerLower(1,:),couplerLower(2,:), 'x-')
    hold on;
    plot(couplerUpper(1,:),couplerUpper(2,:), 'x-')
    hold on;
    plot(couplerEnd(1,:),couplerEnd(2,:), 'x-')
else
    polygon([couplerLower(1,:);couplerLower(2,:)]')
    polygon([couplerUpper(1,:);couplerUpper(2,:)]')
    polygon([couplerEnd(1,:);couplerEnd(2,:)]')
end

% Interaction region
xOff=0;
yOff=vert1Len+vert2Len+(tileLong+tile)+midCurves*2*(tileLong+tile);
horz2Lower = [[xOff+0,yOff+0]; ...
              [xOff-intLen,yOff+0]; ...
              [xOff-intLen,yOff+gap]; ...
              [xOff+0,yOff+gap]]';
horz2Upper = [horz2Lower(1,:);horz2Lower(2,:)+gap+cond];

if not(cleWin) 
    hold on;
    plot(horz2Lower(1,:),horz2Lower(2,:), 'x-')
    hold on;
    plot(horz2Upper(1,:),horz2Upper(2,:), 'x-');
else
    polygon([horz2Lower(1,:);horz2Lower(2,:)]')
    polygon([horz2Upper(1,:);horz2Upper(2,:)]')
end

% Feedline
xOff=tile+feedClearanceXright;
yOff=-(feedDist+feedGap+feedCond+feedGap);
feedLower = [[xOff+0,yOff+0]; ...
             [xOff-feedLen,yOff+0]; ...
             [xOff-feedLen,yOff+feedGap]; ...
             [xOff+0,yOff+feedGap]]';
feedUpper = [feedLower(1,:);feedLower(2,:)+feedGap+feedCond];

if not(cleWin) 
    hold on;
    plot(feedLower(1,:),feedLower(2,:), 'x-')
    hold on; 
    plot(feedUpper(1,:),feedUpper(2,:), 'x-')
else
    polygon([feedLower(1,:);feedLower(2,:)]')
    polygon([feedUpper(1,:);feedUpper(2,:)]')
end

% Figure settings
if not(cleWin) 
    xlim([-intLen-20, feedClearanceXright+tile+20]);  
    ylim([-feedDist-feedGap-feedCond-feedGap-20, intDist-feedDist+gap+cond+gap+20])
end

% Display Info
if midCurves midStatus = 'On'; else midStatus = 'Off';end;
message = sprintf(strcat(...
        'Info, (units: micrometers)', '\n', ...
        '\tTotal length:\t\t', num2str(totalLen), '\n', ...
        '\tCoupler length:\t\t', num2str(couplerLen), '\n', ...
        '\tInteraction length:\t\t', num2str(intLen), '\n', ...
        '\tFeedline-coupler distance:\t', num2str(feedDist), '\n', ...
        '\tFeedline-interaction distance:\t', num2str(intDist), '\n', ...
        '\tLength extension section:\t', midStatus, '\n' ...
        ));
if not(cleWin)
    disp(message)
end
if cleWin && cleWinVerbose
    error(message)
end
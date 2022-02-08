function [vxO,vyO,vxClip,vyClip] = VoronoiLimitRectSquare(X, Y, VoronoiBoundary, ...
                                                clip2boundary, showFigure)

%%***********************************************************************%
%*                        Bounded Voronoi Diagram                       *%
%*                Clips the extended edges of a Voronoi Diagram         *%
%*                to the predefined rectangular or square boundary.     *%
%*                                                                      *%
%* Author: Preetham Manjunatha                                          *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 02/08/2022                                                     *%
%************************************************************************%
%
%************************************************************************%
%
% Usage: [vxO,vyO,vxClip,vyClip] = VoronoiLimitRectSquare(X, Y, VoronoiBoundary)
%        [___]                   = VoronoiLimitRectSquare(___, clip2boundary, showFigure)  
% Inputs:
%
%           X                   - X coordinates of the points
%           Y                   - Y coordinates of the points
%           VoronoiBoundary     - Boundary [xmin, xmax, ymin, ymax]
%           clip2boundary       - Clip the edge points to boundaries limit
%                                 (optional)
%           showFigure          - Show the bounded Voronoi figure and
%                                 original Voronoi figures (optional)
% 
% Outputs: 
%
%           vxO                 - Original Voronoi edges X-coordinates
%           vyO                 - Original Voronoi edges Y-coordinates
%           vxClip              - Clipped Voronoi edges Y-coordinates
%           vyClip              - Clipped Voronoi edges Y-coordinates
%--------------------------------------------------------------------------
% Example 1: Rectangular Boundary
% imX1 = 0; imX2 = 512;
% imY1 = 0; imY2 = 256;
% VoronoiBoundary = [imX1 imX2 imY1 imY2]; % [xmin, xmax, ymin, ymax]
% clip2boundary = 1;
% numPoints = 25;
% X = randi([imX1,imX2], 1, numPoints);
% Y = randi([imY1,imY2], 1, numPoints);
% [vxO,vyO,vxClip,vyClip] = VoronoiLimitRectSquare(X, Y, VoronoiBoundary, clip2boundary, 0);
%
% Example 2: Square Boundary
% imX1 = 0; imX2 = 512;
% imY1 = 0; imY2 = 512;
% VoronoiBoundary = [imX1 imX2 imY1 imY2]; % [xmin, xmax, ymin, ymax]
% clip2boundary = 1;
% numPoints = 25;
% X = randi([imX1,imX2], 1, numPoints);
% Y = randi([imY1,imY2], 1, numPoints);
% [vxO,vyO,vxClip,vyClip] = VoronoiLimitRectSquare(X, Y, VoronoiBoundary, clip2boundary, 0);
%
% Computational time:
% Intel i7-3770 @ 3.4 GHz, 32 GB RAM, Windows 11 
% Random 200 points, boundary [0 512 0 512], 100 iterations
% MATLAB Mapping Toolbox line intersection  : 0.0206 seconds
% Custom parametric line intersection       : 0.0020 seconds
% Overall improvement                       : 9.7100%


%------------------------------------------------------------------------------------------------------------------------
% nargin check
if nargin < 3
    error('Not enough input arguments.');
elseif nargin > 5
    error('Too many input arguments.');
end

if nargin == 3
    %-----------------------
    % Show intersection plot
    showIntersectionPlot = 0;

    % Clip the boundary vertices
    clip2boundary = 1;

    % Show bounded Voronoi diagram 
    showFigure = 1;
end

if nargin == 4
    %-----------------------
    % Show intersection plot
    showIntersectionPlot = 0;

    % Show bounded Voronoi diagram 
    showFigure = 1;
end

if nargin == 5
    %-----------------------
    % Show intersection plot
    showIntersectionPlot = 0;
end

%------------------------------------------------------------------------------------------------------------------------
% Get MAP toolbox license
mapToolBoxLicense = 0;

%------------------------------------------------------------------------------------------------------------------------
% Generate voronoi points
[vxO,vyO] = voronoi(X,Y);
vx = vxO;
vy = vyO;

%------------------------------------------------------------------------------------------------------------------------
% Voronoi boundaries
imX1 = VoronoiBoundary(1);
imX2 = VoronoiBoundary(2);
imY1 = VoronoiBoundary(3);
imY2 = VoronoiBoundary(4);

%% ------------------------------------------------------------------------------------------------------------------------
% Discard out of boundary points
[rowvX,colvX] = find(vx < imX1 | vx > imX2);
[rowvY,colvY] = find(vy < imY1 | vy > imY2);

% Indices to unique values in column 3
[~, indvY] = unique(colvY);
[~, indvX] = unique(colvX);

% Duplicate values
duplicatevX = colvX(setdiff(1:numel(colvX), indvX));
duplicatevY = colvY(setdiff(1:numel(colvY), indvY));

% Valid vertices within the limits
vx(:,[duplicatevX;duplicatevY]) = [];
vy(:,[duplicatevX;duplicatevY]) = [];

%% ------------------------------------------------------------------------------------------------------------------------
% Left edge of a rectangle/square
[row_left,col_left] = find(vx < imX1);
leftEdgeVertsX = vx(:, col_left);
leftEdgeVertsY = vy(:, col_left);

for i = 1:size(leftEdgeVertsX,2)
    if mapToolBoxLicense
        [xiLeftval, yiLeftval] = polyxpoly(leftEdgeVertsX(:,i), leftEdgeVertsY(:,i), [imX1 imX1], [imY1 imY2]);
    else
        [xiLeftval, yiLeftval] = linexline(leftEdgeVertsX(:,i), leftEdgeVertsY(:,i), [imX1 imX1], [imY1 imY2], showIntersectionPlot);
    end

    if isempty(xiLeftval) || isempty(yiLeftval)
        xiLeft(i) = NaN; yiLeft(i) = NaN; 
    else
        xiLeft(i) = xiLeftval;  yiLeft(i) = yiLeftval;
    end
end
[rowNaNLeft, colNaNLeft] = find(isnan(xiLeft));
row_left(colNaNLeft) = [];
col_left(colNaNLeft) = [];
xiLeft(isnan(xiLeft)) = [];
yiLeft(isnan(yiLeft)) = [];

%------------------------------------------------------------------------------------------------------------------------
% Right edge of a rectangle/square
[row_right,col_right] = find(vx > imX2);
rightEdgeVertsX = vx(:, col_right);
rightEdgeVertsY = vy(:, col_right);

for i = 1:size(rightEdgeVertsX,2)  
    if mapToolBoxLicense
        [xiRightval, yiRightval] = polyxpoly(rightEdgeVertsX(:,i), rightEdgeVertsY(:,i), [imX2 imX2], [imY1 imY2]);
    else
        [xiRightval, yiRightval] = linexline(rightEdgeVertsX(:,i), rightEdgeVertsY(:,i), [imX2 imX2], [imY1 imY2], showIntersectionPlot);
    end

    if isempty(xiRightval) || isempty(yiRightval)
        xiRight(i) = NaN; yiRight(i) = NaN; 
    else
        xiRight(i) = xiRightval;  yiRight(i) = yiRightval;
    end
end
[rowNaNRight, colNaNRight] = find(isnan(xiRight));
row_right(colNaNRight) = [];
col_right(colNaNRight) = [];
xiRight(isnan(xiRight)) = [];
yiRight(isnan(yiRight)) = [];

%------------------------------------------------------------------------------------------------------------------------
% Bottom edge of a rectangle/square
[row_bot,col_bot] = find(vy < imY1);
botEdgeVertsX = vx(:, col_bot);
botEdgeVertsY = vy(:, col_bot);

for i = 1:size(botEdgeVertsX,2)
    if mapToolBoxLicense
        [xiBotval, yiBotval] = polyxpoly(botEdgeVertsX(:,i), botEdgeVertsY(:,i), [imX1 imX2], [imY1 imY1]);
    else
        [xiBotval, yiBotval] = linexline(botEdgeVertsX(:,i), botEdgeVertsY(:,i), [imX1 imX2], [imY1 imY1], showIntersectionPlot);
    end
    if isempty(xiBotval) || isempty(yiBotval)
        xiBot(i) = NaN; yiBot(i) = NaN; 
    else
        xiBot(i) = xiBotval;  yiBot(i) = yiBotval;
    end
end
[rowNaNBot, colNaNBot] = find(isnan(xiBot));
row_bot(colNaNBot) = [];
col_bot(colNaNBot) = [];
xiBot(isnan(xiBot)) = [];
yiBot(isnan(yiBot)) = [];

%------------------------------------------------------------------------------------------------------------------------
% Top edge of a rectangle/square
[row_top,col_top] = find(vy > imY2);
topEdgeVertsX = vx(:, col_top);
topEdgeVertsY = vy(:, col_top);

for i = 1:size(topEdgeVertsX,2)
    if mapToolBoxLicense
        [xiTopval, yiTopval] = polyxpoly(topEdgeVertsX(:,i), topEdgeVertsY(:,i), [imX1 imX2], [imY2 imY2]);
    else
        [xiTopval, yiTopval] = linexline(topEdgeVertsX(:,i), topEdgeVertsY(:,i), [imX1 imX2], [imY2 imY2], showIntersectionPlot);
    end
    if isempty(xiTopval) || isempty(yiTopval)
        xiTop(i) = NaN; yiTop(i) = NaN; 
    else
        xiTop(i) = xiTopval;  yiTop(i) = yiTopval;
    end
end
[rowNaNTop, colNaNTop] = find(isnan(xiTop));
row_top(colNaNTop)  = [];
col_top(colNaNTop)  = [];
xiTop(isnan(xiTop)) = [];
yiTop(isnan(yiTop)) = [];

%% ------------------------------------------------------------------------------------------------------------------------
% Stack the arrays to cell
row2clip = {row_left,row_right,row_top,row_bot};
col2clip = {col_left,col_right,col_top,col_bot};
xyIntersections = {[xiLeft; yiLeft]', [xiRight; yiRight]', [xiTop; yiTop]', [xiBot; yiBot]'};

%% ------------------------------------------------------------------------------------------------------------------------
% Replace the values with intersection/updated points
for i = 1:length(row2clip)
    for j = 1:length(row2clip{i})
        vx(row2clip{i}(j), col2clip{i}(j)) = xyIntersections{i}(j,1);
        vy(row2clip{i}(j), col2clip{i}(j)) = xyIntersections{i}(j,2);
    end
end

%% ------------------------------------------------------------------------------------------------------------------------
% Clip vertices values
if clip2boundary
    vxClip = min(max(round(vx),1),imX2);
    vyClip = min(max(round(vy),1),imY2);
else
    vxClip = vx;
    vyClip = vy;
end

%% ------------------------------------------------------------------------------------------------------------------------
% Display figure
if showFigure
    fh = figure();
    fh.WindowState = 'maximized';

    ax1 = subplot(3,1,1);
    h1 = voronoi(X,Y);
    axis equal
    title('Original Voronoi Diagram')
    
    % Dummy image
    I = true(imY2,imX2);
    ax2 = subplot(3,1,2);
    imshow(I)
    hold on
    h2 = voronoi(X,imY2-Y);
    hold off
    axis equal
    title('Original Voronoi Diagram Overlayed on Image')
    
    ax3 = subplot(3,1,3);
    hold on
    rectangle('Position',[imX1 imY1 imX2-imX1 imY2-imY1]')
    for i = 1:size(vxClip,2)
        plot(vxClip(:,i), vyClip(:,i),'Color',rand(1,3), 'LineWidth',3)
    end
    hold off
    axis equal
    title('Bounded Voronoi Diagram')
    linkaxes([ax1 ax2 ax3],'xy')

    exportgraphics(gcf,'clippedVoronoiPlot.png')
end
end

%--------------------------------------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------------------------------------
function [xi,yi] = linexline(L1x, L1y, L2x, L2y, showIntersectionPlot)
%%***********************************************************************%
%*                    Line to line interection point                    *%
%*            Finds the interection of the two line segments.           *%
%*                                                                      *%
%*                                                                      *%
%* Author: Preetham Manjunatha                                          *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 02/08/2022                                                     *%
%************************************************************************%
%
%************************************************************************%
%
% Usage: [xi,yi] = linexline(L1x, L1y, L2x, L2y)
%
% Inputs:
%
%           L1x                     - Line 1 x1 and x2 coordinates [x1, x2]
%           L1y                     - Line 1 y1 and y2 coordinates [y1, y2]
%           L2x                     - Line 2 x1 and x2 coordinates [x3, x4]
%           L2y                     - Line 2 y1 and y2 coordinates [y3, y4]
%           showIntersectionPlot    - show intersection plot (0 or 1)
% 
% Outputs: 
%
%           xi          - Interection point, x coordinate (NaN if no
%                         interesction)
%           yi          - Interection point, y coordinate (NaN if no
%                         interesction)
%--------------------------------------------------------------------------

% Data
x1 = L1x(1);
y1 = L1y(1);
x2 = L1x(2);
y2 = L1y(2);
x3 = L2x(1);
y3 = L2y(1);
x4 = L2x(2);
y4 = L2y(2);

% Line segments intersect parameters
u = ((x1-x3)*(y1-y2) - (y1-y3)*(x1-x2)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));

% Check if intersection exists, if so then store the value
if (u >= 0 && u <= 1.0) && (t >= 0 && t <= 1.0)
    xi = ((x3 + u * (x4-x3)) + (x1 + t * (x2-x1))) / 2; 
    yi = ((y3 + u * (y4-y3)) + (y1 + t * (y2-y1))) / 2;
else
    xi = NaN;
    yi = NaN;
end

if showIntersectionPlot
    % Plot the lines
    plot([x1 x2], [y1 y2], 'LineWidth', 3)
    hold on
    plot([x3 x4], [y3 y4], 'LineWidth', 3)
    
    % Plot intersection points
    plot(x3 + u * (x4-x3), y3 + u * (y4-y3), 'ro', 'MarkerSize', 15)
    plot(x1 + t * (x2-x1), y1 + t * (y2-y1), 'bo', 'MarkerSize', 15)
    hold off
    xlabel('X'); ylabel('Y') 
    grid on

    ax = gca;
    exportgraphics(ax,'LineXPlot.png')
end
end
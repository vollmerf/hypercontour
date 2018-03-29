% File    : test.m
% System  : Matlab/Octave
% Purpose : test program for hypercontour.m
% Author  : Frederick W. Vollmer
% Date    : June 11, 2017 

% get comma delimited test file
[filename, pathname] = uigetfile( {'*.csv'});
rphi = csvread([pathname, filename]);

% create equidistant plot (Elliott plot)
[points, lines, ticks, frame, grid] = hypercontour(rphi, 'deg');

% create high quality stereographic plot
%[points, lines, ticks, frame, grid] = hypercontour(rphi, 'deg,stg', 6, 40, 5, 50);

% create equidistant (log) RfPhi plot
%[points, lines, ticks, frame, grid] = hypercontour(rphi, 'deg,rfp');

% set up figure
figure;
hold on;
axis([-1.0 1.0 -1.0 1.0]);
axis('equal');
axis('off');

% define WBGYR dolormap, image outside plot is white
T = [255,255,255;0,0,255;0,255,0;255,255,0;255,0,0]./255; 
x = [0,64,128,192,255];
cmap = interp1(x/255,T,linspace(0,1,255));

% or use grayscale, cmap = (1-gray), or:
T = [255,255,255;32,32,32]./255; 
x = [0,255];
gmap = interp1(x/255,T,linspace(0,1,255));

% select color map
%colormap(gmap);
colormap(cmap);

% plot grid as gradient map
imagesc(-1:1, -1:1, grid);

% plot ticks, array of (x1, y1, x2, y2)
n = length(ticks);
for i = 1:n
  lx = [ticks(i,1), ticks(i,3)];
  ly = [ticks(i,2), ticks(i,4)]; 
  line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 1);
end

% plot contours, array of (x1, y1, x2, y2)
[n,m] = size(lines);
for i = 1:n
  lx = [lines(i,1), lines(i,3)];
  ly = [lines(i,2), lines(i,4)]; 
  line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 2);
end

% plot frame, array of (x1, y1, x2, y2)
[n,m] = size(frame);
for i = 1:n
  lx = [frame(i,1), frame(i,3)];
  ly = [frame(i,2), frame(i,4)]; 
  line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 2);
end

% plot points, array of (x,y)
px = points(:,1); 
py = points(:,2); 
h = plot(px, py, 'o');
set(h(1),'MarkerEdgeColor','w','MarkerFaceColor','k')

hold off;

% print as eps
%print -deps fig_bw.eps;  
%print -depsc fig_color.eps; 
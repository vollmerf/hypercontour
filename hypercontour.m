function [points, lines, ticks, frame, grid] = hypercontour(rphi, options, rmax, kappa, nlevels, nnodes)
% Plots and contours geological fabric and strain data on hyperbaloidal projections. 
%
% INPUT
% -----
% rphi     : array (R, phi), R = strain ratio (A/B), phi = orientation of long 
%            axis (A) from x.   
% options  : comma separated string with non-default options (default = ''):
%            angle format: 
%              ''       = radians (default)
%              'rad'    = radians
%              'deg'    = degrees
%              'grd'    = gradians
%            projection (polar equidistant is default):
%              ''       = equidistant (log R, Elliott plot)
%              'eqd'    = equidistant 
%              'eqa'    = equal-area 
%              'stg'    = stereographic
%              'ort'    = orthographic
%              'gno'    = gnomic
%              'lin'    = exponential (linear R)
%              'rdl'    = radial
%              'rfp'    = Rf/phi (cylindrical instead of polar)
%            contouring method:
%              ''       = contour (default)
%              'ctr'    = contour
%              'nnm'    = no normalization
%              'nct'    = no contouring (only points will be returned)
%            grid to image interpolation: 
%              ''       = 5 parts (default)
%              'gi0'    = off
%              'gi2'    = 2 parts
%              'gi3'    = 3 parts
%              'gi4'    = 4 parts
%              'gi5'    = 5 parts
%              'gi6'    = 6 parts
%              'gi8'    = 8 parts
%              'gia'    = 10 parts
%            frame and ticks:     
%	             ''       = draw circle and tics  
%	             'ntc'    = draw circle, without tics              
%	             'nfr'    = no frame 
%            grid:
%              ''       = grid     
%              'ngd'    = no grid     
% rmax     : maximum R value on plot, default = 0 (automatic).    
% kappa    : weighting parameter, default = 40. 
% nlevels  : number of levels spaced over the probability density distribution 
%            (pdd), 5 will divide the pdd into 5, giving 4 contour lines 
%            spaced at 20% of the distribution, default = 5. 
% nnodes   : number of grid nodes, higher is more accurate but slower, 30 is 
%            good, but 50 is recommended for final plots, default = 30.
%
% OUTPUT
% ------          
% points : projected data points in unit circle or square as array of 
%          [x,y] = [points(:,1), points(:,2)]
% lines  : projected contour line segments in unit circle or square as array of
%          [x1,y1,x2,y2] = [lines(:,1), lines(:,2), lines(:,3), lines(:,4)].
% ticks  : tick marks as line segments, [x1,y1,x2,y2].
% frame  : bounding circle or square as line segments.
% grid   : grid for display of color gradient: imagesc(-1:1, -1:1, grid).   
% 
% USAGE
% -----
% [points] = hypercontour(rphi, 'deg');
% [points, lines, ticks, frame] = hypercontour(rphi);
% [points, lines, ticks, frame, grid] = hypercontour(rphi);
% [points, lines, ticks, frame] = hypercontour(rphi, 'deg,rfp', 5, 60, 10, 50);
%
% All input parameters except 'data' are optional. The algorithm and function 
% are described in:
%
%   Vollmer, F.W., 2018. Automatic contouring of geological fabric and finite 
%   strain data on the unit hyperboloid. Computers & Geosciences, 
%   https://doi.org/10.1016/j.cageo.2018.03.006
%
% This paper should be referenced in publications or presentations using this 
% or derivative code. See that paper and the files README.md, LICENSE.md, 
% CITATION.md license and additional information. 
%

% END HELP
% -----------------------------------------------------------------------------
%
% File    : hypercontour.m
% Version : 1.0.0.9
% System  : Matlab/Octave
% Author  : Frederick W. Vollmer
% Date    : 29 Mar 2018
% Notice  : Copyright (c) 2017-2018, Frederick W. Vollmer
%
% DESCRIPTION
% -----------
% MATLAB/Octave function for plotting and contouring hyperbaloidal projections 
% of geological fabric and finite strain data density calculations done on the 
% unit hyperbaloid (Vollmer, 2018). Options are given for equidistant 
% (Elliott), equal-area, stereographic, orthographic, exponential, and radial 
% projections, as polar azimuthal or cylindrical (cartesian, RfPhi-type) plots. 
%
% The data must be in a comma delimited csv text file with one (R, phi) pair 
% per line, where R = strain ratio (max/min), phi = orientation of long (max) 
% axis from x. Contours are equally spaced over the probability density 
% distribution. Options are specified with an input string, such as 'deg,ort', 
% see above help for all options.
% 
% Publications or presentations using this or derivative code to produce 
% figures or other content should cite the following paper: 
%  
%   Vollmer, F.W., 2018. Automatic contouring of geological fabric and finite 
%   strain data on the unit hyperboloid. Computers & Geosciences, 
%   https://doi.org/10.1016/j.cageo.2018.03.006
%
% See this paper and the files README.md, LICENSE.md, and CITATION.md for  
% additional information. 
%
% Fabric and finite strain data contouring is also implemented in the 
% standalone program EllipseFit by this author, which is free, has numerous 
% additional functions, and is faster. It runs on Macintosh, Windows, and 
% Linux platforms, and is recommended over this function for non-MATLAB/Octave 
% use. It can be downloaded for free from: 
%
%   www.frederickvollmer.com/ellipsefit
%   www.newpaltz.edu/~vollmerf
%
% Please contact the author for any bug reports or feature requests:
%
% Frederick W. Vollmer
% vollmerf@newpaltz.edu 
% vollmerf@gmail.com
%
%------------------------------------------------------------------------------

%function [points, lines, ticks, frame, grid] = hypercontour(rphi, options, rmax, kappa, nlevels, nnodes)
  global opts;
  switch nargin
    case 1
      options = '';
      rmax = 0.0; % automatic
      kappa = 40.0;
      nlevels = 5;
      nnodes = 30;
    case 2
      rmax = 0.0;
      kappa = 40.0;
      nlevels = 5;
      nnodes = 30;
    case 3
      kappa = 40.0;
      nlevels = 5;
      nnodes = 30;
    case 4
      nlevels = 10;
      nnodes = 30;
    case 5
      nnodes = 30;
    case 6
      nnodes = nnodes;
    otherwise
      return % error
  end  
  if rmax < 1.0
    rmax = ceil(max(rphi(:,1))) + 1.0;
  end
  opts = getOptions(options);
  if nargout < 5 % no grid
    opts.grid = 0
  end
  if nargout < 4 % no frame
    opts.frame = 0
  end
  if nargout < 3 % no ticks
    opts.ticks = 0
  end
  if nargout < 2 % no lines
    opts.contour = 0
  end
  if nargout < 1 
    return % error
  end
  if (opts.angfmt == 1) % degrees
    f = pi/180.0;
  elseif (opts.angfmt == 2) % gradians
    f = pi/200.0;
  else % radians
    f = 1.0;
  end
  for i = 1:length(rphi)
    rphi(i,2) = rphi(i,2) * f; % to radians
    [points(i,1), points(i,2)] = rPhiToXYUnit(rphi(i,1), rphi(i,2), rmax);  
  end
  if opts.contour
    [grid, lines] = contour(rphi, rmax, kappa, nlevels, nnodes);
    if opts.grid
      grid = processGrid(grid, rmax);
    end
  end
  ticks = drawTicks(rmax);
  frame = drawFrame();
end

function [grid] = processGrid(grid, rmax)
  global opts;
  grid = grid';
  [n, m] = size(grid);
  [x y] = meshgrid(1:n);
  if (opts.interp < 0.0)
    zi = grid;
  else
    [xi yi] = meshgrid(1:opts.interp:n);
    zi = interp2(x,y,grid,xi,yi);
  end
  [ni, mi] = size(zi);
  [xi yi] = meshgrid(1:ni);
  if ~opts.rfp
    r = 0.5 * (ni-1);
    r2 = r * r;
    % clip to circle, but NaN is implementation depependent
    % so use 0, and colormap starting with white
    %zi((xi - r - 1).^2 + (yi - r - 1).^2 > r2) = NaN;
    zi((xi - r - 1).^2 + (yi - r - 1).^2 > r2) = 0.0;
  end
  grid = zi;
end

% rToZeta - projects R to zeta. Ref: Yamaji, 2008. 
function [z] = rToZeta(r)
  global opts;
  switch opts.proj
    case 0 % equidistant (Elliott)
      z = log(r);  
    case 1 % equal-area 
      t = sqrt(r);
      z = t - 1.0/t;
    case 2 % stereographic
      t = sqrt(r);
      s = 1.0 / t;
      z = 2.0 * (t - s) / (t + s);
    case 3 % orthographic
      z = 0.5 * (r - 1.0/r);
    case 4 % gnomic
      t = r * r;
      z = (t-1)/(t+1);
    case 5 % linear (exponential)
      z = r - 1;
    case 6 % cylindrical (radial) 
      z = 0.5 * (r + 1.0/r) - 1;
  end
end  

% zetaToR - back projects zeta to R.
function [r] = zetaToR(z)
  global opts;
  switch opts.proj
    case 0 % equidistant (Elliott)
      r = exp(z); 
    case 1 % equal-area 
      t = z + sqrt(z*z + 4.0);
      r = t * t * 0.25;
    case 2 % stereographic
      t = 0.5 * z;
      r = (1.0 + t)/(1.0 - t);
    case 3 % orthographic
      r = z + sqrt(z * z + 1);
    case 4 % gnomic
      t = 0.0;
      if z < 0.99 % 0..1, 1 -> inf, 0.99 -> 199
        t = sqrt((1.0+z)/(1.0-z));
      end
      if t < 50.001 % cap r
        r = t; 
      else 
        r = 0.0;
      end
    case 5 % linear (exponential)
      r = z + 1.0;
    case 6 % cylindrical (radial) 
      t = z + 1.0;
      r = t + sqrt(t * t - 1.0);
  end
end    

% rPhiToXY - projects R, phi to cartesian coordinates of unit hyperbaloidal
% projection. Maps to [-1..-1, +1..+1] to overlie unit image.                                              
function [x, y] = rPhiToXYUnit(r, phi, rmax)
  global opts; 
  z = rToZeta(r);
  zm = rToZeta(rmax);
  s = z / zm;
  if opts.rfp
    p = phi;
    if p < -0.5 * pi
      p = p + pi;
    elseif p > 0.5 * pi
      p = p - pi;
    end       
    x = 2.0 * p/pi;
    y = 2.0 * s - 1.0;
  else % polar
    x = s * cos(2.0 * phi);
    y = s * sin(2.0 * phi);
  end
end

% xYToRPhi - back projects cartesian coordinates of hyperbaloidal projection.
% Not scaled from unit plot.
function [r, phi] = xYToRPhi(x, y, rmax)
  global opts; 
  zm = rToZeta(rmax);
  if opts.rfp 
    z = (y + zm) * 0.5; 
    r = zetaToR(z);     
    phi = x * (0.5 * pi / zm);    
    if phi < -0.5 * pi
      phi = phi + pi;
    elseif phi > 0.5 * pi
      phi = phi - pi;
    end       
  else % polar
    t = sqrt(x*x + y*y);
    r = zetaToR(t); 
    phi = 0.5 * atan2(y, x); 
  end
end

% rhoPsiToH - set as a hyperbolic position vector from rho and psi. For strain 
% ellipses: rho = ln(R), psi = 2 phi. Ref: Yamaji, 2008. }
function [h] = rhoPsiToH(rho, psi)
  s = sinh(rho);
  h(1) = cosh(rho);
  h(2) = s * cos(psi);
  h(3) = s * sin(psi);
end

% rPhiToH - converts R, phi to hyperbaloidal point. 
function [h] = rPhiToH(r, phi) 
  if r < 1.0 % error
    rho = 0.0;
  else
    rho = log(r);
  end
  psi = 2.0 * phi;
  h = rhoPsiToH(rho, psi);
end

% dotH - hyperbolic inner product. Ref: Yamaji, 2008, eqn 4. 
function [d] = dotH(a, b)
  d = -a(1) * b(1) + a(2) * b(2) + a(3) * b(3);
end

% lineCircleInt - determine intersection parameters for line segment and 
% circle. Adopted from Rankin 1989, p.220.                               
function [t1, t2, visible] = lineCircleInt(x1, y1, x2, y2, xc, yc, r) 
  visible = 0; % FALSE
  t1 = 0.0;
  t2 = 1.0;
  dx = x2-x1; 
  dy = y2-y1; 
  dxc = x1-xc; 
  dyc = y1-yc;
  a = dx*dxc + dy*dyc; 
  b = dx*dx + dy*dy; 
  c = dxc*dxc + dyc*dyc - r*r;
  disc = a*a - b*c;
  if ((disc > 0.0) && (abs(b) > 1e-9)) 
    d = sqrt(disc);
    t1 = (-a + d)/b; 
    t2 = (-a - d)/b;
    if (t1 > t2) 
      t = t1; 
      t1 = t2; 
      t2 = t; 
    end
    visible = 1; % TRUE
  end
end

% clipLineCircle - clip line segment to circle. 
function [cx1, cy1, cx2, cy2, visible] = clipLineCircle(xc, yc, r, x1, y1, x2, y2) 
  cx1 = x1;
  cy1 = y1;
  cx2 = x2;
  cy2 = y2;
  visible = 0; % FALSE 
  if (((x1 < xc-r) && (x2 < xc-r)) || ((x1 > xc+r) && (x2 > xc+r)))
    return;                  
  end
  if (((y1 < yc-r) && (y2 < yc-r)) || ((y1 > yc+r) && (y2 > yc+r)))
    return;                  
  end
  [t1, t2, vis] = lineCircleInt(x1,y1,x2,y2,xc,yc,r);  
  if (vis == 0)
    return;  
  end
  if ((t2 < 0.0) || (t1 > 1.0))
    visible = 0; % FALSE 
    return;
  end
  if (t1 > 0.0) 
    cx1 = x1 + (x2-x1) * t1; 
    cy1 = y1 + (y2-y1) * t1; 
  end
  if (t2 < 1.0)  
    cx2 = x1 + (x2-x1) * t2; 
    cy2 = y1 + (y2-y1) * t2; 
  end
  visible = 1; % TRUE
end

% gridHyper - calculate a grid for contouring. 
%   Input:
%     rphi           = array (R, phi) ellipse axial ratios
%     kappa          = weighting parameter
%     nnodes         = number of grid nodes, n, in x and y
%     opts.normalize = normalize by n
%   Output:
%     z              = matrix of z values at the nxn grid nodes.
function [z] = gridHyper(rphi, rmax, kappa, nnodes)
  global opts;
  n = length(rphi); % number of data points
  if n < 2 
    return; % error
  end
  z = zeros(nnodes, nnodes);
  s = rToZeta(rmax);
  dx = (2.0 * s) / (nnodes-1);
  dy = dx;
  if opts.normalize
    f = kappa / (n^(1.0/3.0));
  else
    f = kappa;
  end
  % form the data vectors to save time
  h = zeros(n,3);
  for i = 1:n
    h(i,:) = rPhiToH(rphi(i,1), rphi(i,2));
  end
  x = -s;
  for i = 1:nnodes 
    y = -s;
    for j = 1:nnodes
      zsum = 0.0;
      % back-project node to hyberbolic surface
      [rn, pn] = xYToRPhi(x, y, rmax);
      hn = rPhiToH(rn, pn);
      hn = hn * -1.0; % -a for dotH
      for k = 1:n % sum weights
        d = dotH(hn, h(k,:));
        zt = exp(f * (1.0 - d)); % cumulative distribution
        zsum = zsum + zt;
      end % k
      z(i,j) = zsum;
      y = y + dy;
    end % j
    x = x + dx;
  end % i
end

% interpolate - determine linear interpolation point between two nodes. 
% Adopted from Vollmer, 1995.
function [x, y, bool] = interpolate(x1, y1, z1, x2, y2, z2, z0)
  dz1 = z0-z1; 
  dz2 = z0-z2;
  if (dz1 == 0.0) 
    x = x1; 
    y = y1; 
    bool = 1;
  elseif (dz2 == 0.0) 
    x = x2; 
    y = y2; 
    bool = 0;
  elseif (((dz1 > 0.0) && (dz2 > 0.0)) || ((dz1 < 0.0) && (dz2 < 0.0))) 
    x = 0.0; 
    y = 0.0; 
    bool = 0; % FALSE
  else
    dz = z2-z1;
    t = dz1/dz;
    x = x1 + (x2-x1) * t; 
    y = y1 + (y2-y1) * t;
    bool = 1; % TRUE
  end
end

% contourGrid - output one contour level by linear interpolation among grid 
% nodes. Adopted from Vollmer, 1995.
function [lines] = contourGrid(lines, x1, y1, x2, y2, grid, level)
  [ng,mg] = size(grid);
  dnx = (x2-x1)/(ng-1.0); 
  dny = (y2-y1)/(mg-1.0);
  gy1 = y1;
  nx = x1;
  for i = 1:ng-1
    ny = gy1;
    nxp = nx + dnx;
    for j = 1:mg-1
      nyp = ny + dny;
      z1 = grid(i,j); 
      z2 = grid(i+1,j);
      z3 = grid(i+1,j+1); 
      z4 = grid(i,j+1);
      found = 0;
      [x1,y1,bool] = interpolate(nx,ny,z1,nxp,ny,z2,level);
      if bool
        found = found+1;
      end
      [x2,y2,bool] = interpolate(nxp,ny,z2,nxp,nyp,z3,level);
      if bool
        found = found+2;
      end
      [x3,y3,bool] = interpolate(nxp,nyp,z3,nx,nyp,z4,level);
      if bool
        found = found+4;
      end
      [x4,y4,bool] = interpolate(nx,nyp,z4,nx,ny,z1,level);
      if bool
        found = found+8;
      end
      switch (found) 
        case  3 
          lines = cLineOut(lines,x1,y1,x2,y2); 
        case  5
          lines = cLineOut(lines,x1,y1,x3,y3); 
        case  9
          lines = cLineOut(lines,x1,y1,x4,y4); 
        case  6 
          lines = cLineOut(lines,x2,y2,x3,y3); 
        case 10 
          lines = cLineOut(lines,x2,y2,x4,y4); 
        case 12 
          lines = cLineOut(lines,x3,y3,x4,y4); 
        case 15
          d1 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)); 
          d2 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
          d3 = sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4)); 
          d4 = sqrt((x4-x1)*(x4-x1) + (y4-y1)*(y4-y1));
          if ((d1+d3) < (d2+d4)) 
            lines = cLineOut(lines,x1,y1,x2,y2); 
            lines = cLineOut(lines,x3,y3,x4,y4);      
          else 
            lines = cLineOut(lines,x2,y2,x3,y3); 
            lines = cLineOut(lines,x1,y1,x4,y4);
          end
      end % switch
      ny = nyp;
    end % j 
    nx = nxp;
  end % i
end

% lineOut - output a line segment 
function [lines] = lineOut(lines, x1, y1, x2, y2) 
  lines = [lines; [x1,y1,x2,y2]];
end

% cLineOut - output a line segment clipped to current projection. 
function [lines] = cLineOut(lines, x1, y1, x2, y2) 
  global opts;
  if opts.rfp
    lines = [lines; [x1,y1,x2,y2]];
  else
    [cx1, cy1, cx2, cy2, visible] = clipLineCircle(0.0, 0.0, 1.0, x1, y1, x2, y2); 
    if (visible)
      lines = [lines; [cx1,cy1,cx2,cy2]];
    else
      lines = lines;
    end
  end
end

% contour - grids data and outputs contours. 
function [grid, lines] = contour(rphi, rmax, kappa, nlevels, nnodes)
  global opts; 
  grid = gridHyper(rphi, rmax, kappa, nnodes);  
  zmax = max(max(grid));
  x1 = -1.0; 
  y1 = -1.0;
  x2 = 1.0; 
  y2 = 1.0;
  lines = zeros(0,4);
  zinc = zmax/nlevels;
  level = 0.0;
  for i = 1:nlevels-1
    level = level + zinc;
    lines = contourGrid(lines, x1, y1, x2, y2, grid, level);
  end
end

% drawCircle - output a circle, adopted from Rodgers and Adams, 1976, p. 216. 
function [lines] = drawCircle(lines, x, y, radius, n)
  ainc = 2.0 * pi/n;
  c1 = cos(ainc); 
  s1 = sin(ainc);
  x1 = x + radius; 
  y1 = y;
  for i = 0:n
    x2 = x + (x1-x)*c1 - (y1-y)*s1; 
    y2 = y + (x1-x)*s1 + (y1-y)*c1;
    lines = lineOut(lines, x1,y1,x2,y2);
    x1 = x2; 
    y1 = y2;
  end
end

% drawTicks - output projection ticks. 
function [ticks] = drawTicks(rmax)
  global opts;
  if opts.ticks
    ticks = zeros(0,4);
    if opts.rfp
      ts = 0.05;
      ticks = lineOut(ticks, 0.0, 1.0, 0.0, 1.0-ts);
      ticks = lineOut(ticks, 0.0,-1.0, 0.0, -1.0+ts);
      %ticks = lineOut(ticks, -0.5, 1.0, -0.5, 1.0-ts);
      %ticks = lineOut(ticks, 0.5,-1.0, 0.5, -1.0+ts);
      r = 2;
      zm = rToZeta(rmax); 
      while r < rmax
        z = rToZeta(r); 
        t = 2.0 * (z/zm) - 1.0;
        ticks = lineOut(ticks, -1.0, t, -1.0+ts, t);
        ticks = lineOut(ticks, 1.0, t, 1.0-ts, t);
        r = r + 1;
      end
    else % polar
      ts = 0.025;
      %ticks = lineOut(ticks, 1.0, 0.0, 1.0-ts, 0.0);
      %ticks = lineOut(ticks, -1.0, 0.0, -1.0+ts, 0.0);
      %ticks = lineOut(ticks, 0.0, 1.0, 0.0, 1.0-ts);
      %ticks = lineOut(ticks, 0.0,-1.0, 0.0, -1.0+ts);
      ticks = lineOut(ticks, 0.0,-1.0,0.0,1.0);
      ticks = lineOut(ticks,-1.0, 0.0,1.0,0.0);
      r = 1;
      zm = rToZeta(rmax); 
      while r < rmax
        z = rToZeta(r); 
        t = z/zm;
        ticks = lineOut(ticks, t, 0.0+ts, t, 0.0-ts);
        ticks = lineOut(ticks, 0.0+ts, t, 0.0-ts, t);
        if z > 0
          ticks = lineOut(ticks, -t, 0.0+ts, -t, 0.0-ts);
          ticks = lineOut(ticks, 0.0+ts, -t, 0.0-ts, -t);
        end
        r = r + 1;
      end
    end
  end 
end

% drawFrame - output projection frame. 
function [frame] = drawFrame()
  global opts;
  if opts.frame
    frame = zeros(0,4);
    if opts.rfp
      frame = lineOut(frame, -1.0, -1.0, 1.0, -1.0);
      frame = lineOut(frame, 1.0, -1.0, 1.0, 1.0);
      frame = lineOut(frame, 1.0, 1.0, -1.0, 1.0);
      frame = lineOut(frame, -1.0, 1.0, -1.0, -1.0);
    else % polar
      frame = drawCircle(frame, 0.0, 0.0, 1.0, 360);
    end
  end 
end

function [bool] = hasOption(options, option)
  a = strfind(options, option);
  bool = (size(a) > 0);
end

function [opts] = getOptions(options)
  global opts;
  if hasOption(options, 'rad') % radians
    opts.angfmt = 0;
  elseif hasOption(options, 'deg') % degrees
    opts.angfmt = 1;
  elseif hasOption(options, 'grd') % gradians
    opts.angfmt = 2;
  else % radians
    opts.angfmt = 0;
  end  
  if hasOption(options, 'eqd') % equidistant (log R, Elliott plot)
    opts.proj = 0;
  elseif hasOption(options, 'eqa') % equal-area    
    opts.proj = 1;  
  elseif hasOption(options, 'stg') % stereographic   
    opts.proj = 2;  
  elseif hasOption(options, 'ort') % orthographic   
    opts.proj = 3;  
  elseif hasOption(options, 'gno') % gnomic   
    opts.proj = 4;  
  elseif hasOption(options, 'lin') % exponential (linerar R)
    opts.proj = 5;  
  elseif hasOption(options, 'rdl') % radial
    opts.proj = 6;  
  else % equidistant (Elliott plot)
    opts.proj = 0;
  end
  if hasOption(options, 'rfp') 
    opts.rfp = 1;
  else % polar
    opts.rfp = 0;
  end
  if hasOption(options, 'nfr')
    opts.frame = 0;
  else
    opts.frame = 1;
  end  
  if hasOption(options, 'ntc') 
    opts.ticks = 0;
  else
    opts.ticks = 1;
  end
  if hasOption(options, 'ngd')
    opts.grid = 0;
  else
    opts.grid = 1;
  end
  if hasOption(options, 'nnm')
    opts.normalize = 0;
  else
    opts.normalize = 1;  
  end
  if hasOption(options, 'nct')
    opts.contour = 0;
  elseif hasOption(options, 'ctr')
    opts.contour = 1;
  else
    opts.contour = 1;  
  end
  if hasOption(options, 'gi0')
    opts.interp = -1.0;
  elseif hasOption(options, 'gi2')
    opts.interp = 0.5;
  elseif hasOption(options, 'gi3')
    opts.interp = 1.0/3.0;
  elseif hasOption(options, 'gi4')
    opts.interp = 0.25;
  elseif hasOption(options, 'gi5')
    opts.interp = 0.2;
  elseif hasOption(options, 'gi6')
    opts.interp = 1.0/6.0;
  elseif hasOption(options, 'gi8')
    opts.interp = 0.125;
  elseif hasOption(options, 'gia')
    opts.interp = 0.1;
  else
    opts.interp = 0.2;  
  end
end
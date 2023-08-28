function [tout,rout] = ig_rose(theta,x,convert2percent)
%IG_ROSE		- rose
%ROSE   Angle histogram plot.
%   ROSE(THETA) plots the angle histogram for the angles in THETA.  
%   The angles in the vector THETA must be specified in radians.
%
%   ROSE(THETA,N) where N is a scalar, uses N equally spaced bins 
%   from 0 to 2*PI.  The default value for N is 20.
%
%   ROSE(THETA,X) where X is a vector, draws the histogram using the
%   bins specified in X.
%
%   [T,R] = ROSE(...) returns the vectors T and R such that 
%   POLAR(T,R) is the histogram.  No plot is drawn.
%
%   H = ROSE(...) returns a vector of line handles.
%
%   See also HIST, POLAR, COMPASS.

%   Clay M. Thompson 7-9-91
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.9 $  $Date: 1998/09/01 14:47:32 $

if isstr(theta)
        error('Input arguments must be numeric.');
end
theta = rem(rem(theta,2*pi)+2*pi,2*pi); % Make sure 0 <= theta <= 2*pi
if nargin==2,
  x = (0:19)*pi/10+pi/20;

elseif nargin==3,
  if isstr(x)
        error('Input arguments must be numeric.');
  end
  if length(x)==1,
    x = (0:x-1)*2*pi/x + pi/x;
  else
    x = sort(rem(x(:)',2*pi));
  end

end
if isstr(x) | isstr(theta) | isstr(convert2percent)
        error('Input arguments must be numeric.');
end


% Determine bin edges and get histogram
edges = sort(rem([(x(2:end)+x(1:end-1))/2 (x(end)+x(1)+2*pi)/2],2*pi));
edges = [edges edges(1)+2*pi];
nn = histc(rem(theta+2*pi-edges(1),2*pi),edges-edges(1));
nn(end-1) = nn(end-1)+nn(end);
nn(end) = [];



% Form radius values for histogram triangle
if min(size(nn))==1, % Vector
  nn = nn(:); 
end

if convert2percent,
	total_count = sum(nn); % Calculate the total count
    nn = (nn / total_count) * 100;
end


[m,n] = size(nn);
mm = 4*m;
r = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

% Form theta values for histogram triangle from triangle centers (xx)
zz = edges;

t = zeros(mm,1);
t(2:4:mm) = zz(1:m);
t(3:4:mm) = zz(2:m+1);


if nargout<2
  h = polar(t,r);
  if nargout==1, tout = h; end
  return
end

if min(size(nn))==1,
  tout = t'; rout = r';
else
  tout = t; rout = r;
end


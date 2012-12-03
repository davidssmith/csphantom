function I = csphantom (n, bandwidth, signalToNoiseRatio, b1scale, imgReal)
%%CSPHANTOM  Test image for compressive sensing reconstructions.
%
%  IMG = CSPHANTOM(N) returns an NxN phantom with various features designed
%  to fail a Fourier-based compressed sensing reconstruction.
%
%  IMG = CSPHANTOM(N,BW) returns an NxN phantom with sampling bandwidth
%  equal to BW*N. This more accurately simulates what occurs when a
%  physical phantom is sampled in an MRI scanner. The default value of BW
%  is Inf. BW can be thought of as the intrinsic phantom pixel size divided 
%  by the Fourier measurement resolution. This parameter affects the
%  sharpness of small features and the level of Gibbs ringing visible.
%
%  Copyright 2012 David S. Smith <david.smith@vanderbilt.edu>.
%
%  This file is part of CSPHANTOM.
% 
%  CSPHANTOM is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  CSPHANTOM is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with CSPHANTOM.  If not, see <http://www.gnu.org/licenses/>.
% 


if nargin < 1; n = 256; end  % default size if none given
if nargin < 2; bandwidth = Inf; end
if nargin < 3; signalToNoiseRatio = Inf; end
if nargin < 4; b1scale = 1; end
if nargin < 5; imgReal = true; end


if ~isinf(bandwidth)
  n = n * bandwidth; % temporarily upsample the phantom creation
end

if n < 8, error('Smallest sensible phantom is 8 x 8.'); end

% quadrant sizes - slightly less than one fourth of the area, to give room 
% for a border
quadrantSize = round((n - round(0.1*n))/2);

% Set up background shading: low amplitude sinc, a la B1 inhomogeneities
bgLevel = 0.5;  % central background gray level
[X, Y] = meshgrid(1:(2*quadrantSize), 1:(2*quadrantSize));
R = sqrt((X-quadrantSize).^2 + (Y-quadrantSize).^2);
BG = bgLevel * cos(b1scale*R/(0.4*pi*quadrantSize));

% go back to quadrant coordinates
[X, Y] = meshgrid(1:quadrantSize, 1:quadrantSize);

% QUADRANT II: low-contrast circles
nc = 4;  % # rows and columns of circles
Q2 = BG(1:quadrantSize,1:quadrantSize);
inc = quadrantSize / (nc+1);
for j = 1:nc
	rad = quadrantSize / (11 + j^2);
	for k = 1:nc
		val  = 0.9*bgLevel - 0.1/k;
		rctr = round(j*inc);
		cctr = round(k*inc);
		Z2 = (X-cctr).^2 + (Y-rctr).^2 <= rad^2;
		Q2(Z2) = val;
	end
end

% QUADRANT I: diagonal ramp
Q1 = BG(1:quadrantSize,(quadrantSize+1):(2*quadrantSize));
scale = quadrantSize/5;
ctr = round(quadrantSize/2);
Z = abs(X-ctr) + abs(Y-ctr);
Z1 = Z <= scale;
Z2 = Z > scale & Z <= 2*scale;
Q1(Z1) = Q1(Z1) * 2 .* (1 - Z(Z1) / scale);
Q1(Z2) = Q1(Z2) .* (Z(Z2) / scale - 1);
d = round(3*quadrantSize/16);
r = [d quadrantSize-d];
c = [d quadrantSize-d];
w = round(d/4);
for j = 1:2
	for k = 1:2
		Q1(r(j)-w:r(j)+w,c(k)-w:c(k)+w) = Q1(r(j)-w:r(j)+w,c(k)-w:c(k)+w)*(1+k*j/64);
	end
end

% QUADRANT III: Gaussian bumps and quadratic hole
ctr = round(quadrantSize/2);
scale = quadrantSize / 4;
R2 = (X-ctr).^2 + (Y-ctr).^2;
Q3 = BG((quadrantSize+1):(2*quadrantSize),1:quadrantSize);
Q3withHole = Q3 .* R2 / scale^2;
Q3 = min(Q3,Q3withHole);
ctr = round([quadrantSize/5 4*quadrantSize/5]);
scale = scale / 6;
for r = 1:2
	for c = 1:2
		Q3 = Q3 + (3-r)*c*0.1*bgLevel * exp(-((X-ctr(r)).^2+(Y-ctr(c)).^2)/scale^2);
	end
end

% QUADRANT IV: line pairs
Q4 = BG((quadrantSize+1):(2*quadrantSize),(quadrantSize+1):(2*quadrantSize));
borderThick = round(quadrantSize/8); % define a border around the line pair region
lineLength = max(1,borderThick); 
Q4floor = 0.1;
%
% HORIZONTAL LINES
%
lineSep = 1;
lineThick = 1;
nlines = 3;
sepIncr = 1;
if ~isinf(bandwidth)
  lineSep = lineSep*round(bandwidth);
  lineThick = lineThick*round(bandwidth);
  sepIncr = sepIncr*round(bandwidth);
end
col = 1;
row = 1 + round(0.5*borderThick) + lineLength;
while (row + lineThick) <= quadrantSize - borderThick
  for t = 1:nlines
    Q4(row:row+lineThick-1,col:col+lineLength-1) = Q4floor;
    row = row + lineThick + lineSep;
    if (row + lineThick) > quadrantSize - borderThick, break; end
  end
  row = row + lineSep + lineThick;
	lineSep = lineSep + sepIncr;
  lineThick = lineThick + sepIncr;
end
%
% VERTICAL LINES
%
lineSep = 1;
lineThick = 1;
sepIncr = 1;
if ~isinf(bandwidth)
  lineSep = lineSep*round(bandwidth);
  lineThick = lineThick*round(bandwidth);
  sepIncr = sepIncr*round(bandwidth);
end
col = 1 + round(0.5*borderThick) + lineLength;
row = 1;
while (col + lineThick) <= quadrantSize - borderThick  % vertical lines
  %Q4(row:row+lineLength-1,col:col+lineThick-1) = Q4floor;
  for t = 1:nlines
    Q4(row:row+lineLength-1,col:col+lineThick-1) = Q4floor;
    col = col + lineThick + lineSep;
    if (col + lineThick) > quadrantSize - borderThick, break; end
  end
	col = col + lineSep + lineThick;
	lineSep = lineSep + sepIncr;
  lineThick = lineThick + sepIncr;
end
%
% CIRCLES
%
nlines = 3;
ctr = round(quadrantSize/2);  % concentric circles
R = sqrt((X-ctr).^2 + (Y-ctr).^2);
circThick = 1;
circGap = 3;
circIncr = 1;
if ~isinf(bandwidth)
  circThick = circThick*round(bandwidth);
  circGap = circGap*round(bandwidth);
  circIncr = circIncr*round(bandwidth);
end
r = circThick + circGap;
white = false;
while (ctr + r) < (quadrantSize - 1.5*borderThick)
  for t = 1:nlines
    rim = abs(R - r) <= 2*circThick;
    if white
      Q4(rim) = Q4(rim).*(1+Q4floor*exp(-2*(R(rim)-r).^2/circThick^2));
    else
      Q4(rim) = Q4(rim).*(1-exp(-2*(R(rim)-r).^2/circThick^2));
    end
    r = r + circThick + circGap;
    if (ctr + r) > (quadrantSize - 1.2*borderThick), break; end
  end
  if white
    white = false;
  else
    white = true;
  end
  r = r + circGap;
end





% pad area of signal with 'air'
npadf = (n - 2*quadrantSize) / 2;
npad = round(npadf);
I = zeros(n,n);
I(npad+1:npad+2*quadrantSize,npad+1:npad+2*quadrantSize) = [Q2 Q1; Q3 Q4];


% downsample to specified sampling bandwidth
if ~isinf(bandwidth)
  F = fft2(I);
  n = n / bandwidth;  % truncate Fourier data
  F = F([1:ceil(n/2) end-floor(n/2)+1:end],[1:ceil(n/2) end-floor(n/2)+1:end]);
  % ringing filter
  W = ifftshift(genRiesz([n n],0.6,0.6));
  %W = ifftshift(tukey(n)*tukey(n)');
  F = F .* W; % apply window
  I = abs(ifft2(F)) / bandwidth^2;
end

% finally add some Gaussian complex noise, if SNR .ne. Inf
if ~isinf(signalToNoiseRatio)
  mu = 0;
  sigma = sqrt(0.5^2 / (signalToNoiseRatio^2 + 2)) / sqrt(2);
  I = I + sigma*randn(size(I)) + 1i*sigma*randn((size(I)));
end

if imgReal
  I = abs(I);
end


function w = tukey (N)
%%WINDOW  Tukey window
alpha = 0.8;
w = ones(N,1);
n = -(N-1)/2 : -alpha*N/2;
L = length(n);
w(1:L) = 0.5*(1+cos(pi*(abs(n)-alpha*N/2)/((1-alpha)*N/2)));
w(N : -1 : N-L+1) = w(1:L);


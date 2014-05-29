function I = csphantom (n, bandwidth, signalToNoiseRatio, b1scale, img_real)
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
if nargin < 5; img_real = true; end


if ~isinf(bandwidth)
  n = n * bandwidth; % temporarily upsample the phantom creation
end

if n < 8, error('Smallest sensible phantom is 8 x 8.'); end

% quadrant sizes - slightly less than one fourth of the area, to give room 
% for a border
quadrant_size = round((n - round(0.1*n))/2);

% Set up background shading: low amplitude sinc, a la B1 inhomogeneities
bglevel = 0.5;  % central background gray level
[X, Y] = meshgrid(1:(2*quadrant_size), 1:(2*quadrant_size));
R = sqrt((X-quadrant_size).^2 + (Y-quadrant_size).^2);
BG = bglevel * besselj(0,b1scale*2*R/max(R(:)));

% go back to quadrant coordinates
[X, Y] = meshgrid(1:quadrant_size, 1:quadrant_size);

% QUADRANT II: low-contrast circles
nc = 4;  % # rows and columns of circles
Q2 = BG(1:quadrant_size,1:quadrant_size);
inc = quadrant_size / (nc+1);
for j = 1:nc
	rad = quadrant_size / (11 + j^2);
	for k = 1:nc
		val  = 0.9*bglevel - 0.1/k;
		rctr = round(j*inc);
		cctr = round(k*inc);
		Z2 = (X-cctr).^2 + (Y-rctr).^2 <= rad^2;
		Q2(Z2) = val;
	end
end

% QUADRANT I: diagonal ramp
Q1 = BG(1:quadrant_size,(quadrant_size+1):(2*quadrant_size));
scale = quadrant_size/5;
ctr = round(quadrant_size/2);
Z = abs(X-ctr) + abs(Y-ctr);
Z1 = Z <= scale;
Z2 = Z > scale & Z <= 2*scale;
Q1(Z1) = Q1(Z1) * 2 .* (1 - Z(Z1) / scale);
Q1(Z2) = Q1(Z2) .* (Z(Z2) / scale - 1);
d = round(3*quadrant_size/16);
r = [d quadrant_size-d];
c = [d quadrant_size-d];
w = round(d/4);
for j = 1:2
	for k = 1:2
		Q1(r(j)-w:r(j)+w,c(k)-w:c(k)+w) = Q1(r(j)-w:r(j)+w,c(k)-w:c(k)+w)*(1+k*j/64);
	end
end

% QUADRANT III: Gaussian bumps and quadratic hole
ctr = round(quadrant_size/2);
scale = quadrant_size / 4;
R2 = (X-ctr).^2 + (Y-ctr).^2;
Q3 = BG((quadrant_size+1):(2*quadrant_size),1:quadrant_size);
Q3withHole = Q3 .* R2 / scale^2;
Q3 = min(Q3,Q3withHole);
ctr = round([quadrant_size/5 4*quadrant_size/5]);
scale = scale / 6;
for r = 1:2
	for c = 1:2
		Q3 = Q3 + (3-r)*c*0.1*bglevel * exp(-((X-ctr(r)).^2+(Y-ctr(c)).^2)/scale^2);
	end
end

% QUADRANT IV: line pairs
Q4 = BG((quadrant_size+1):(2*quadrant_size),(quadrant_size+1):(2*quadrant_size));
border_thick = round(quadrant_size/8); % define a border around the line pair region
line_length = max(1,border_thick); 
Q4floor = 0.1;
%
% HORIZONTAL LINES
%
line_sep = 1;
line_thick = 1;
nlines = 3;
sep_incr = 1;
if ~isinf(bandwidth)
  line_sep = line_sep*round(bandwidth);
  line_thick = line_thick*round(bandwidth);
  sep_incr = sep_incr*round(bandwidth);
end
col = 1;
row = 1 + round(0.5*border_thick) + line_length;
while (row + line_thick) <= quadrant_size - border_thick
  for t = 1:nlines
    Q4(row:row+line_thick-1,col:col+line_length-1) = Q4floor;
    row = row + line_thick + line_sep;
    if (row + line_thick) > quadrant_size - border_thick, break; end
  end
  row = row + line_sep + line_thick;
	line_sep = line_sep + sep_incr;
  line_thick = line_thick + sep_incr;
end
%
% VERTICAL LINES
%
line_sep = 1;
line_thick = 1;
sep_incr = 1;
if ~isinf(bandwidth)
  line_sep = line_sep*round(bandwidth);
  line_thick = line_thick*round(bandwidth);
  sep_incr = sep_incr*round(bandwidth);
end
col = 1 + round(0.5*border_thick) + line_length;
row = 1;
while (col + line_thick) <= quadrant_size - border_thick  % vertical lines
  %Q4(row:row+line_length-1,col:col+line_thick-1) = Q4floor;
  for t = 1:nlines
    Q4(row:row+line_length-1,col:col+line_thick-1) = Q4floor;
    col = col + line_thick + line_sep;
    if (col + line_thick) > quadrant_size - border_thick, break; end
  end
	col = col + line_sep + line_thick;
	line_sep = line_sep + sep_incr;
  line_thick = line_thick + sep_incr;
end
%
% CIRCLES
%
nlines = 3;
ctr = round(quadrant_size/2);  % concentric circles
R = sqrt((X-ctr).^2 + (Y-ctr).^2);
circ_thick = 1;
circ_gap = 3;
circIncr = 1;
if ~isinf(bandwidth)
  circ_thick = circ_thick*round(bandwidth);
  circ_gap = circ_gap*round(bandwidth);
  circIncr = circIncr*round(bandwidth);
end
r = circ_thick + circ_gap;
white = false;
while (ctr + r) < (quadrant_size - 1.5*border_thick)
  for t = 1:nlines
    rim = abs(R - r) <= 2*circ_thick;
    if white
      Q4(rim) = Q4(rim).*(1+Q4floor*exp(-2*(R(rim)-r).^2/circ_thick^2));
    else
      Q4(rim) = Q4(rim).*(1-exp(-2*(R(rim)-r).^2/circ_thick^2));
    end
    r = r + circ_thick + circ_gap;
    if (ctr + r) > (quadrant_size - 1.2*border_thick), break; end
  end
  if white
    white = false;
  else
    white = true;
  end
  r = r + circ_gap;
end

% pad area of signal with 'air'
npadf = (n - 2*quadrant_size) / 2;
npad = round(npadf);
I = zeros(n);
I(npad+1:npad+2*quadrant_size,npad+1:npad+2*quadrant_size) = [Q2 Q1; Q3 Q4];

% downsample to specified sampling bandwidth
if ~isinf(bandwidth)
  F = fft2(I);
  n = n / bandwidth;  % truncate Fourier data
  F = F([1:ceil(n/2) end-floor(n/2)+1:end],[1:ceil(n/2) end-floor(n/2)+1:end]);
  % ringing filter
  W = ifftshift(tukey(n)*tukey(n)');
  F = F .* W; % apply window
  I = abs(ifft2(F)) / bandwidth^2;
end

% finally add some Gaussian complex noise, if SNR .ne. Inf
if ~isinf(signalToNoiseRatio)
  sigma = sqrt(0.125 / (signalToNoiseRatio^2 + 1));
  I = I + sigma*randn(size(I)) + 1i*sigma*randn((size(I)));
end

if img_real, I = abs(I); end


function w = tukey (N)
%%WINDOW  Tukey window
alpha = 0.8;
w = ones(N,1);
n = -(N-1)/2 : -alpha*N/2;
L = length(n);
w(1:L) = 0.5*(1+cos(pi*(abs(n)-alpha*N/2)/((1-alpha)*N/2)));
w(N : -1 : N-L+1) = w(1:L);


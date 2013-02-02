function example
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

%% BASIC USAGE
I = csphantom;
subplot(211);
imagesc(I);
colormap(gray);
axis image;
title 'basic usage'

%% FANCY USAGE with customized settings
n = 512;
bw = 4;
snr = 20;
b1scale = 0.8;
imgReal = false;
I = csphantom(n,bw,snr,b1scale,imgReal);
subplot(212);
imshow(abs(I));
colormap(gray);
title 'advanced usage'
axis image;


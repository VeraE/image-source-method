function random_shifts = create_random_shifts(N,shift)
%CREATE_RANDOM_SHIFTS creates random shifts of up to the value 'shift' for image
% sources with an order > 3 [BM09]
%
%   Input parameters:
%       N               :  index for finite summation in eq. (10) in [AB79]
%       shift           :  max. shift in m
%
%   Output parameters:
%       random shifts   :  random shifts for image sources / m
%                          [2 x 2 x 2 x 2*N+1 x 2*N+1 x 2*N+1 x 3], last
%                          dimension for values for x, y and z
%
%   [AB79]  J. B. Allen and D. A. Berkley (1979): Image source method for
%           efficiently simulating small-room acoustics. J. Acoust. Soc. Am.
%           65(4):943-950
%   [BM09]  C. Borﬂ and R. Martin (2009): An Improved Parametric Model for
%           Perception-Based Design of Virtual Acoustics. Proc. of the 35th Int.
%           Audio Eng. Soc. Conf.

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2020 Vera Erbes                                              *
%                                                                            *
% Permission is hereby granted, free of charge, to any person obtaining a copy  *
% of this software and associated documentation files (the "Software"), to deal *
% in the Software without restriction, including without limitation the rights  *
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
% copies of the Software, and to permit persons to whom the Software is         *
% furnished to do so, subject to the following conditions:                      *
%                                                                               *
% The above copyright notice and this permission notice shall be included in    *
% all copies or substantial portions of the Software.                           *
%                                                                               *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     *
% THE SOFTWARE.                                                                 *
%********************************************************************************

%create random shifts for all sources
randx = (rand(2,2,2,2*N+1,2*N+1,2*N+1)*2-1)*shift;
randy = (rand(2,2,2,2*N+1,2*N+1,2*N+1)*2-1)*shift;
randz = (rand(2,2,2,2*N+1,2*N+1,2*N+1)*2-1)*shift;

%matrices of summation indices for image source order calculation
[u,v,w,l,m,n] = ndgrid(0:1,0:1,0:1,-N:N,-N:N,-N:N);
%order of image sources
order = abs(l-u)+abs(m-v)+abs(n-w)+abs(l)+abs(m)+abs(n);
%no randomisation for image sources up to order 3
idx = find(order<=3);
randx(idx) = 0;
randy(idx) = 0;
randz(idx) = 0;

%concatenate random shifts
random_shifts = cat(7,randx,randy,randz);
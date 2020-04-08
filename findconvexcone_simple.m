function [idx,weights] = findconvexcone_simple(x0,xs)
%This is a simplified version of the function findconvexcone() that is contained
%in the Sound Field Synthesis Toolbox, see description and license below. The
%simplification was realised to make it more efficient when used together with an
%implementation of the image source model for selection of HRTF directions, but as
%a consequence the function now only works for full-spherical point clouds (which
%is also sensible when using it for the image source model).
%
%FINDCONVEXCONE selects up to 3 points from x0 with xs in their conic span
%
%   Usage: [idx,weights] = findconvexcone(x0,xs)
%
%   Input parameters:
%       x0          - point cloud on a sphere around the origin / m [nx3]
%       xs          - desired direction as point in R^3 / m [1x3]
%
%   Output parameters:
%       idx         - row indices of N points in x0 [Nx1]
%                     where N is 1,2 or 3
%       weights     - weights [Nx1]
%
%   FINDCONVEXCONE(x0,xs) returns 1,2 or 3 row indices into x0 and non-negative
%   weights w1, ..., w3 such that w1*x1 + w2*x2 + w3*x3 with
%   [x1; x2; x3] == x0(idx,:) composes the point inside the triangle spanned
%   by x1, x2, x3.
%
%   x1...x3 are selected from the convex hull in R3.
%   Various precautions are taken to make this well-behaved in most cases.
%
%   (If all x0 and xs have unit norm this is VBAP.)
%
%   See also: findnearestneighbour, test_interpolation_point_idx

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Prepare Grid (see local functions below) ========================
% Normalise x0 and xs to lie on unit sphere as only direction is relevant
xs = xs./norm(xs,2);
radii = vector_norm(x0,2);
x0 = x0./repmat(radii,[1,size(x0,2)]);


%% ===== Computation =====================================================
% Delaunay triangulation of convex hull
simplices = convhulln(x0);

% Find x0 with smallest angle to xs
[~,most_aligned_point] = ...
    max(vector_product(x0,repmat(xs,size(x0,1),1),2));

% The simplices at "most aligned point" are the most likely candidates,
% put them at the beginning of the list
mask = logical(sum(simplices==most_aligned_point,2));
simplices = [simplices(mask,:); simplices(~mask,:) ];

% One of these simplices contains xs
for n = 1:size(simplices,1);
    A = x0(simplices(n,:),:);
    weights = xs/A;
    weights(abs(weights)<1e-10) = 0;
    if all(weights >= 0) % non-negative weights == convex combination
        idx = simplices(n,:);
        break;
    end
end

% Normalise weights
weights = weights/sum(weights);

[weights,order] = sort(weights.','descend');
idx = idx(order).';
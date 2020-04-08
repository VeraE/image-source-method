function h = calculate_ISM(s,r,L,beta,N,fs,receiver_type,cutIS,...
    downsampling,random_shifts,listener_dir,predelay,listener_view)
%CALCULATE_ISM calculates the impulse response for a rectangular room with the
% image source model [AB79,BEW18]
%
%   Input parameters:
%   mandatory:
%       s               :  source position in R^3 / m [1x3]
%       r               :  receiver position in R^3 / m [1x3]
%       L               :  room dimensions / m [1x3]
%       beta            :  reflection coefficients of the 6 walls, can be given
%                          in multiple versions X to be calculated at once [3x2xX]
%       N               :  index for finite summation in eq. (10) in [AB79]
%       fs              :  sampling rate / samples/s
%       receiver_type   :  'monaural' or 'binaural'
%       cutIS           :  1: Image sources that do not form a complete grid after
%                          a certain max. tau are left out, 0: all image sources
%                          are calculated
%       downsampling    :  1: h is downsampled to 44.1 kHz (monaural) or to 
%                          original sampling rate of HRIRs (binaural), 0: no
%                          downsampling
%       random_shifts   :  random shifts for image sources / m
%                          [2 x 2 x 2 x 2*N+1 x 2*N+1 x 2*N+1 x 3], last
%                          dimension for values for x, y and z
%   optional (for binaural receiver):
%       listener_dir    :  contains HRIRs as SOFA file (www.sofaconventions.org)
%                          [AES15], source positions have to be given in
%                          [deg,deg,m]
%       predelay        :  predelay contained in HRIRs / samples
%       listener_view   :  vector indicating listener view direction, only azimuth
%                          is taken into account! / m [1x3]
%
%   Output parameters:
%       h               :  impulse response(s), [lengthxX] (monaural) or
%                          [lengthx2xX] (binaural)
%
%   [AB79]   J. B. Allen and D. A. Berkley (1979): Image source method for
%            efficiently simulating small-room acoustics. J. Acoust. Soc. Am.
%            65(4):943-950
%   [AES15]  AES69-2015: AES standard for file exchange - Spatial acoustic data
%            file format, Audio Eng. Soc., Inc
%   [BEW18]  F. Brinkmann, V. Erbes and S. Weinzierl (2018): Extending the closed
%            form image source model for source directivity. Proc. of the 44th
%            German Annual Conf. on Acoustics (DAGA)

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

%definitions
c = 343; %speed of sound / m/s

%number of image sources
number_of_IS = 2^3*(2*N+1)^3;

%matrices of summation indices
[u,v,w,l,m,n] = ndgrid(0:1,0:1,0:1,-N:N,-N:N,-N:N);

%order of image sources
order = abs(l-u)+abs(m-v)+abs(n-w)+abs(l)+abs(m)+abs(n);

%entries of vectors pointing from receiver to image sources
phatx = (1-2*u)*s(1)+2*l*L(1)-r(1);
phaty = (1-2*v)*s(2)+2*m*L(2)-r(2);
phatz = (1-2*w)*s(3)+2*n*L(3)-r(3);

%randomise image source positions
phatx = phatx+random_shifts(:,:,:,:,:,:,1);
phaty = phaty+random_shifts(:,:,:,:,:,:,2);
phatz = phatz+random_shifts(:,:,:,:,:,:,3);

%absolute positions of image sources
px = phatx+r(1);
py = phaty+r(2);
pz = phatz+r(3);

%distances between image sources and receiver
d = sqrt(phatx.^2+phaty.^2+phatz.^2);

%distances between image sources and source (for cutIS)
dsource = sqrt((px-s(1)).^2+(py-s(2)).^2+(pz-s(3)).^2);

%number of beta versions
X = size(beta,3);

%amplitudes of image sources
A = zeros([size(u) X]); %preallocation
otherdims = repmat({':'},[1 ndims(u)]);
for k = 1:X
    A(otherdims{:},k) = ...
        beta(1,1,k).^abs(l-u).*beta(2,1,k).^abs(m-v).*beta(3,1,k).^abs(n-w)...
        .*beta(1,2,k).^abs(l).*beta(2,2,k).^abs(m).*beta(3,2,k).^abs(n)...
        ./(4*pi*d);
end

%delay in samples, quantised to integer samples, Matlab indexing
tau = round(d/c*fs)+1;

%vectorise multidimensional arrays
tau = tau(:);
A = reshape(A,[length(tau) X]);
phatx = phatx(:);
phaty = phaty(:);
phatz = phatz(:);
px = px(:);
py = py(:);
pz = pz(:);
d = d(:);
dsource = dsource(:);
order = order(:);

if cutIS
    %find max. distance around source until image sources still form a complete
    %grid
    [~,idx] = min(L);
    dcut = 2*N*L(idx);
    %delete all image sources outside the complete grid
    idx = find(dsource>dcut);
    tau(idx) = [];
    A(idx,:) = [];
    phatx(idx) = [];
    phaty(idx) = [];
    phatz(idx) = [];
    px(idx) = [];
    py(idx) = [];
    pz(idx) = [];
    d(idx) = [];
    dsource(idx) = [];
    order(idx) = [];
end

%calculate RIRs & BRIRs
if strcmp(receiver_type,'monaural')
    %maximum delay in samples for length of impulse response
    taumax = max(tau);
    h = zeros(taumax,X); %preallocation
    for k = 1:length(tau)
        %add amplitudes to impulse response
        h(tau(k),:) = h(tau(k),:)+A(k,:);
    end
    if downsampling
        %get resampling factor
        [p,q] = rat(fs/44100);
        %zeropad to not lose samples when resampling last dirac pulse
        %(number 10 seems to be hardcoded in Matlab)
        h = [h; zeros(ceil(10*p/q),X)];
        h = resample(h,q,p)*p/q;
    end
elseif strcmp(receiver_type,'binaural')    
    %resampling of HRIRs
    %get resampling factor
    [p,q] = rat(fs/listener_dir.Data.SamplingRate);
    %resample all HRIRs (and keep sofa data format)
    hrirs(:,1,:) = resample(squeeze(listener_dir.Data.IR(:,1,:))',p,q)';
    hrirs(:,2,:) = resample(squeeze(listener_dir.Data.IR(:,2,:))',p,q)';
    
    %maximum delay in samples for length of impulse response
    length_hrirs = size(hrirs,3);
    taumax = max(tau)+length_hrirs-1-round(predelay*p/q);
    
    %shift tau by predelay
    tau_shift = tau-round(predelay*p/q);
    
    %calculate incidence angles...
    phi_r = atan2(phaty,phatx); %in [-pi,pi]
    theta_r = asin(phatz./d); %in [-pi/2,pi/2]
    %...and construct image source directions relative to listener
    ISpositions = [phi_r,theta_r,ones(length(tau),1)];
    
    %turn HRIR source positions according to listener view direction
    %calculate listener view azimuth
    phi_view = atan2(listener_view(2),listener_view(1)); %in [-pi,pi]
    %turn HRIR source positions and ensure correct range of angles in rad
    HRIR_SourcePositions(:,1) = ...
        mod(listener_dir.SourcePosition(:,1)*pi/180+phi_view,2*pi);
    HRIR_SourcePositions(:,2) = listener_dir.SourcePosition(:,2)*pi/180;
    HRIR_SourcePositions(:,3) = listener_dir.SourcePosition(:,3);
    
    %convert source positions to cartesian coordinates for interpolation and
    %nearest neighbour search
    %conversion of HRIR source positions (radius normalised)
    [HRIR_SourcePositions(:,1),HRIR_SourcePositions(:,2),...
        HRIR_SourcePositions(:,3)] = sph2cart(...
        HRIR_SourcePositions(:,1),HRIR_SourcePositions(:,2),...
        ones(size(HRIR_SourcePositions,1),1));
    %conversion of image source directions (radius already normalised)
    [ISpositions(:,1),ISpositions(:,2),ISpositions(:,3)] = sph2cart(...
        ISpositions(:,1),ISpositions(:,2),ISpositions(:,3));
    
    %find image sources up to order 3 for interpolation of HRIRs (rest with
    %nearest neighbour search)
    idx = find(order<=3);
    tau_order3 = tau(idx);
    tau_shift_order3 = tau_shift(idx);
    ISpositions_order3 = ISpositions(idx,:);
    A_order3 = A(idx,:);
    %delete image sources up to order 3 in main vectors
    tau(idx) = [];
    tau_shift(idx) = [];
    ISpositions(idx,:) = [];
    A(idx,:) = [];
    
    %image sources calculated with interpolation of HRIRs
    h = zeros(taumax,2,X); %preallocation
    %struct for SFS Toolbox
    conf = struct;
    conf.ir.useinterpolation = true;
    conf.ir.interpolationmethod = 'timedomain';
    for k = 1:length(tau_order3)
        %interpolation
        [idx,weights] = findconvexcone_simple(...
            HRIR_SourcePositions,ISpositions_order3(k,:));
        hrirs_selected = squeeze(...
            interpolate_ir(hrirs(idx,:,:),weights,conf))';
        
        %weight selected HRIRs by corresponding amplitudes for all beta versions
        hrirs_selected_X = zeros([size(hrirs_selected) X]); %preallocation
        for kk = 1:X
            hrirs_selected_X(:,:,kk) = hrirs_selected*A_order3(k,kk);
        end
        %add selected and weighted HRIRs to h
        h(tau_shift_order3(k):tau_shift_order3(k)+length_hrirs-1,:,:) = ...
            h(tau_shift_order3(k):tau_shift_order3(k)+length_hrirs-1,:,:) + ...
            hrirs_selected_X;
    end
    
    %image sources calculated with nearest neighbour search of HRIRs
    for k = 1:length(tau)
        %find nearest neighbour
        %calculate distance between points
        distance = sqrt(sum(abs(...
            bsxfun(@minus,HRIR_SourcePositions.',ISpositions(k,:).')).^2,1));
        %find minimum distance (first occurring index)
        [~,idx] = min(distance);
        %selected HRIRs
        hrirs_selected = squeeze(hrirs(idx,:,:))';
        
        %weight selected HRIRs by corresponding amplitudes for all beta versions
        hrirs_selected_X = zeros([size(hrirs_selected) X]); %preallocation
        for kk = 1:X
            hrirs_selected_X(:,:,kk) = hrirs_selected*A(k,kk);
        end
        %add selected and weighted HRIRs to h
        h(tau_shift(k):tau_shift(k)+length_hrirs-1,:,:) = ...
            h(tau_shift(k):tau_shift(k)+length_hrirs-1,:,:) + ...
            hrirs_selected_X;
    end
    %Point source factor 4*pi is supposed to be contained in the measured HRIRs
    %and is therefore corrected:
    h = h*4*pi;
    
    %downsample final impulse response
    if downsampling
        h_tmp = zeros(ceil((size(h,1)+ceil(10*p/q))/p*q),2,X); %preallocation
        for k = 1:X
            %zeropad to not lose samples when resampling last dirac pulse
            %(number 10 seems to be hardcoded in Matlab)
            h_tmp(:,:,k) = resample([h(:,:,k); zeros(ceil(10*p/q),2)],q,p);
        end
        h = h_tmp;
    end
end
%This script demonstrates the use of the function calculate_ISM(), confer the
%information given in the function header for details on the parameters.
%
%You need a full-spherical HRTF dataset given in the format SOFA
%(www.sofaconventions.org) [AES15], e.g. download the FABIAN HRTF dataset for a
%head-above-torso orientation of 0° at http://doi.org/10.14279/depositonce-5718.3
%[BLW+17] and save it as 'HRIR_dataset.sofa'.
%
%For higher orders of N, this script may take a considerable amount of time to
%execute depending on your system, so be careful when testing different
%parameters.
%
%   [AES15]   AES69-2015: AES standard for file exchange - Spatial acoustic data
%             file format, Audio Eng. Soc., Inc
%   [BLW+17]  F. Brinkmann, A. Lindau, S. Weinzierl, S. van de Par, M. Müller-
%             Trapet, R. Opdam and M. Vorländer (2017): A High Resolution and
%             Full-Spherical Head-Related Transfer Function Database for Different
%             Head-Above-Torso Orientations. J. Audio Eng. Soc. 65(10):841-848

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

close all; clear; clc

% SFS_start
% SOFAstart

%% Specify geometry and parameters
s = [4.93 5.37 1.6];
r = [4.93 2.87 1.6];
L = [10 7 3];
beta(:,:) = ones(3,2)*.9;
N = 1;
fs = 441000;
cutIS = 1;
downsampling = 1; 

%% Create random shifts for image source orders greater than 3
shift = 1; %max. shift / m
random_shifts = create_random_shifts(N,shift);

%% Calculate RIR
receiver_type = 'monaural';

%calculate ir
RIR = calculate_ISM(...
    s,r,L,beta,N,fs,receiver_type,cutIS,downsampling,random_shifts);

%plot
figure
    plot((0:size(RIR,1)-1)/44100,RIR)
    grid
    xlabel('time / s'), ylabel('amplitude')
    title('RIR')

%% Calculate BRIRs
receiver_type = 'binaural';
%load HRIR dataset given in SOFA
listener_dir = SOFAload('HRIR_dataset.sofa');
%predelay of HRIRs until start of impulse response for accurate time alignment 
predelay = 26; %/ samples
%listener view
azimuth = 90*pi/180; %receiver looking towards the source in this example
listener_view = [cos(azimuth) sin(azimuth) zeros(size(azimuth))];

%calculate irs
BRIRs = calculate_ISM(...
    s,r,L,beta,N,fs,receiver_type,cutIS,downsampling,random_shifts,...
    listener_dir,predelay,listener_view);

%plot
figure
    plot((0:size(BRIRs,1)-1)/listener_dir.Data.SamplingRate,BRIRs)
    grid
    xlabel('time / s'), ylabel('amplitude')
    legend('left','right')
    title('BRIRs')
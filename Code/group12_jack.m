clear
close all
clc

Tw = 30; % s - Lenght of each window, smaller => More averages, less freq res
overlap = .5; % -% like in lab 4

% I think we'll measure for a long time a random input response
% Then we split the response in pieces of T lenght, apply windows and go in
% frequency domain
% Compute expected value of power and cross spectra and estimate the FRF

% Don't know which of these
Data = load('modal_extraction.mat');
input = Data.input(:); % Vector of the overall measurement
output = Data.output(:);

% [n,p]=uigetfile('.mat','choose the data file');
% cd(p);
% load(n);

fs = Data.fs;
dt = 1/fs;
N = length(output); % Any vector which is the long measurement

t_vec = 0:dt:(N-1)*dt; % Global time vector

figure
plot (t_vec, output) % This should plot the entire response in time
grid on

N_samples = round(Tw/dt);             % Number of elements of each split
N_overlap = floor(overlap*N_samples); % Numer of elements overlapping

win = hanning(N_samples);

[Gxx,freq_xx]=autocross(input,input,fs,N_samples,N_overlap,win);
[Gyy,freq_yy]=autocross(output,output,fs,N_samples,N_overlap,win); 
[Gxy,freq_xy]=autocross(input,output,fs,N_samples,N_overlap,win); 
[Gyx,freq_yx]=autocross(output,input,fs,N_samples,N_overlap,win); 
coherence = abs(Gxy).^2./(Gxx.*Gyy);

figure
subplot(1,2,1)
semilogy(freq_xx, Gxx)
grid on

subplot(1,2,2)
semilogy(freq_yy, Gyy)
grid on

figure
subplot(2,1,1)
semilogy(freq_xy, abs(Gxy))
grid on
subplot(2,1,2)
plot(freq_xy, angle(Gxy))
% Low range, not loglog, we can change it

% Estimators - idk which would be better

H1 = Gxy./Gxx;
H2 = Gyy./Gyx;

figure
subplot(2,2,1)
loglog(freq_xx, abs(H1))
grid on

subplot(2,2,3)
plot(freq_xx, angle(H1))
grid on

subplot(2,2,2)
loglog(freq_xx, abs(H2))
grid on

subplot(2,2,4)
plot(freq_xx, angle(H2))

figure
loglog(freq_xx, abs(H1))
hold on
loglog(freq_xx, abs(H2))
yyaxis("left")
semilogx(freq_xx, coherence)
% From these plots we should be able to choose the correct estimator
%% Modal extraction










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

% As we are dealing with only first 4 modes of a simple beam the HP of well 
% separeted peaks is likely respected -> Sdof method shuld be reliable

%Gjk(omega) = Ajk / (-omega^2+ 1i*2*csi + omega_i^2) + RLjk/omega^2 + RHjk
% x = [om_i, csi_i, A_i, RL_i, RH_i];

G_num = @(omega,x) x(3) ./ (-omega.^2+ 1i*2*x(2)*omega*x(1) + x(1)^2) + x(4)./omega.^2 + x(5);

G_exp = H2; 
% or G_exp = H1 based on which one is the best esimatr

[~, all_peaks] = findpeaks(abs(G_exp),'MinPeakProminence', 0.2);

n_modes = 4; 
res_locs = all_peaks(1:n_modes);

om_nat = freq_xx(res_locs) * 2*pi;
dphase = diff(unwrap(angle(G_exp))) ./ diff(freq_xx*2*pi)';
f_max = freq_xx(end);
f_res = f_max/length(freq_xx);

x = zeros(5,n_modes);

for ii = 1:n_modes

    % interval
    range = 0.5; % hz
    range_idx(:,ii)= res_locs(ii) - round(range/f_res) : res_locs(ii) + round(range*f_res);
    fs_int(:,ii) = freq_xx(range_idx(:,ii));
    om_int(:,ii) = fs_int(:,ii) * 2*pi;

    %initial guesses
    x0(1) = om_nat(ii); %om_i
    x0(2)= -1 ./ (dphase(res_locs(ii)) * x0(1)); %csi0_i;
    x0(3) = G_exp(res_locs(ii))*1i*x0(2)*x0(1).^2*2; %A_i
    x0(4) = 0; %RL_i
    x0(5) = 0; %RH_i

    % error function
    err = @(x) [
        real(G_exp(range_idx(:,ii),1) - G_num(om_int(:,ii),x));  
        imag(G_exp(range_idx(:,ii),1) - G_num(om_int(:,ii),x)); 
        ];

    x(:,ii) = lsqnonlin(err,x0,[],[],[]);

end

%% graphical check
for ii = 1 : n_modes
    %amplitudes
    %experim FRF 
    subplot(2,n_modes,ii)
    semilogy(fs_int(:,ii),abs(G_exp(range_idx(:,ii))))
    hold on
    %extraolated numerical FRF
    semilogy(fs_int(:,ii),abs(G_num(om_int(:,ii),x(:,ii))),'or');
    grid on
    tit = sprintf('Modo %d - f = %.2f Hz', ii, x(1,ii)/2/pi);
    title(tit)

    %phases
    %experim FRF 
    subplot(2,n_modes,ii+n_modes)
    semilogy(fs_int(:,ii),angle(G_exp(range_idx(:,ii))))
    hold on
    %extraolated numerical FRF
    semilogy(fs_int(:,ii),angle(G_num(om_int(:,ii),x(:,ii))),'or');
    grid on

end 










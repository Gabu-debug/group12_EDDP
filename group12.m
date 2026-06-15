% Uploading data - .mat file?
data = load('modal_extraction.mat');

% Assembling variables
s = data.sensibility; 
raw_data = data.measures/s;  % Displacements, by a vibrometer i think

% Do we have 1 measuring point ??
%  we can repeat measurements more than 1 time

min_d = 0; % Looking at the data we chose a threshold, to discard
% the initial time of nothing

fs = data.fs;
dt = 1/fs;

for ii=1:length(raw_data(:,1)) % ! if they are column is 1,:
    index = find(d(ii,:) > min_d,'1','first');
    d(ii).data = d.raw_data(ii, index:end); % each data can have different 
    % dimension
    d(ii).t = 0:dt:dt*length(d(ii).data);
end

% We now have n repitition of the random signal we gave
% Problem: if the signal is random, we can't average to take out the noise

%% Solving for Circuit Component Values - Targeting Mode 1
% Values supplied to us:
R2 = 10E3; % Ohms
Ccap = 40E-9; % Farad, C,hat in NC circuit diagram.
Rcap = 0.9E6; % Ohms, R,hat in NC circuit diagram.
Beta = 0.7;
Cp1 = 40E-9; % Farad, mode 1 equivalent capacitance - we may get an updated value once we arrive.
Req = -50E6; % Ohms, equivalent resistance of Rs and R~ in parallel.

% Need R1, Rsh, and Rs.
C2 = Cp1/Beta;  % Capacitance of NC circuit. (Eq. 37)
R1 = (R2*Ccap)/C2;  % R1 from NC circuit. (Eq. 54)
Rtilda = R1*Rcap/R2;  % Ohms, R~ (Eq. 55)
Rs = Rtilda*Req/(Req + Rtilda);  % Ohms, balancing resistor. (Eq. 59 rearranged)

omega1 = "?";  % Open Circuit Nat. Freq. Mode 1: Needs to be solved for from modal analysis.
omegaHat1 = "?";  % Closed Circuit Nat. Freq. Mode 1: Needs to be solved for from modal analysis.
k1 = sqrt( omegaHat1^2/omega1^2 - 1 );  % Mode 1 coupling factor.
Rsh = (1-Beta2)/(omega1*Cp1) * sqrt(2/( k1^2/(1-Beta) + 2 -2*k1*sqrt(Beta/(1-Beta) ) ) );

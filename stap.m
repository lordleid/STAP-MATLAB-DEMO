%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:             Fully adaptive space-time processing
%                      for the side-looking case
% Description:         single target and two targets cases
% References:   
% [1] James Ward, "Space-Time Adaptive Processing for Airborne
%     Radar". MIT Lincoln Lab tech report 1015, 1994.
% [2] J. R. Guerci, "Space-Time Adaptive Processing for Radar".
%     Artech House, 2003.
% Conclusion:
% $Revision:	     $1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear;
lambda = 0.03; % wavelength, not neccessary in array signal processing
d = lambda/2; % interelement spacing

N = 10; % the number of elements in the ULA
M = 12; % the number of pulses per CPI

CNR = 50;   % dB Clutter to Noise Ratio  
SNR = 0;    % dB Signal to Noise Ratio 
JNR = 30;   % dB Jammer to Noise Ratio

noisePower = 1;
clutterPower = noisePower * 10^(CNR/10);
tgtPower = noisePower * 10^(SNR/10);
jammerPower = noisePower * 10^(JNR/10);

beta = 1; % Eq.(3.2) in [1],  control the way the clutter filling the angle Doppler plane


%% (1) Rc, clutter covariance
% ---the ground clutter ridge, see Eq.(3.2) in [1] ----
% the clutter Doppler shift induced on the clutter patch on the iso-range ring situated at the angle 'theta'
% fd_ClutterNormalized = beta*spatialFreq_normalized, where
% fd_ClutterNormalized = fd_Clutter/PRF;  % fd_Clutter is the Doppler frequency of clutter
% spatialFreq_normalized = d*sin(theta)/lambda;

% (debug1)--- for debug ---
No = 250;       % k-th clutter bins
sintheta = linspace(-1, 1, No);
clutterSpatialFreq_normalized = d./lambda*sintheta;
%% (debug1)--- for practice ---
%clutterAzimuth = -90:1:90; % since backlobe rejection
%clutterSpatialFreq_normalized = d./lambda*sind(clutterAzimuth); % d./lambda = 0.5

fd_ClutterNormalized = beta*clutterSpatialFreq_normalized; % fd_ClutterNormalized
Rc = complex(zeros(M*N));
V = zeros(M*N, length(clutterSpatialFreq_normalized));
for k = 1:length(clutterSpatialFreq_normalized)
	a_clutter = exp(-1j*2*pi*clutterSpatialFreq_normalized(k)*[0: N - 1].'); % spatial steering vector
	b_clutter = exp(-1j*2*pi*fd_ClutterNormalized(k)*[0:M - 1].');
	v_clutter = sqrt(clutterPower)*kron(b_clutter, a_clutter);	% space-time steering vector
	V(:, k) = v_clutter;
	Rc = Rc + v_clutter*v_clutter';
end
Rc = Rc./length(clutterSpatialFreq_normalized);


%% (2) Rn, noise covariance
% Noise signals decorrelate from pulse-to-pulse
% With this assumption, noise covariance matrix is
Rn = noisePower*eye(M*N);


%% (3) Rt, target covariance
tgtAzimuth = 0; 
tgtSpatialFreq_normalized = d./lambda*sind(tgtAzimuth);
fd_tgtNormalized = -0.1;  % the domain of fd_tgtNormalized is [-0.5, 0.5]
a_tgt = exp(-1j*2*pi*tgtSpatialFreq_normalized*[0: N - 1].'); % spatial steering vector
b_tgt = exp(-1j*2*pi*fd_tgtNormalized*[0:M - 1].');
v_tgt = sqrt(tgtPower)*kron(b_tgt, a_tgt);	% space-time steering vector
Rt = v_tgt*v_tgt';


%% (4) Rj, jammer covariance
jammerAzimuth = -30;
jammerSpatialFreq_normalized = d./lambda*sind(jammerAzimuth);
%fd_jammerNormalized = beta*jammerSpatialFreq_normalized; % debug
fd_jammerNormalized = 0.1; % fd_ClutterNormalized
a_jammer = exp(-1j*2*pi*jammerSpatialFreq_normalized*[0: N - 1].'); % spatial steering vector
b_jammer = exp(-1j*2*pi*fd_jammerNormalized*[0:M - 1].');
v_jammer = sqrt(jammerPower)*kron(b_jammer, a_jammer);	% space-time steering vector
Rj = v_jammer*v_jammer';


%% (5) R, the total space-time interference matrix due to clutter, jammer, noise
R = Rc + Rj + Rn;
[U, S, V] = svd(R);
figure; plot(10*log10(diag(S))); xlabel('Eigenvalue index'); ylabel('Eigenvalues, dB');


%% (6) space-time DBF spectrum which refers to as power spectrum density, Eq.(3.16)
% debug:Rc1 = Rc; Rj1 = Rj; Rn1 = Rn; Rt11 = Rt; save data Rc1 Rj1 Rn1 Rt11;
R_total = Rc + Rj + Rn + Rt;
inv4R_total = inv(R_total);

% ---> note(1): the following grid will not work well!
%azimuthGrid = linspace(-90, 90);%-90:1:90; % since backlobe rejection
%SpatialFreqGrid_normalized = d./lambda*sind(azimuthGrid); % d./lambda = 0.5
% <--- note(1): do a good job this way!
sintheta = linspace(-1, 1);
SpatialFreqGrid_normalized = d./lambda*sintheta; % d./lambda = 0.5

fdGrid_Normalized = beta*SpatialFreqGrid_normalized; % fd_Normalized
a_Grid = exp(-1j*2*pi*[0: N - 1].'*SpatialFreqGrid_normalized);
b_Grid = exp(-1j*2*pi*[0:M - 1].'*fdGrid_Normalized);
P_dbf = complex(zeros(length(fdGrid_Normalized), length(SpatialFreqGrid_normalized))); % DBF spectrum of the total return from the beamformer
P_capon = complex(zeros(length(fdGrid_Normalized), length(SpatialFreqGrid_normalized))); % Capon spectrum of the total return from the beamformer
for iSpatialFreq = 1:length(SpatialFreqGrid_normalized)
	for j_fdNormalized = 1:length(fdGrid_Normalized)
		v = kron(b_Grid(:, j_fdNormalized), a_Grid(:, iSpatialFreq));
		P_dbf(j_fdNormalized, iSpatialFreq) = v'*R_total*v; % Eq.(3.16) in [1]
		P_capon(j_fdNormalized, iSpatialFreq) = 1./(v'*inv4R_total*v); % Eq.(3.18) in [1]
	end
end
% Display the total return DBF spectrum
figure; 
imagesc(SpatialFreqGrid_normalized, fdGrid_Normalized, 10*log10(abs(P_dbf)))
set(gca,'ydir','normal');colorbar; xlabel('normalized spatial frequency'); ylabel('normalized Doppler');
title('Total DBF spectrum before STAP Detection of target, clutter, noise & jammer');
figure;
surf(SpatialFreqGrid_normalized, fdGrid_Normalized, 10*log10(abs(P_dbf)))
shading interp;colorbar; xlabel('normalized spatial frequency'); ylabel('normalized Doppler');
title('Total DBF spectrum before STAP Detection of target, clutter, noise & jammer');

% Display the total Capon spectrum, MVDR spctrum
figure; 
imagesc(SpatialFreqGrid_normalized, fdGrid_Normalized, 10*log10(abs(P_capon)))
set(gca,'ydir','normal');colorbar; xlabel('normalized spatial frequency'); ylabel('normalized Doppler');
title('Total Capon spectrum before STAP Detection of target, clutter, noise & jammer');
figure;
surf(SpatialFreqGrid_normalized, fdGrid_Normalized, 10*log10(abs(P_capon)))
shading interp;colorbar; xlabel('normalized spatial frequency'); ylabel('normalized Doppler');
title('Total Capon spectrum before STAP Detection of target, clutter, noise & jammer');


%% (7) Calculate optimal weights
inv_R = inv(R); % R = Rc + Rj + Rn;
wopt = inv_R * v_tgt;
% wopt = inv_R * v_tgt/(v_tgt'*inv_R * v_tgt); % Eq.(3.43) in [1]

%% (8) the optimal space-time output response
sintheta = linspace(-1, 1);
SpatialFreqGrid_normalized = d./lambda*sintheta; % d./lambda = 0.5
fdGrid_Normalized = beta*SpatialFreqGrid_normalized; % fd_Normalized
a_Grid = exp(-1j*2*pi*[0: N - 1].'*SpatialFreqGrid_normalized);
b_Grid = exp(-1j*2*pi*[0:M - 1].'*fdGrid_Normalized);
Y = complex(zeros(length(fdGrid_Normalized), length(SpatialFreqGrid_normalized))); % Capon spectrum of the total return from the beamformer
for iSpatialFreq = 1:length(SpatialFreqGrid_normalized)
	for j_fdNormalized = 1:length(fdGrid_Normalized)
		v = kron(b_Grid(:, j_fdNormalized), a_Grid(:, iSpatialFreq));
		Y(j_fdNormalized, iSpatialFreq) = wopt'*v; % 
	end
end
% Display the total output power spectrum, DBF spctrum
figure; 
imagesc(SpatialFreqGrid_normalized, fdGrid_Normalized, 10*log10(abs(Y).^2))
set(gca,'ydir','normal');colorbar; xlabel('normalized spatial frequency'); ylabel('normalized Doppler');
title('STAP Detection of target & jammer; clutter removed');
figure; 
mesh(SpatialFreqGrid_normalized, fdGrid_Normalized, 10*log10(abs(Y).^2))
colorbar; xlabel('normalized spatial frequency'); ylabel('normalized Doppler');
title('STAP Detection of target & jammer; clutter removed');

%% (9) The SINR loss defined in Eq.(120) in [2]
sintheta = linspace(-1, 1, 181);
SpatialFreqGrid_normalized = d./lambda*sintheta; % d./lambda = 0.5
fdGrid_Normalized = beta*SpatialFreqGrid_normalized; % fd_Normalized
a_Grid = exp(-1j*2*pi*[0: N - 1].'*SpatialFreqGrid_normalized);
b_Grid = exp(-1j*2*pi*[0:M - 1].'*fdGrid_Normalized);
SINR_loss = zeros(length(fdGrid_Normalized), length(SpatialFreqGrid_normalized)); % Capon spectrum of the total return from the beamformer
for iSpatialFreq = 1:length(SpatialFreqGrid_normalized)
	for j_fdNormalized = 1:length(fdGrid_Normalized)
		v = kron(b_Grid(:, j_fdNormalized), a_Grid(:, iSpatialFreq));
		SINR_current = v'*inv_R*v; 
		SNRopt = v'*v; % the optimal SNR, Eq.(119) in [2]
		SINR_loss(j_fdNormalized, iSpatialFreq) = SINR_current/SNRopt; % 
	end
end

figure; plot(SpatialFreqGrid_normalized, 10*log10(abs(SINR_loss(:, 91))));
xlabel('Normalized Doppler'); ylabel('SINR loss, dB');


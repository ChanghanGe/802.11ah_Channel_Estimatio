clear 
clc

% waveform configuration
cfgS1G = wlanS1GConfig;
cfgS1G.ChannelBandwidth = 'CBW2'; % 2 MHz channel bandwidth
cfgS1G.Preamble = 'Long';        % Short preamble
cfgS1G.NumTransmitAntennas = 1;   % 2 transmit antennas
cfgS1G.NumSpaceTimeStreams = 1;   % 2 space-time streams
cfgS1G.APEPLength = 256;          % APEP length in bytes
cfgS1G.MCS = 0;                   % 64-QAM rate-5/6

% STF sequence
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f , 64);
sts_t = sts_t(1:16);

% Indices for accessing each field within the time-domain packet
fieldInd = wlanFieldIndices(cfgS1G);

%ofdm information
ofdmInfo = wlanS1GOFDMInfo('S1G-Data',cfgS1G);

symbolLength = 16; % In samples
lenLSTF = symbolLength*10; % Length of 10 L-STF symbols 
 
%symbol transmitted
load('/Users/changhange/Desktop/samples_bpsk.mat');
tx_signal = waveStruct.waveform;
clear waveStruct

% Read binary file and decode it to complex double
rx_signal = read_complex_binary('./receivedaod.bin',100000000000000000000);

% cross = xcorr(rx_signal,sts_t);
% cross = cross(length(rx_signal):end); 
% ind=find (abs(cross)>0.9*max(abs(cross)));

cross = xcorr(rx_signal,sts_t);
cross = cross(length(rx_signal):end); 
ind=find(abs(cross)>0.6*max(abs(cross)));

%Packet Detection and Synchronization
% initial_offset = 0;
% [startOffset,M] = wlanPacketDetect(rx_signal, 'CBW5', initial_offset, 0.99);
% offset = initial_offset + startOffset;

%First Batch
rxSync_batch = zeros(700,7120,4);
for i = 1:700
    rxSync_batch(i,:,:) = reshape(rx_signal(ind(1)+(i-1)*28480+1:ind(1)+i*28480,:),7120,4);
end

%Demodulation and Channel Estimation
% First Antenna
rx_ant1_batch_chanEst = zeros(700,56,1);
rx_ant2_batch_chanEst = zeros(700,56,1);
rx_ant3_batch_chanEst = zeros(700,56,1);
rx_ant4_batch_chanEst = zeros(700,56,1);
for i = 1:700
    rx_ant1_batch_chanEst(i,:,:) = ChannelEst_Demod_11ah(rxSync_batch(i,:,1)', fieldInd, cfgS1G, ofdmInfo);
    rx_ant2_batch_chanEst(i,:,:) = ChannelEst_Demod_11ah(rxSync_batch(i,:,2)', fieldInd, cfgS1G, ofdmInfo);
    rx_ant3_batch_chanEst(i,:,:) = ChannelEst_Demod_11ah(rxSync_batch(i,:,3)', fieldInd, cfgS1G, ofdmInfo);
    rx_ant4_batch_chanEst(i,:,:) = ChannelEst_Demod_11ah(rxSync_batch(i,:,4)', fieldInd, cfgS1G, ofdmInfo);
end
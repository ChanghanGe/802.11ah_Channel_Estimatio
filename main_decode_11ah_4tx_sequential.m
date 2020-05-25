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
num_packet = 400;

%ofdm information
ofdmInfo = wlanS1GOFDMInfo('S1G-Data',cfgS1G);

symbolLength = 16; % In samples
lenLSTF = symbolLength*10; % Length of 10 L-STF symbols 
 
%symbol transmitted
%load('samples_bpsk.mat');
%tx_signal = waveStruct.waveform;
%clear waveStruct
tx_signal = read_complex_binary('/Users/changhange/Downloads/antenna/antenna4_051920_bin.mat',1000000000000000000000);
[~, txPSDU] = ChannelEst_Demod_11ah(tx_signal(21361:end,:), fieldInd, cfgS1G, ofdmInfo);

% Read binary file and decode it to complex double
rx_signal = read_complex_binary('./receivedaod/receivedaod.bin.aa',1000000000000000000000);

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
rxSync_batch = complex(zeros(num_packet,7120,4));
for i = 1:num_packet
    rxSync_batch(i,:,:) = reshape(rx_signal(ind(1)+(i-1)*28480+1:ind(1)+i*28480,:),7120,4);
end

%Demodulation and Channel Estimation
% First Antenna
rx_batch_chanEst = complex(zeros(num_packet,56,1));
rx_batch_PSDU = int64(zeros(num_packet,2064,4));
for i = 1:num_packet
    for j = 1:4
        [rx_batch_chanEst(i,:,j), rx_batch_PSDU(i,:,j)] = ChannelEst_Demod_11ah(rxSync_batch(i,:,j)', fieldInd, cfgS1G, ofdmInfo);
    end
end
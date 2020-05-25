clear 
clc

%% 802.11ah packet structure
% waveform configuration
cfgS1G = wlanS1GConfig;
cfgS1G.ChannelBandwidth = 'CBW2'; % 2 MHz channel bandwidth
cfgS1G.Preamble = 'Long';        % Long preamble
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
 
%% signal transmitted
tx_signal = read_complex_binary('/Users/changhange/Desktop/samples_80211ah_2MHz_NHT.bin');
[~, txPSDU] = ChannelEst_Demod_11ah(tx_signal, fieldInd, cfgS1G, ofdmInfo);
txPSDU = int64(txPSDU);

%% signal received
% Read binary file and decode it to complex double
rx_signal_1 = read_complex_binary('/Users/changhange/Downloads/data_0521/antenna_wo_external_lo/antenna1.bin',10000000000000000);
rx_signal_2 = read_complex_binary('/Users/changhange/Downloads/data_0521/antenna_wo_external_lo/antenna2.bin',10000000000000000);
rx_signal_3 = read_complex_binary('/Users/changhange/Downloads/data_0521/antenna_wo_external_lo/antenna3.bin',10000000000000000);
rx_signal_4 = read_complex_binary('/Users/changhange/Downloads/data_0521/antenna_wo_external_lo/antenna4.bin',10000000000000000);

num_packet = 1000;
num_ant = 4;

%% Packet Detection and Synchronization
% ind_1 = wlanPacketDetect_80211ah(rx_signal_1,cfgS1G.ChannelBandwidth,0,0.5);
% ind_2 = wlanPacketDetect_80211ah(rx_signal_2,cfgS1G.ChannelBandwidth,0,0.5);
% ind_3 = wlanPacketDetect_80211ah(rx_signal_3,cfgS1G.ChannelBandwidth,0,0.5);
% ind_4 = wlanPacketDetect_80211ah(rx_signal_4,cfgS1G.ChannelBandwidth,0,0.5);

cross_1 = xcorr(rx_signal_1,tx_signal);
cross_1 = cross_1(length(rx_signal_1):end); 
ind_1=find(abs(cross_1)>0.8*max(abs(cross_1)));

cross_2 = xcorr(rx_signal_2,tx_signal);
cross_2 = cross_2(length(rx_signal_2):end); 
ind_2=find(abs(cross_2)>0.8*max(abs(cross_2)));

cross_3 = xcorr(rx_signal_3,tx_signal);
cross_3 = cross_3(length(rx_signal_3):end); 
ind_3=find(abs(cross_3)>0.8*max(abs(cross_3)));

cross_4 = xcorr(rx_signal_4,tx_signal);
cross_4 = cross_4(length(rx_signal_4):end); 
ind_4=find(abs(cross_4)>0.8*max(abs(cross_4)));

rxSync_batch = complex(zeros(num_packet,7120,4));
for i = 1:num_packet
    rxSync_batch(i,:,1) = rx_signal_1(ind_1(1)+(i-1)*7120+1:ind_1(1)+i*7120,:)';
    rxSync_batch(i,:,2) = rx_signal_2(ind_2(1)+(i-1)*7120+1:ind_2(1)+i*7120,:)';
    rxSync_batch(i,:,3) = rx_signal_1(ind_3(1)+(i-1)*7120+1:ind_3(1)+i*7120,:)';
    rxSync_batch(i,:,4) = rx_signal_2(ind_4(1)+(i-1)*7120+1:ind_4(1)+i*7120,:)';
end

%% Demodulation and Channel Estimation
rx_batch_chanEst = complex(zeros(num_packet,56,4));
rx_batch_PSDU = int64(zeros(num_packet,2064,4));
for i = 1:num_packet
    for j = 1:num_ant
        [rx_batch_chanEst(i,:,j), rx_batch_PSDU(i,:,j)] = ChannelEst_Demod_11ah(rxSync_batch(i,:,j)', fieldInd, cfgS1G, ofdmInfo);
    end
end

%% Evaluate Packet Error Rate (PER)
PE = zeros(1,num_ant);
for i = 1:num_packet
    for j = 1:num_ant
        packetError = any(biterr(txPSDU,rx_batch_PSDU(i,:,j)'));
        PE(1,j) = PE(1,j)+packetError;
    end
end

PER = PE./(ones(1,num_ant)*num_packet);
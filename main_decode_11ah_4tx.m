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
rx_signal = read_complex_binary('./receivedaod/receivedaod.bin.ah',1000000000);

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

rxSync_ant1 = rx_signal(ind(1)+1:ind(1)+7120,:);
rxSync_ant2 = rx_signal(ind(1)+7121:ind(1)+7120*2,:);
rxSync_ant3 = rx_signal(ind(1)+7120*2+1:ind(1)+7120*3,:);
rxSync_ant4 = rx_signal(ind(1)+7120*3+1:ind(1)+7120*4,:);

%Demodulation      
% LTF demodulation and channel estimation
% Demodulate S1G-LTF1
txLTF1 = tx_signal(fieldInd.S1GLTF1(1):fieldInd.S1GLTF1(2),:);
demodLTF1_tx = wlanS1GDemodulate(txLTF1,'S1G-LTF1',cfgS1G);

rxLTF1_ant1 = rxSync_ant1(fieldInd.S1GLTF1(1):fieldInd.S1GLTF1(2),:);
demodLTF1_ant1 = wlanS1GDemodulate(rxLTF1_ant1,'S1G-LTF1',cfgS1G);

rxLTF1_ant2 = rxSync_ant2(fieldInd.S1GLTF1(1):fieldInd.S1GLTF1(2),:);
demodLTF1_ant2 = wlanS1GDemodulate(rxLTF1_ant2,'S1G-LTF1',cfgS1G);

rxLTF1_ant3 = rxSync_ant3(fieldInd.S1GLTF1(1):fieldInd.S1GLTF1(2),:);
demodLTF1_ant3 = wlanS1GDemodulate(rxLTF1_ant3,'S1G-LTF1',cfgS1G);

rxLTF1_ant4 = rxSync_ant4(fieldInd.S1GLTF1(1):fieldInd.S1GLTF1(2),:);
demodLTF1_ant4 = wlanS1GDemodulate(rxLTF1_ant4,'S1G-LTF1',cfgS1G);

%Channel Estimation
if cfgS1G.NumSpaceTimeStreams>1
    % Use S1G-LTF1 and S1G-LTF2N for channel estimation
    rxLTF2N = rxSync(fieldInd.S1GLTF2N(1):fieldInd.S1GLTF2N(2),:);
    demodLTF2N = wlanS1GDemodulate(rxLTF2N,'S1G-LTF2N',cfgS1G);
    chanEst = s1gLTFChannelEstimate([demodLTF1 demodLTF2N],cfgS1G);
else
    % Use only S1G-LTF1 for channel estimation
    chanEst_tx = s1gLTFChannelEstimate(demodLTF1_tx,cfgS1G);
    chanEst_ant1 = s1gLTFChannelEstimate(demodLTF1_ant1,cfgS1G);
    chanEst_ant2 = s1gLTFChannelEstimate(demodLTF1_ant2,cfgS1G);
    chanEst_ant3 = s1gLTFChannelEstimate(demodLTF1_ant3,cfgS1G);
    chanEst_ant4 = s1gLTFChannelEstimate(demodLTF1_ant4,cfgS1G);
end

% Estimate Noise Level
noiseVarEst_tx = helperNoiseEstimate(demodLTF1_tx);
noiseVarEst_ant1 = helperNoiseEstimate(demodLTF1_ant1);
noiseVarEst_ant2 = helperNoiseEstimate(demodLTF1_ant2);
noiseVarEst_ant3 = helperNoiseEstimate(demodLTF1_ant3);
noiseVarEst_ant4 = helperNoiseEstimate(demodLTF1_ant4);

% Extract S1G-Data field
rxData_tx = tx_signal(fieldInd.S1GData(1):fieldInd.S1GData(2),:);
rxData_ant1 = rxSync_ant1(fieldInd.S1GData(1):fieldInd.S1GData(2),:);
rxData_ant2 = rxSync_ant2(fieldInd.S1GData(1):fieldInd.S1GData(2),:);
rxData_ant3 = rxSync_ant3(fieldInd.S1GData(1):fieldInd.S1GData(2),:);
rxData_ant4 = rxSync_ant4(fieldInd.S1GData(1):fieldInd.S1GData(2),:);

% OFDM demodulation
demodSym_tx = wlanS1GDemodulate(rxData_tx,'S1G-Data',cfgS1G);
demodSym_ant1 = wlanS1GDemodulate(rxData_ant1,'S1G-Data',cfgS1G);
demodSym_ant2 = wlanS1GDemodulate(rxData_ant2,'S1G-Data',cfgS1G);
demodSym_ant3 = wlanS1GDemodulate(rxData_ant3,'S1G-Data',cfgS1G);
demodSym_ant4 = wlanS1GDemodulate(rxData_ant4,'S1G-Data',cfgS1G);

% Extract data subcarriers from demodulated symbols and channel
% estimate
demodDataSym_tx = demodSym_tx(ofdmInfo.DataIndices,:,:);
demodDataSym_ant1 = demodSym_ant1(ofdmInfo.DataIndices,:,:);
demodDataSym_ant2 = demodSym_ant2(ofdmInfo.DataIndices,:,:);
demodDataSym_ant3 = demodSym_ant3(ofdmInfo.DataIndices,:,:);
demodDataSym_ant4 = demodSym_ant4(ofdmInfo.DataIndices,:,:);

chanEstData_tx = chanEst_tx(ofdmInfo.DataIndices,:,:);
chanEstData_ant1 = chanEst_ant1(ofdmInfo.DataIndices,:,:);
chanEstData_ant2 = chanEst_ant2(ofdmInfo.DataIndices,:,:);
chanEstData_ant3 = chanEst_ant3(ofdmInfo.DataIndices,:,:);
chanEstData_ant4 = chanEst_ant4(ofdmInfo.DataIndices,:,:);

% MMSE frequency domain equalization
[eqDataSym_tx,csi_tx] = helperSymbolEqualize(demodDataSym_tx,chanEstData_tx,noiseVarEst_tx); 
[eqDataSym_ant1,csi_ant1] = helperSymbolEqualize(demodDataSym_ant1,chanEstData_ant1,noiseVarEst_ant1); 
[eqDataSym_ant2,csi_ant2] = helperSymbolEqualize(demodDataSym_ant2,chanEstData_ant2,noiseVarEst_ant2); 
[eqDataSym_ant3,csi_ant3] = helperSymbolEqualize(demodDataSym_ant3,chanEstData_ant3,noiseVarEst_ant3); 
[eqDataSym_ant4,csi_ant4] = helperSymbolEqualize(demodDataSym_ant4,chanEstData_ant4,noiseVarEst_ant4); 

% Recover PSDU bits
txPSDU = s1gDataBitRecover(eqDataSym_tx,noiseVarEst_tx,csi_tx,cfgS1G);
rxPSDU_ant1 = s1gDataBitRecover(eqDataSym_ant1,noiseVarEst_ant1,csi_ant1,cfgS1G);
rxPSDU_ant2 = s1gDataBitRecover(eqDataSym_ant2,noiseVarEst_ant2,csi_ant2,cfgS1G);
rxPSDU_ant3 = s1gDataBitRecover(eqDataSym_ant3,noiseVarEst_ant3,csi_ant3,cfgS1G);
rxPSDU_ant4 = s1gDataBitRecover(eqDataSym_ant4,noiseVarEst_ant4,csi_ant4,cfgS1G);
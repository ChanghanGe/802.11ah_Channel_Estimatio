% Generated by MATLAB(R) 9.8 (R2020a) and WLAN Toolbox 3.0 (R2020a).
% Generated on: 25-May-2020 00:59:50
clear
clc

%% Generate 802.11ah Waveform
% 802.11ah configuration:
s1gGCfg = wlanS1GConfig('ChannelBandwidth', 'CBW2', ...
    'Preamble', 'Short', ...
    'NumUsers', 1, ...
    'NumTransmitAntennas', 1, ...
    'NumSpaceTimeStreams', [1], ...
    'SpatialMapping', 'Direct', ...
    'STBC', false, ...
    'MCS', 0, ...
    'APEPLength', 256, ...
    'GuardInterval', 'Long', ...
    'PartialAID', 37, ...
    'UplinkIndication', false, ...
    'Color', 0, ...
    'TravelingPilots', false, ...
    'ResponseIndication', 'None', ...
    'RecommendSmoothing', true);

numPackets = 1;
% input bit source:
in = ones(2064, 1);


% waveform generation:
waveform = wlanWaveformGenerator(in, s1gGCfg, ...
    'NumPackets', numPackets, ...
    'IdleTime', 0, ...
    'ScramblerInitialization', 93, ...
    'WindowTransitionTime', 1e-07);

Fs = wlanSampleRate(s1gGCfg); 								 % sample rate of waveform

%% Visualize 802.11ah Waveform
% Spectrum Analyzer
spectrum = dsp.SpectrumAnalyzer('SampleRate', Fs);
spectrum(waveform);
release(spectrum);

write_complex_binary(waveform,'./txdata/samples_80211ah_2MHz_NHT_Short_preamble_0525.bin');



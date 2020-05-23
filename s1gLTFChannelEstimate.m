function est = s1gLTFChannelEstimate(rxSym,cfgS1G,varargin)
%s1gLTFChannelEstimate Channel estimation using S1G-LTF1 and S1G-LTF2N 
%   EST = s1gLTFChannelEstimate(RXSYM,CFGS1G) returns the estimated channel
%   between all space-time streams and receive antennas using S1G-LTF1 and
%   S1G-LTF2N for a S1G short preamble packet. The channel estimate
%   includes the effect of the applied spatial mapping matrix and cyclic
%   shifts at the transmitter.
%
%   Estimation can be performed using only the S1G-LTF1, or using both the
%   S1G-LTF1 and S1G-LTF2N fields. To estimate the channel for all
%   space-time streams, both fields are required.
%
%   EST is an array characterizing the estimated channel for the data and
%   pilot subcarriers. The size depends on whether S1G-LTF1 and S1G-SLTF2N
%   are used for estimation, or only S1G-LTF1. When both fields are used,
%   EST is a complex Nst-by-Nsts-by-Nr array characterizing the estimated
%   channel for the data and pilot subcarriers, where Nst is the number of
%   occupied subcarriers, Nsts is the total number of space-time streams,
%   and Nr is the number of receive antennas. If channel estimation is
%   performed using only the S1G-LTF1 field, EST is a complex
%   Nst-by-1-by-Nr array. The singleton dimension corresponds to the single
%   transmitted stream in S1G-LTF1, which includes the combined cyclic
%   shifts if multiple transmit antennas are used.
%
%   RXSYM is a complex Nst-by-Nsym-by-Nr array containing demodulated
%   concatenated S1G-LTF1 and S1G-LTF2N OFDM symbols. Nsym is the number of
%   demodulated S1G-LTF symbols. If RXSYM contains only the S1G-LTF1
%   demodulated field, channel estimation is performed using this field.
%
%   CFGS1G is a packet format configuration object of type <a href="matlab:help('wlanS1GConfig')">wlanS1GConfig</a>. 
%
%   EST = s1gLTFChannelEstimate(...,SPAN) performs frequency smoothing by
%   using a moving average filter across adjacent subcarriers to reduce the
%   noise on the channel estimate. The span of the filter in subcarriers,
%   SPAN, must be odd. If adjacent subcarriers are highly correlated
%   frequency smoothing will result in significant noise reduction, however
%   in a highly frequency selective channel smoothing may degrade the
%   quality of the channel estimate.
 
%   Copyright 2017-2018 The MathWorks, Inc.

%#codegen

narginchk(2,3);

% Validate the format configuration object
%validateattributes(cfgS1G,{'wlanS1GConfig'},{'scalar'},mfilename,'S1G format configuration object');
%coder.internal.errorIf(~strcmp(packetFormat(cfgS1G),'S1G-Short'),'wlan:s1gLTFChannelEstimate:OnlyShortPreambleSupported');
chanBW = cfgS1G.ChannelBandwidth;
cpType = cfgS1G.GuardInterval;

% Validate symbol type
validateattributes(rxSym,{'single','double'},{'3d'},mfilename,'S1G LTF OFDM symbol(s)');
[numST,numSym,numRx] = size(rxSym);
numSTS = sum(cfgS1G.NumSpaceTimeStreams);
if isempty(rxSym)
    est = zeros(numST,numSTS,numRx);
    return;
end

% Verify number of subcarriers to estimate
% Get OFDM configuration
[cfgOFDM,dataInd] = wlan.internal.s1gOFDMConfig(chanBW,cpType,'LTF',numSTS);
FFTLen = cfgOFDM.FFTLength;
[ind, sortedDataPilotIdx] = sort([cfgOFDM.DataIndices; cfgOFDM.PilotIndices]);
k = ind-FFTLen/2-1; % Active subcarrier frequency index
coder.internal.errorIf(numST~=numel(ind),'wlan:wlanChannelEstimate:IncorrectNumSC',numel(ind),numST);
% Get S1G-LTF sequences
% Nltf is the number of LTF fields (LTF1 + LTF2N)
numLTF = wlan.internal.numVHTLTFSymbols(numSTS);
coder.internal.errorIf((size(rxSym,2)<2)||(size(rxSym,2)<(2+numLTF-1)),'wlan:s1gLTFChannelEstimate:InvalidRxSymInputSize',numel(ind),(2+numLTF-1));

if nargin == 3
    enableFreqSmoothing = true;
    span = varargin{1};
else    
    enableFreqSmoothing = false;
end

if numSym==2
    % Only LTF1 channel estimation
    
    % Get S1G-LTF sequences
    % In SHORT preamble, Nltf is unknown and there is only 1 space-time stream
    % to estimate the channel from LTF1.
    ltf = wlan.internal.vhtltfSequence(cfgS1G.ChannelBandwidth, 1);

    % If one space time stream then use LS estimation directly
    ls = rxSym./repmat(ltf(ind,:,:), 1, size(rxSym,2) ,numRx); % Least-square estimate   
    est = mean(ls,2); % Average over the symbols
    
    % Perform frequency smoothing
    if enableFreqSmoothing
        est = frequencySmoothing(est,chanBW,span);
    end
else
    % LTF1 and LTF2N estimation

    % Verify enough symbols to estimate
    coder.internal.errorIf(numSym<numLTF,'wlan:wlanChannelEstimate:NotEnoughSymbols',numSTS,numLTF,numSym);

    % Get cyclic shifts applied at transmitter (they affect when numSTS>1)
    csh = wlan.internal.getCyclicShiftVal('S1G',numSTS,wlan.internal.cbwStr2Num(chanBW));

    % One symbol per LTF. Average 2 symbols of LTF1
    rxSymLTF1 = rxSym(:,1:2,:);
    rxSymLTF2N = rxSym(:,2+(1:(numLTF-1)),:); 
    rxSymLTF1 = mean(rxSymLTF1,2);
    symLTF = cat(2,rxSymLTF1,rxSymLTF2N);

    % Perform channel estimation for data carrying subcarriers as we
    % must interpolate the pilots
    estData = wlan.internal.vhtltfEstimate(symLTF(dataInd,:,:),chanBW,numSTS,cfgOFDM.DataIndices);

    % Undo cyclic shift for each STS before averaging and interpolation
    estData = wlan.internal.cyclicShiftChannelEstimate(estData,-csh,FFTLen,k(dataInd));

    % Estimate pilot subcarriers
    estPilots = pilotInterpolation(estData,FFTLen,cfgOFDM.DataIndices,cfgOFDM.PilotIndices);

    % Combine data and pilots into one container
    allSubcarriers = [estData; estPilots];
    est = allSubcarriers(sortedDataPilotIdx,:,:);
    
    % Perform frequency smoothing
    if enableFreqSmoothing
        est = frequencySmoothing(est,chanBW,span);
    end
    
    % Re-apply cyclic shift after averaging and interpolation
    if numSTS>1  
        est = wlan.internal.cyclicShiftChannelEstimate(est,csh,FFTLen,k);
    end
end

end

function chanEst = frequencySmoothing(chanEst,chanBW,span)
% Smooth segments between DC gaps
    switch chanBW
        case 'CBW2'
            numGroups = 1;
        case 'CBW4'
            numGroups = 2;
        case 'CBW8'
            numGroups = 2;
        otherwise % 'CBW16'
            numGroups = 4;
    end
    groupSize = size(chanEst,1)/numGroups;
    for i = 1:numGroups
        idx = (1:groupSize)+(i-1)*groupSize;
        chanEst(idx,:,:) = wlan.internal.frequencySmoothing(chanEst(idx,:,:),span);
    end
end

function estPilots = pilotInterpolation(estData,Nfft,dataIndices,pilotIndices)
% Interpolate over the pilot locations

    numSTS = size(estData,2);
    numRxAnts = size(estData,3);

    % Construct full FFT size to allow us to interpolate over DC nulls
    est = complex(ones(Nfft,numSTS,numRxAnts),ones(Nfft,numSTS,numRxAnts));
    est(dataIndices,:,:) = estData;
       
    % Interpolate over missing parts of the waveform in magnitude and
    % phase (as opposed to real and imaginary)
    magPart = interp1(dataIndices,abs(est(dataIndices,:,:)),1:Nfft);
    phasePart = interp1(dataIndices,unwrap(angle(est(dataIndices,:,:))),1:Nfft);
    [realPart,imagPart] = pol2cart(phasePart,magPart);
    estInterp = complex(realPart,imagPart);
    if isrow(estInterp)
        est = estInterp(:,:,1).';
    else
        est = estInterp;
    end
    
    % Extract pilots
    estPilots = est(pilotIndices,:,:);

end
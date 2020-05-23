function dataBits = s1gDataBitRecover(eqSym,noiseVar,csi,cfgS1G)
%s1gDataBitRecover Recover data bits from S1G Data field
%
%   DATABITS = s1gDataBitRecover(EQSYM,NOISEVAREST,CSI,CFGS1G) recovers the
%   data bits from the Data field of an S1G transmission. Only a
%   single-user transmission is supported.
%
%   DATABITS is an int8 column vector of length 8*CFGS1G.PSDULength
%   containing the recovered information bits.
%
%   EQSYM is the equalized S1G-Data field OFDM symbols, specified as a
%   Nsd-by-Nsym-by-Nss matrix of real or complex values. Nsd is the number
%   of data subcarriers in the S1G-Data field, Nsym is the number of OFDM
%   symbols, and Nss is the number of spatial-streams.
%
%   NOISEVAREST is the noise variance estimate, specified as a nonnegative
%   scalar.
%
%   CSI is a Nsd-by-Nss matrix containing the channel state information per
%   subcarrier and spatial stream.
%
%   CFGS1G is the format configuration object of type <a href="matlab:help('wlanS1GConfig')">wlanS1GConfig</a>, which 
%   specifies the parameters for the S1G format.

%   Copyright 2017-2018 The MathWorks, Inc.

%#codegen

% Validate configuration object
validateattributes(cfgS1G,{'wlanS1GConfig'},{'scalar'},mfilename,'cfgS1G');
cfgInfo = validateConfig(cfgS1G,'MCS');
coder.internal.errorIf(cfgS1G.NumUsers>1,'wlan:s1gDataBitRecover:MultiUserNotSupported'); % SU only

% Validate noise variance
validateattributes(noiseVar,{'double'},{'real','scalar','nonnegative','finite'},mfilename,'noiseVar');

% Rate parameters
chanBW = cfgS1G.ChannelBandwidth;
mcsTable = wlan.internal.getRateTable(cfgS1G);
% Set up some implicit configuration parameters
numBPSCS = mcsTable.NBPSCS(1);     % Number of coded bits per single carrier per spatial stream
numCBPS = mcsTable.NCBPS(1);       % Number of coded bits per OFDM symbol
rate = mcsTable.Rate(1);           % Coding rate
numES = mcsTable.NES(1);           % Number of encoded streams (BCC encoders)
numSS = mcsTable.Nss(1);           % Number of spatial streams
numSeg = strcmp(chanBW,'CBW16')+1; % Number of segments (BCC interleaver blocks): 2 for 'CBW16', 1 for the rest
% Number of coded bits per OFDM symbol, per spatial stream, per segment (BCC interleaver block)
numCBPSSI = numCBPS/numSS/numSeg;

% Validate data input (now we know numSS)
validateattributes(eqSym,{'double'},{'3d','finite'},mfilename,'eqSym'); 
% If input eqDataSym is empty then do not attempt to decode; return empty
ofdmInfo = wlanS1GOFDMInfo('S1G-Data',cfgS1G);
Nsd = size(ofdmInfo.DataIndices,1); % Number of data carrying subcarriers
coder.internal.errorIf((size(eqSym,1)~=Nsd)||(size(eqSym,3)~=numSS),'wlan:s1gDataBitRecover:InvalidDataSymInputSize',Nsd,cfgInfo.NumDataSymbols,numSS);
if isempty(eqSym) % Return no bits if no OFDM symbols passed
    dataBits = zeros(0,1,'int8');
    return;
end
coder.internal.errorIf(size(eqSym,2)<cfgInfo.NumDataSymbols,'wlan:s1gDataBitRecover:InvalidDataSymInputSize',Nsd,cfgInfo.NumDataSymbols,numSS);

% Validate CSI
validateattributes(csi,{'double'},{'2d','finite'},mfilename,'csi'); 
coder.internal.errorIf(any(size(csi)~=[Nsd numSS]),'wlan:s1gDataBitRecover:InvalidInputCSISize',Nsd,cfgS1G.NumSpaceTimeStreams(1));

% Segment parsing of symbols
parserOut = wlanSegmentParseSymbols(eqSym(:,1:cfgInfo.NumDataSymbols,:),chanBW); % [Nsd/Nseg Nsym Nss Nseg]
csiParserOut = wlanSegmentParseSymbols(reshape(csi,[],1,numSS),chanBW);          % [Nsd/Nseg 1 Nss Nseg]

% Constellation demapping
qamDemodOut = wlanConstellationDemap(parserOut,noiseVar,numBPSCS); % [Ncbpssi Nsym Nss Nseg]

% Apply CSI and concatenate OFDM symbols in the first dimension
csiAppOut = applyCSI(qamDemodOut,csiParserOut,numBPSCS); % [(Ncbpssi*Nsym),Nss,Nseg]

% BCC Deinterleaving
deintlvrOut = wlanBCCDeinterleave(csiAppOut,'VHT',numCBPSSI,chanBW); % [(Ncbpssi*Nsym),Nss,Nseg]

% Segment deparsing of bits
segDeparserOut = wlanSegmentDeparseBits(deintlvrOut,chanBW,numES,numCBPS,numBPSCS); % [(Ncbpss*Nsym),Nss]

% Stream deparsing
streamDeparserOut = wlanStreamDeparse(segDeparserOut(:,:),numES,numCBPS,numBPSCS); % [(Ncbps*Nsym/Nes),Nes]

% Repetition decoding for MCS 10
if cfgS1G.MCS(1) == 10
    combinedData = combineRepeatedMCS10(streamDeparserOut);
else
    combinedData = streamDeparserOut;
end

% Channel decoding for BCC
numTailBits = 6;
chanDecOut = wlanBCCDecode(combinedData,rate);
% Channel decoder deparser
preDescramBits = reshape(chanDecOut(1:end-numTailBits,:)',[],1);

% Derive initial state of the scrambler and descramble
scramInit = wlan.internal.scramblerInitialState(preDescramBits(1:7));
if ~all(scramInit==0) 
    descramBits = wlanScramble(preDescramBits(1:8+8*cfgS1G.PSDULength(1)),scramInit);
else % Skip scrambling as recovered seed is invalid (0)
    descramBits = preDescramBits(1:8+8*cfgS1G.PSDULength(1));
end

% Discard service bits
dataBits = descramBits(9:end);

end

function y = applyCSI(x,csi,numBPSCS)
    % Apply CSI and concatenate OFDM symbols in the first dimension
    [~,NumDataSymbols,numSS,numSeg] = size(x);
    tmp = bsxfun(@times, ...
        reshape(x,numBPSCS,[],NumDataSymbols,numSS), ...
        reshape(csi,1,[],1,numSS));
    y = reshape(tmp,[],numSS,numSeg);
end

function y = combineRepeatedMCS10(x)
    % Combine repeated bits for MCS 10
    numRep = 2;
    bpsPreRep = 12;
    reshapedX = reshape(x,bpsPreRep*numRep,[]);
    s = [1; 0; 0; 0; 0; 1; 0; 1; 0; 1; 1; 1]; % IEEE P802.11ah/D5.0 Eqn 24-45
    in = reshapedX(bpsPreRep+1:end,:);
    
    % Convert LLR into 0 and 1, before xor
    buffer = ones(size(in));
    buffer(in<0) = 0; % Negative LLRs are converted to zeros
    xorBuffer = bsxfun(@xor,buffer,s);
    xorin = in.*(-2*(buffer~=xorBuffer)+1); % Convert bits to +/-1 and multiply with input

    % Sum repeated LLRs
    y = reshape(reshapedX(1:bpsPreRep,:)+xorin,[],1);
end

function [startOffset,M] = wlanPacketDetect_80211ah(x, chanBW, varargin)

narginchk(2,4);
nargoutchk(0,2);

% Check if M is requested
if nargout==2
    M = [];
end

validateattributes(x, {'double'}, {'2d','finite'}, mfilename, 'signal input');
startOffset = [];

if isempty(x)
    return;
end

if nargin == 2
    threshold = 0.5; % Default value
    offset = 0;      % Default value
elseif nargin == 3
    threshold = 0.5;
    validateattributes(varargin{1}, {'double'}, {'integer','scalar','>=',0}, mfilename, 'offset');
    offset = varargin{1};
else
    validateattributes(varargin{1}, {'double'}, {'integer','scalar','>=',0}, mfilename, 'offset');
    validateattributes(varargin{2}, {'double'}, {'real','scalar','>',0,'<=',1}, mfilename, 'threshold');
    offset = varargin{1};
    threshold = varargin{2};
end

% Validate Offset value
coder.internal.errorIf(offset>size(x,1) - 1, 'wlan:shared:InvalidOffsetValue');

% Validate Channel Bandwidth
chanBW = wlan.internal.validateParam('S1GVHTCHANBW', chanBW, mfilename);

Td = 0.8e-5; % Time period of a short training symbol for 20MHz
                
switch chanBW
    case 'CBW4'
        symbolLength = Td/(1/4e6);
    case {'CBW8'}
        symbolLength = Td/(1/8e6);
    case {'CBW16'}
        symbolLength = Td/(1/16e6);
    otherwise % 'CBW5', 'CBW10', 'CBW20'
        symbolLength = Td/(1/2e6); % In samples
end

lenLSTF = symbolLength*10; % Length of 10 L-STF symbols 
lenHalfLSTF = lenLSTF/2;   % Length of 5 L-STF symbols
inpLength = (size(x,1) - offset); 

% Append zeros to make the input equal to multiple of L-STF/2
if inpLength<=lenHalfLSTF
    numPadSamples = lenLSTF - inpLength;
else
    numPadSamples = lenHalfLSTF*ceil(inpLength/lenHalfLSTF) - inpLength;
end
padSamples = zeros(numPadSamples, size(x,2));

% Process the input waveform in blocks of L-STF length. The processing
% blocks are offset by half the L-STF length.
numBlocks = (inpLength + numPadSamples)/lenHalfLSTF;

if nargout==2
% Define decision statistics vector
DS = coder.nullcopy(zeros(size(x,1) + length(padSamples) - offset -2*symbolLength + 1, 1));
    if numBlocks > 2
        for n=1:numBlocks-2
            % Update buffer
            buffer = x((n-1)*lenHalfLSTF + (1:lenLSTF) + offset, :);
            [startOffset, out] = correlateSamples(buffer, symbolLength, lenLSTF, threshold);

            DS((n-1)*lenHalfLSTF + 1:lenHalfLSTF*n, 1) = out(1:lenHalfLSTF);

            if ~(isempty(startOffset))
                % Packet detected
                startOffset = startOffset + (n-1)*lenHalfLSTF;
                DS((n-1)*lenHalfLSTF + (1:length(out)), 1) = out;
                % Resize decision statistics
                M = DS(1:(n-1)*lenHalfLSTF + length(out));
                return;
            end
        end
        % Process last block of data
        blkOffset = lenHalfLSTF*(numBlocks-2);
        buffer = [x(blkOffset + 1 + offset:end, :); padSamples];
        [startOffset, out] = correlateSamples(buffer, symbolLength, lenLSTF, threshold);
            if ~(isempty(startOffset))
                startOffset = startOffset + blkOffset; % Packet detected
            end
        DS(blkOffset + 1:end, 1) = out;
        M = DS(1:end-length(padSamples)); 
    else
        buffer = [x(offset + 1:end, :); padSamples];
        [startOffset, out] = correlateSamples(buffer, symbolLength, lenLSTF, threshold);
        M = out;
    end
else
    if numBlocks > 2
        for n=1:numBlocks-2
            buffer = x((n-1)*lenHalfLSTF + (1:lenLSTF) + offset, :); % Update buffer
            startOffset = correlateSamples(buffer, symbolLength, lenLSTF, threshold);

            if ~(isempty(startOffset))
                startOffset = startOffset + (n-1)*lenHalfLSTF; % Packet detected
                return;
            end
        end
    % Process last block of data
    blkOffset = lenHalfLSTF*(numBlocks-2);
    buffer = [x(blkOffset + 1 + offset:end, :); padSamples];
    startOffset = correlateSamples(buffer, symbolLength, lenLSTF, threshold);
        if ~(isempty(startOffset))
            startOffset = startOffset + blkOffset; % Packet detected
        end
    else
        buffer = [x(offset + 1:end, :); padSamples];
        startOffset = correlateSamples(buffer, symbolLength, lenLSTF, threshold); 
    end
end

end

function [packetStart,Mn] = correlateSamples(rxSig, symbolLength, lenLSTF, threshold)
%   Estimate the start offset of the preamble of the receive WLAN packet,
%   using auto-correlation method [1,2].

%   [1] OFDM Wireless LANs: A Theoretical and Practical Guide 1st Edition
%       by Juha Heiskala (Author),John Terry Ph.D. ISBN-13:978-0672321573
%   [2] OFDM Baseband Receiver Design for Wireless Communications by
%       Tzi-Dar Chiueh, Pei-Yun Tsai. ISBN: 978-0-470-82234-0

correlationLength = lenLSTF - (symbolLength*2);
pNoise = eps; % Adding noise to avoid the divide by zero
weights = ones(symbolLength, 1);
index = 1:correlationLength + 1;

packetStart = []; % Initialize output

% Shift data for correlation
rxDelayed = rxSig(symbolLength + 1:end , :); % Delayed samples
rx = rxSig(1:end-symbolLength, :);        % Actual samples

% Sum output on multiple receive antennas
C = sum(filter(weights, 1,(conj(rxDelayed).*rx)), 2);
CS = C(symbolLength:end)./symbolLength;

% Sum output on multiple receive antennas
P = sum(filter(weights, 1, (abs(rxDelayed).^2+abs(rx).^2)/2)./symbolLength, 2);

PS = P(symbolLength:end) + pNoise;

Mn = abs(CS).^2./PS.^2;
N = Mn > threshold;

if (sum(N) >= symbolLength*1.5)
    found = index(N);
    packetStart = found(1) - 1;
    % Check the relative distance between peaks relative to the first
    % peak. If this exceed three times the symbol length then the
    % packet is not detected.
    if sum((found(2:symbolLength) - found(1))>symbolLength*3)
        packetStart = [];
    end
end

end
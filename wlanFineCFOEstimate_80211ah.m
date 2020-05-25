function foffset = wlanFineCFOEstimate_80211ah(in, chanBW, varargin)
    
    % Validate number of arguments
    narginchk(2,3);

    % Validate the channel bandwidth
    chanBW = wlan.internal.validateParam('S1GVHTCHANBW',chanBW,mfilename);

    % Optional correlation offset
    if nargin>2
        corrOffset = varargin{1};
        validateattributes(corrOffset,{'numeric'},{'scalar','>=',0,'<=',1},mfilename,'CORROFFSET');
    else
        corrOffset = 0.75;
    end
    
    if any(strcmp(chanBW,{'CBW1', 'CBW2', 'CBW5','CBW10','CBW20'}))
    % Same FFT length for 1/2/5/10/20 MHz
        num20 = 1;
    else
        num20 = wlan.internal.cbwStr2Num(chanBW)/20;
    end
    
    FFTLen = 64*num20;
    Nltf = 160*num20;   % Number of samples in L-LTF
    fs = wlan.internal.cbwStr2Num(chanBW)*1e6;

    % Extract L-LTF or as many samples as we can
    lltf = in(1:min(Nltf,end),:);

    % Fine CFO estimate assuming one repetition per FFT period (2 OFDM symbols)
    M = FFTLen;             % Number of samples per repetition
    GI = FFTLen/2;          % Guard interval length
    S = M*2;                % Maximum useful part of L-LTF (2 OFDM symbols)
    N = size(lltf,1);       % Number of samples in the input

    % We need at most S samples
    offset = round(corrOffset*GI);
    use = lltf(offset+(1:min(S,N-offset)),:);
    foffset = wlan.internal.cfoEstimate(use,M).*fs/M;
end
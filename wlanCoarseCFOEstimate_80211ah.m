function foffset = wlanCoarseCFOEstimate_80211ah(in,chanBW,varargin)

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
        % Same FFT length for 5/10/20 MHz
        num20 = 1;
    else
        num20 = wlan.internal.cbwStr2Num(chanBW)/20;
    end
    FFTLen = 64*num20;
    Nstf = 160*num20; % Number of samples in L-STF
    fs = wlan.internal.cbwStr2Num(chanBW)*1e6;

    % Extract L-STF or as many samples as we can
    lstf = in(1:min(Nstf,end),:);

    % Coarse CFO estimate assuming 4 repetitions per FFT period
    M = FFTLen/4;           % Number of samples per repetition
    GI = FFTLen/4;          % Guard interval length
    S = M*9;                % Maximum useful part of L-STF
    N = size(lstf,1);       % Number of samples in the input

    % We need at most S samples
    offset = round(corrOffset*GI);
    use = lstf(offset+(1:min(S,N-offset)),:);
    foffset = wlan.internal.cfoEstimate(use,M).*fs/M;
end
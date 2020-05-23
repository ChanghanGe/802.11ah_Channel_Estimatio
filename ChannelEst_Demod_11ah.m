function [chanEst, rxPSDU] =  ChannelEst_Demod_11ah(rxSync, fieldInd, cfgS1G, ofdmInfo)
    rxLTF1 = rxSync(fieldInd.S1GLTF1(1):fieldInd.S1GLTF1(2),:);
    demodLTF1 = wlanS1GDemodulate(rxLTF1,'S1G-LTF1',cfgS1G);

    %Channel Estimation
    if cfgS1G.NumSpaceTimeStreams>1
        % Use S1G-LTF1 and S1G-LTF2N for channel estimation
        rxLTF2N = rxSync(fieldInd.S1GLTF2N(1):fieldInd.S1GLTF2N(2),:);
        demodLTF2N = wlanS1GDemodulate(rxLTF2N,'S1G-LTF2N',cfgS1G);
        chanEst = s1gLTFChannelEstimate([demodLTF1 demodLTF2N],cfgS1G);
    else
        % Use only S1G-LTF1 for channel estimation
        chanEst = s1gLTFChannelEstimate(demodLTF1,cfgS1G);
    end

    % Estimate Noise Level
    noiseVarEst = helperNoiseEstimate(demodLTF1);

    % Extract S1G-Data field
    rxData = rxSync(fieldInd.S1GData(1):fieldInd.S1GData(2),:);

    % OFDM demodulation
    demodSym = wlanS1GDemodulate(rxData,'S1G-Data',cfgS1G);

    % Extract data subcarriers from demodulated symbols and channel
    % estimate
    demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);
    chanEstData = chanEst(ofdmInfo.DataIndices,:,:);

    % MMSE frequency domain equalization
    [eqDataSym,csi] = helperSymbolEqualize(demodDataSym,chanEstData,noiseVarEst); 

    % Recover PSDU bits
    rxPSDU = s1gDataBitRecover(eqDataSym,noiseVarEst,csi,cfgS1G);
end 
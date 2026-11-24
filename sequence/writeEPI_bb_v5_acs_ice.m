function [seq, sys, res, gro, seq_acs_ice] = writeEPI_bb_v5_acs_ice(voxelSize, N, TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, type,...
    filename, use_gre_acs, gMax, sMax, dwell, save_path, rf_freqOffset, rfTB2use, rfDur2use, use_trig, delay_after_trig, trig, use_gz_excitation, T_gzReph, varargin)
% function [gro, adc] = writeSMSEPI(voxelSize, N, TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, varargin)
%
% Create .seq file for SMS-EPI or 3D EPI fMRI.
%
% The recommended use of this function to define a few frames, 
% then use the 'runs' (repeat) feature on the scanner to repeat this .seq file.
% 
% Uses skippedcaipi_sampling.py (R. Stirnberg) to define the ky-kz CAIPI sampling pattern.
%
% Inputs:
%   voxelSize    [1 3]    m
%   N            [1 3]    image matrix size 
%   TE           [1]      s
%   alpha        [1]      flip angle (degrees)
%   mb           [1]      multiband/SMS factor
%   pf_ky        [1]      partial Fourier factor
%   Ry           [1]      ky undersampling factor
%   Rz           [1]      kz undersampling factor
%   caipiShiftZ  [1]      caipi shift along z for every ky blip
%   nFrames      [1]      number of time frames (image volumes)
%   nDummyFrames [1]      number of frames w/o data acquisition to reach steady state
%   type         string   'SMS' or '3D'
%
% Keyword-argument input options:
%   gro          struct   readout gradient struct. Default: create a new one and return from this function
%   adc          struct   acquisition struct. Default: create a new one and return from this function
%
% Outputs:
%   gro 
%   adc

if nargin < 16
    gMax = 40;
end
if nargin < 17
    sMax = 110;
end
if nargin < 18
    dwell = 4e-6;
end
if nargin < 19
    save_path = pwd;
end
if nargin < 20
    rf_freqOffset = 0;
end
if nargin < 21
    rfTB2use = 30;
end
if nargin < 22
    rfDur2use = 8e-3;
end

[nx ny nz] = deal(N(1), N(2), N(3));


if strcmp(type, '3D') 
    mb = 1;
    np = nz / Rz;
else
    np = nz/mb;   % number of excitations/partitions (sets of SMS slices)
end
    

% parse input options
arg.gro = [];
arg.adc = [];
arg.fatsatOn = true;  
arg.simulateSliceProfile = true;    % simulate SMS profile and display
arg.gzPreOn = true;                 % prephase along kz by -mb/2*deltak
arg.plot = false;
arg.segmentRingdownTime = 120e-6;   % segment ringdown time for Pulseq on GE
arg.doRefScan = true;               % do EPI ghost reference scan in frame 1
arg.caipiPythonPath = '~/github/HarmonizedMRI/3DEPI/caipi/';

arg = toppe.utils.vararg_pair(arg, varargin);

warning('OFF', 'mr:restoreShape');

% Define experimental parameters
sys = mr.opts('maxGrad', gMax, 'gradUnit','mT/m', ...
              'maxSlew', sMax, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 20e-6, ...
              'adcRasterTime', 2e-6, ...
              'gradRasterTime', 10e-6, ...
              'blockDurationRaster', 10e-6, ...
              'B0', 3.0);

% reduced slew for spoilers
sys2 = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 60, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', sys.rfDeadTime, ...
              'rfRingdownTime', sys.rfRingdownTime, ...
              'adcDeadTime', sys.adcDeadTime, ...
              'gradRasterTime', sys.gradRasterTime, ...
              'blockDurationRaster', sys.blockDurationRaster, ...
              'B0', sys.B0);

fov = voxelSize .* [nx ny nz];       % FOV (m)

if strcmp(type, 'SMS')
    slThick = fov(3)/nz;
else
    % 3d encoding: 85% of FOV is targeted, 15% probably oversampled
    % slThick = 0.85*fov(3);

    slThick = 0.95 * fov(3);
end

% this needs to change to allow shorter echo spacings
if nargin < 18
    dwell = 4e-6;                   % ADC sample time (s). For GE, must be multiple of 2us. For Siemens, must be multiple of 0.1 us and at least 1 us
end

etl = 2*ceil(pf_ky*ny/Ry/2);    % echo train length. even
% BB: i think this may need to be floor instead of ceil, since it gives an error
% when computing: kyInds = double(indices((end-etl+1):end, 2));
% make sure that matrix size is divisible by (Ry*2) to avoid this

nCyclesSpoil = 2;               % number of spoiler cycles, along x and z
rfSpoilingInc = 117;            % RF spoiling increment (degrees)

if strcmp(type, 'SMS')
    rfTB  = 6;          % RF pulse time-bandwidth product
    rfDur = 8e-3;       % RF pulse duration (s)
else
    % rfTB  = 8;          % RF pulse time-bandwidth product
    % rfDur = 2e-3;       % RF pulse duration (s)

    rfTB = rfTB2use;
    rfDur = rfDur2use;
end


% Create fat sat pulse 
fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3.5;     % time-bandwidth product
fatsat.dur     = 8.0;     % pulse duration (ms)

% RF waveform in Gauss
wav = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, toppe.systemspecs(), ...
    'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'ftype', 'min', ...
    'writeModFile', false);

% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
rfp = rf2pulseq(wav, 4e-6, sys.rfRasterTime);

% Create pulseq object
% Try to account for the fact that makeArbitraryRf scales the pulse as follows:
% signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi);
flip = fatsat.flip/180*pi;

rfsat = mr.makeArbitraryRf(rfp, ...
    flip*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'system', sys);
rfsat.signal = rfsat.signal/max(abs(rfsat.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
rfsat.freqOffset = -440;  % Hz


% Create excitation pulse
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', 0.7*sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.25);

sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)

if strcmp(type, 'SMS')
    ftype = 'ls';
else
    ftype = 'min';
end

fig_num=201;

% Jon's version of slab select's gzRF already includes rephaser
% [rf, gzRF, freq] = getsmspulse_v0(alpha, slThick, rfTB, rfDur, ...
%     mb, sliceSep, sysGE, sys, fig_num, ...
%     'doSim', arg.simulateSliceProfile, ...      % Plot simulated SMS slice profile
%     'type', 'st', ...                           % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
%     'ftype', ftype);                            % filter design. 'ls' = least squares

% freq = freq/sliceSep*voxelSize(3); % frequency offset for z shift of voxelSize(3)

% to match 3d-gre RF design, use sinc pulse:
% Create alpha-degree slice selection pulse and gradient
[rf, gzRF] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDur,...
    'SliceThickness',slThick,'apodization',0.42,'timeBwProduct',rfTB,'system',sys);

% refocusing Gz gradient -> this is required for the sinc pulse version
gzReph = mr.makeTrapezoid('z','Area',-gzRF.area/2,'Duration',T_gzReph,'system',sys);



% Get CAIPI sampling pattern (for one shot/echo train)

% create caipi.mat, and load it
pyFile = [arg.caipiPythonPath 'skippedcaipi_sampling.py'];


% BB Hack: use the created python environment from vs code
% pyCmd = sprintf('/cluster/berkin/berkin/miniconda/install/envs/skip_caipi_env/bin/python %s %d %d %d %d %d %d', pyFile, ny, nz, Ry, Rz, caipiShiftZ, 1);

pyCmd = sprintf('python3 %s %d %d %d %d %d %d', pyFile, ny, nz, Ry, Rz, caipiShiftZ, 1);


if system(pyCmd) ~= 0
    fprintf('Open a terminal and run the following python command:\n\t%s\n', pyCmd);
    input('\tWhen done, press Enter to continue');
end

load caipi

% kz and ky indeces (multiples of deltak)
kyInds = double(indices((end-etl+1):end, 2));
kzInds = double(indices((end-etl+1):end, 1));

% ky/kz encoding blip amplitude along echo train (multiples of deltak)
kyStep = diff(kyInds);
kzStep = diff(kzInds);

kyStepMax = max(abs(kyStep));
kzStepMax = max(abs(kzStep));


% Define readout gradients and ADC event

deltak = 1./fov;


% Start with the blips
% Wrap in trap4ge.m so all times are on 20us boundary for accurate interpolation to 4us for GE
commonRasterTime = 20e-6;   
gyBlip = trap4ge(mr.makeTrapezoid('y', sys, 'Area', kyStepMax*deltak(2)), commonRasterTime, sys);
gzBlip = trap4ge(mr.makeTrapezoid('z', sys, 'Area', kzStepMax*deltak(3)), commonRasterTime, sys);

blipDuration = max(mr.calcDuration(gyBlip), mr.calcDuration(gzBlip));
maxBlipArea = max(gyBlip.area, gzBlip.area);

% Readout trapezoid
if isempty(arg.gro) 
    systmp = sys;

    % this seems to be limiting the echo spacing
    systmp.maxGrad = deltak(1)/dwell;  % to ensure >= Nyquist sampling
    
    gro = trap4ge(mr.makeTrapezoid('x', systmp, 'Area', nx*deltak(1) + maxBlipArea), commonRasterTime, systmp);
else
    gro = arg.gro;
end

% ADC event
% Number of readout samples must be multiple of 4 (TODO: check if this is actually needed)
% because of manually input dwell time and divisibility considerations,
% read FOV may end up bigger than speficied, which should be fine 
 
if isempty(arg.adc) 
    TreadTmp = mr.calcDuration(gro) - blipDuration;
    numSamplesTmp = ceil(TreadTmp/dwell);
    numSamples = ceil(numSamplesTmp/4)*4;
    Tread = numSamples*dwell;
    
    adc = mr.makeAdc(numSamples, sys, ...
        'Duration', Tread, ...
        'Delay', blipDuration/2);  % + dwell/2? TODO
else
    adc = arg.adc;
end

% Split blips at block boundary
[gyBlipUp, gyBlipDown] = mr.splitGradientAt(gyBlip, mr.calcDuration(gyBlip)/2);
gyBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gyBlipUp);
gyBlipDown.delay = 0;

[gzBlipUp, gzBlipDown] = mr.splitGradientAt(gzBlip, mr.calcDuration(gzBlip)/2);
gzBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gzBlipUp);
gzBlipDown.delay = 0;

% prephasers and spoilers
gxPre = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2), ...
    commonRasterTime, sys);
Tpre = mr.calcDuration(gxPre) + 4*commonRasterTime;
gyPre = trap4ge(mr.makeTrapezoid('y', sys, ...
    'Area', (kyInds(1)-ny/2)*deltak(2), ... 
    'Duration', Tpre-4*commonRasterTime), ...  % make a bit shorter than Tpre to ensure duration doesn't exceed Tpre after trap4ge
    commonRasterTime, sys);

if strcmp(type, 'SMS')
    area = -floor(mb/2)*deltak(3);
else
    area = nz/2*deltak(3);
end

gzPre = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', area, ...
    'Duration', Tpre-commonRasterTime), ...   % make < Tpre to ensure duration doesn't exceed Tpre
    commonRasterTime, sys);

gxSpoil = mr.makeTrapezoid('x', sys2, ...
    'Area', -nx*deltak(1)*nCyclesSpoil);

gzSpoil = mr.makeTrapezoid('z', sys2, ...
    'Area', nx*deltak(1)*nCyclesSpoil);


% Calculate delays to achieve desired TE and TR
kyIndAtTE = find(kyInds-ny/2 == min(abs(kyInds-ny/2)));

% this version assumes gzRF comes with rephasing lobe:
% minTE = mr.calcDuration(gzRF) - mr.calcDuration(rf)/2 - rf.delay + mr.calcDuration(gxPre) + ...
        % (kyIndAtTE-0.5) * mr.calcDuration(gro);

% this version now includes Gz rephaser's timing
minTE = mr.calcDuration(gzRF) + mr.calcDuration(gzReph) - mr.calcDuration(rf)/2 - rf.delay + mr.calcDuration(gxPre) + ...
        (kyIndAtTE-0.5) * mr.calcDuration(gro);


if TE < minTE
    error(sprintf('Requested TE < minimum TE (%d)', minTE));
end

TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;

% minTR = arg.fatsatOn*(mr.calcDuration(rfsat) + mr.calcDuration(gxSpoil)) + ...
%     mr.calcDuration(gzRF) + TEdelay + mr.calcDuration(gxPre) + ...
%     etl*mr.calcDuration(gro) + mr.calcDuration(gxSpoil);

% same as above, include gz rephaser in the timing
minTR = arg.fatsatOn*(mr.calcDuration(rfsat) + mr.calcDuration(gxSpoil)) + ...
    mr.calcDuration(gzRF) + mr.calcDuration(gzReph) + TEdelay + mr.calcDuration(gxPre) + ...
    etl*mr.calcDuration(gro) + mr.calcDuration(gxSpoil);

% new addition in v3:
if strcmp(type, '3D')
    minTR = minTR + mr.calcDuration(gzPre);
end


if TR < np*minTR
    error(sprintf('Requested TR < minimum TR (%d),  slice/partition min TR = %d', np*minTR, minTR));
end

TRdelay = round((TR/np-minTR-arg.segmentRingdownTime)/sys.blockDurationRaster) * sys.blockDurationRaster;

if TRdelay < 0
    assert('error: TR delay is negative')    
else
    disp(['TR delay: ', num2str(1e3*TRdelay), ' ms'])
end

% Interleaved slice ordering for SMS
if strcmp(type, 'SMS')
    % if even number of shots, swap last two to reduce slice cross-talk
    IP = [1:2:np 2:2:np];
    if ~mod(np,2)
        % for np = even, change order of last two partitions/shots
        l = length(IP);
        IP = IP([1:(l-2) l l-1]);
    end
else
    IP = -nz/2:Rz:nz/2-1;
end


% Assemble sequence
seq = mr.Sequence(sys);           

% separate seq file with flipped phase encoding table for ice recon for acs
seq_acs_ice = mr.Sequence(sys);           


%--------------------------------------------------------------------------
% GRE acs parameters
%--------------------------------------------------------------------------

if use_gre_acs
    bandwidth = 250;        % Hz/pixel

    % when Tread is defined as 1/bandwidth, there is "an assertion failed"
    % message when writing the sequence at the end (for bandwidth=300)
    % maybe Tread becomes 3.33... ms, which causes some problem?
    Tread_acs = 1 / bandwidth;  % readout duration in sec

    Tpre = 3e-3;            % duration of Gz/Gy gradient blips
    Ndummy = 50;            % dummy shots to enter steady state


    %--------------------------------------------------------------------------
    % calculate timing: ACS
    %--------------------------------------------------------------------------

    TE_acs = 15e-3;
    TR_acs = 25e-3;
%     TE_acs = 7e-3;
%     TR_acs = 21e-3;

    % Ny_acs = 32;
    % Nz_acs = 32;
    % increase resolution for online recon
    Ny_acs = 48;
    Nz_acs = 48;


    rf_acs = rf;

    % include frequency offset in the acs
    rf_acs.freqOffset = rf_freqOffset;  % Hz


    gzRF_acs = gzRF;

    % need to include Gz rephaser
    gzReph_acs = gzReph;


    % Define other gradients and ADC events
    deltak = 1 ./ fov;

    gx_acs = mr.makeTrapezoid('x', sys, 'FlatArea', nx*deltak(1), 'FlatTime', Tread_acs);    % readout gradient

    gxPre_acs = mr.makeTrapezoid('x', sys, 'Area', -gx_acs.area/2, 'Duration', Tpre);        % Gx prewinder

    gxSpoil_acs = mr.makeTrapezoid('x', sys, 'Area', gx_acs.area, 'Duration', Tpre);         % Gx spoiler
    % gxSpoil_acs = mr.makeTrapezoid('x', sys, 'Area', 2*gx_acs.area, 'Duration', Tpre);         % Gx spoiler -> try doubling the spoiler area


    % ACS adc object: match the number of samples in the epi readout
    % adc = mr.makeAdc(numSamples, sys, 'Duration', Tread, 'Delay', blipDuration/2);  % + dwell/2? TODO
    adc_acs = mr.makeAdc(numSamples, 'Duration', gx_acs.flatTime, 'Delay', gx_acs.riseTime);            % adc event


    areaY = ((0:ny-1)-ny/2)*deltak(2);
    areaZ = ((0:nz-1)-nz/2)*deltak(3);


    % update timing with Gz rephaser
    % delayTE_acs = ceil( (TE_acs - mr.calcDuration(rf) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre_acs)  ...
        % - mr.calcDuration(gx_acs)/2)/seq.gradRasterTime)*seq.gradRasterTime;
    delayTE_acs = ceil( (TE_acs - mr.calcDuration(rf) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre_acs) - mr.calcDuration(gzReph_acs) ...
        - mr.calcDuration(gx_acs)/2)/seq.gradRasterTime)*seq.gradRasterTime;

    if delayTE_acs < 0
        disp(['error: acs delay TE is negative, ', num2str(delayTE_acs*1e3), ' ms'])    
        assert(delayTE_acs>0)    
    else
        disp(['acs delay TE: ', num2str(1e3*delayTE_acs), ' ms'])
    end

    % update here as well
    % delayTR_acs = ceil((TR_acs - mr.calcDuration(rf) - mr.calcDuration(gxPre_acs) ...
        % - mr.calcDuration(gx_acs) - mr.calcDuration(gxSpoil_acs) - delayTE_acs)/seq.gradRasterTime)*seq.gradRasterTime;
    delayTR_acs = ceil((TR_acs - mr.calcDuration(rf) - mr.calcDuration(gxPre_acs) - mr.calcDuration(gzReph_acs) ...
        - mr.calcDuration(gx_acs) - mr.calcDuration(gxSpoil_acs) - delayTE_acs)/seq.gradRasterTime)*seq.gradRasterTime;

    if delayTR_acs < 0
        disp(['error: acs delay TR is negative, ', num2str(delayTR_acs*1e3), ' ms'])    
        assert(delayTR_acs>0)  
    else
        disp(['acs delay TR: ', num2str(delayTR_acs*1e3), ' ms'])
    end


    %--------------------------------------------------------------------------
    % ACS data: dummies
    %--------------------------------------------------------------------------

    gyPre_acs = mr.makeTrapezoid('y', 'Area', areaY(floor(ny/2)), 'Duration', Tpre);
    gyReph_acs = mr.makeTrapezoid('y', 'Area', -areaY(floor(ny/2)), 'Duration', Tpre);

    rf_phase = 0;
    rf_inc = 0;

    % Drive magnetization to steady state
    for iY = 1:Ndummy
        if use_trig
            seq.addBlock(trig);
            seq.addBlock(mr.makeDelay(delay_after_trig));
        end

        % excitation pulse
        rf_acs.phaseOffset = rf_phase/180*pi;
        
        if use_gz_excitation
            seq.addBlock(rf_acs,gzRF_acs);
            seq.addBlock(gzReph_acs);                  
        else
            seq.addBlock(rf_acs);
            % add a delay (of 5 ms) for Jason's ac/dc
            seq.addBlock(mr.makeDelay(T_gzReph));
        end

        rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);

        % Gradients    
        seq.addBlock(gxPre_acs, gyPre_acs);                  % add Gx pre-winder, go to desired ky
        seq.addBlock(mr.makeDelay(delayTE_acs));    % add delay needed before the start of readout

        seq.addBlock(gx_acs);                           % add readout Gx

        seq.addBlock(gyReph_acs, gxSpoil_acs);               % add Gx spoiler, and go back to DC in ky
        seq.addBlock(mr.makeDelay(delayTR_acs))     % add delay to the end of TR
    end


    temp = 1:ny;
    iY_acs_indices = temp(1+end/2-Ny_acs/2:end/2+Ny_acs/2);
   
    temp = 1:nz;
    iZ_acs_indices = temp(1+end/2-Nz_acs/2:end/2+Nz_acs/2);
    
    % Make trapezoids for inner loop to save computation
    for iY = iY_acs_indices
        gyPre_acs(iY) = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);
        gyReph_acs(iY) = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre);
    end

    %--------------------------------------------------------------------------
    % ACS data: Loop over phase encodes and define sequence blocks
    %--------------------------------------------------------------------------

    rf_phase = 0;
    rf_inc = 0;

    for iZ = iZ_acs_indices
        % Gz blips to go desired kz, and to come back to DC in kz
        gzPre_acs = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre);  

        % change the name of this negative blip to avoid overlap with Gz
        % rephaser after the RF
        gzRephBlip_acs = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',Tpre);

        disp(['kZ: ', num2str(iZ), ' kY: ', num2str(iY_acs_indices(1)), ' to ', num2str(iY_acs_indices(end))])

        for iY = iY_acs_indices
            if use_trig
                seq.addBlock(trig);
                seq.addBlock(mr.makeDelay(delay_after_trig));
            end

            % excitation pulse
            rf_acs.phaseOffset = rf_phase/180*pi;
            adc_acs.phaseOffset = rf_phase/180*pi;

            if use_gz_excitation
                seq.addBlock(rf_acs,gzRF_acs);
    
                seq.addBlock(gzReph_acs);                
            else
                seq.addBlock(rf_acs);
    
                % add a delay (of 5 ms) for Jason's ac/dc
                seq.addBlock(mr.makeDelay(T_gzReph));
            end

            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);

            % Encoding
            seq.addBlock(gxPre_acs, gyPre_acs(iY), gzPre_acs);        % Gz, Gy blips, Gx pre-winder
            seq.addBlock(mr.makeDelay(delayTE_acs));    % delay until readout

            seq.addBlock(gx_acs, adc_acs);                       % readout

            seq.addBlock(gyReph_acs(iY), gzRephBlip_acs, gxSpoil_acs);% -Gz, -Gy blips, Gx spoiler

            seq.addBlock(mr.makeDelay(delayTR_acs));
        end
    end

    disp(['ACS scan time: ', num2str((Ndummy + Ny_acs * Nz_acs) * TR_acs), ' sec'])



    %--------------------------------------------------------------------------
    % dummies for separate gre acs scan with ice recon
    %--------------------------------------------------------------------------

    gyPre_acs = mr.makeTrapezoid('y', 'Area', areaY(floor(ny/2)), 'Duration', Tpre);
    gyReph_acs = mr.makeTrapezoid('y', 'Area', -areaY(floor(ny/2)), 'Duration', Tpre);

    rf_phase = 0;
    rf_inc = 0;

    for iY = 1:Ndummy
        if use_trig
            seq_acs_ice.addBlock(trig);
            seq_acs_ice.addBlock(mr.makeDelay(delay_after_trig));
        end

        % excitation pulse
        rf_acs.phaseOffset = rf_phase/180*pi;
        
        if use_gz_excitation
            seq_acs_ice.addBlock(rf_acs,gzRF_acs);
            seq_acs_ice.addBlock(gzReph_acs);                  
        else
            seq_acs_ice.addBlock(rf_acs);
            % add a delay (of 5 ms) for Jason's ac/dc
            seq_acs_ice.addBlock(mr.makeDelay(T_gzReph));
        end

        rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);

        % Gradients   
        seq_acs_ice.addBlock(gxPre_acs, gyPre_acs);                  % add Gx pre-winder, go to desired ky
        seq_acs_ice.addBlock(mr.makeDelay(delayTE_acs));    % add delay needed before the start of readout

        seq_acs_ice.addBlock(gx_acs);                           % add readout Gx

        seq_acs_ice.addBlock(gyReph_acs, gxSpoil_acs);               % add Gx spoiler, and go back to DC in ky
        seq_acs_ice.addBlock(mr.makeDelay(delayTR_acs))     % add delay to the end of TR
    end


    % this flipping was necessary to get the dicoms to work, otherwise they
    % were flipped:
    iY_acs_indices_ice = flipdim(iY_acs_indices,2);
    iZ_acs_indices_ice = flipdim(iZ_acs_indices,2);
    % this requires twix data to be flipped back so that it is consistent
    % with the EPI imaging data: img_acs_unflipped = circshift(flip(flip(img_acs,2),3),[0,1,1]);

    seq_acs_ice.addBlock(mr.makeLabel('SET','REV', 1)); % left-right swap fix (needed for 1.4.0)

    % for 3d imaging, use PAR
    seq_acs_ice.addBlock(mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','PAR', 0)); 

    for iY = iY_acs_indices_ice
        gyPre_acs(iY) = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);
        gyReph_acs(iY) = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre);
    end

    rf_phase = 0;
    rf_inc = 0;

    for iZ = iZ_acs_indices_ice
        % Gz blips to go desired kz, and to come back to DC in kz
        gzPre_acs = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre);  

        % change the name of this negative blip to avoid overlap with Gz
        % rephaser after the RF
        gzRephBlip_acs = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',Tpre);

        disp(['kZ: ', num2str(iZ), ' kY: ', num2str(iY_acs_indices_ice(1)), ' to ', num2str(iY_acs_indices_ice(end))])

        for iY = iY_acs_indices_ice
            if use_trig
                seq_acs_ice.addBlock(trig);
                seq_acs_ice.addBlock(mr.makeDelay(delay_after_trig));
            end

            % excitation pulse
            rf_acs.phaseOffset = rf_phase/180*pi;
            adc_acs.phaseOffset = rf_phase/180*pi;

            if use_gz_excitation
                seq_acs_ice.addBlock(rf_acs,gzRF_acs);
    
                seq_acs_ice.addBlock(gzReph_acs);                
            else
                seq_acs_ice.addBlock(rf_acs);
    
                % add a delay (of 5 ms) for Jason's ac/dc
                seq_acs_ice.addBlock(mr.makeDelay(T_gzReph));
            end

            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);

            % Encoding
            seq_acs_ice.addBlock(gxPre_acs, gyPre_acs(iY), gzPre_acs);        % Gz, Gy blips, Gx pre-winder
            seq_acs_ice.addBlock(mr.makeDelay(delayTE_acs));    % delay until readout

            seq_acs_ice.addBlock(gx_acs, adc_acs);                       % readout

            seq_acs_ice.addBlock(gyReph_acs(iY), gzRephBlip_acs, gxSpoil_acs);% -Gz, -Gy blips, Gx spoiler

            spoilBlockContents = {mr.makeDelay(delayTR_acs)}; 

            if iY ~= iY_acs_indices_ice(end)  
                spoilBlockContents = [spoilBlockContents {mr.makeLabel('INC','LIN', 1)}];
            else
                spoilBlockContents = [spoilBlockContents {mr.makeLabel('SET','LIN', 0), mr.makeLabel('INC','PAR', 1)}];
            end
    
            seq_acs_ice.addBlock(spoilBlockContents{:});
        end
    end

end



%--------------------------------------------------------------------------
% EPI imaging data
%--------------------------------------------------------------------------

% temporal frame loop
rf_phase = 0;
rf_inc = 0;
msg = [];

for ifr = (1-nDummyFrames):nFrames

    for ii = 1:length(msg)
        fprintf('\b');
    end
    msg = sprintf('Frame %d of %d     ', ifr, nFrames);
    fprintf(msg);

    % dummy shots before turning on ADC, to reach steady state 
    isDummyShot = ifr < 1;

    % First frame is EPI calibration/reference scan (blips off)
    isRefShot = ifr == 1 & arg.doRefScan;

    % Segment ID (see 'Pulseq on GE' manual)
    segmentID = 3;

    if isDummyShot
        segmentID = 1;
    end

    if isRefShot
        segmentID = 2;
    end

    blipsOn = ~isDummyShot & ~isRefShot;

    % slice offset (partition/SMS group) loop
    for p = IP
        if strcmp(type, 'SMS')
            rf.freqOffset = round((p-np/2-1)*freq);  % frequency offset (Hz) for slice shift
        else
            % rf.freqOffset = 0;
            rf.freqOffset = rf_freqOffset;  % Hz
        end

        % Must label the first block in segment with segment ID. See Pulseq on GE manual.
        if arg.fatsatOn
            % fat sat
            rfsat.phaseOffset = rf_phase/180*pi;

            % seq.addBlock(rfsat, mr.makeLabel('SET', 'LIN', segmentID));
            seq.addBlock(rfsat);

            seq.addBlock(gxSpoil, gzSpoil);
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);

            if use_trig
                seq.addBlock(trig);
                seq.addBlock(mr.makeDelay(delay_after_trig));
            end

            % excitation pulse
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;

            if use_gz_excitation
                seq.addBlock(rf,gzRF);

                % Gz rephaser which was missing 
                seq.addBlock(gzReph);                  
            else
                seq.addBlock(rf);

                % add a delay (of 5 ms) for Jason's ac/dc
                seq.addBlock(mr.makeDelay(T_gzReph));
            end

            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
        else
            % don't do fat sat
            if use_trig
                seq.addBlock(trig);
                seq.addBlock(mr.makeDelay(delay_after_trig));
            end

            % excitation pulse
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;

            if use_gz_excitation    
                % seq.addBlock(rf, gzRF, mr.makeLabel('SET', 'LIN', segmentID));
                seq.addBlock(rf, gzRF);

                % Gz rephaser which was missing 
                seq.addBlock(gzReph);                  
            else
                % seq.addBlock(rf, mr.makeLabel('SET', 'LIN', segmentID));
                seq.addBlock(rf);

                % add a delay (of 5 ms) for Jason's ac/dc
                seq.addBlock(mr.makeDelay(T_gzReph));
            end

            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
        end

        % TE delay
        seq.addBlock(mr.makeDelay(TEdelay));

        % Readout
        %seq.addBlock(gxPre, gyPre, mr.scaleGrad(gzPre, 1 - 2/nz*(shot-1)));
        if strcmp(type, 'SMS')
            if caipiShiftZ > 0 & arg.gzPreOn
                seq.addBlock(gxPre, mr.scaleGrad(gyPre, blipsOn), mr.scaleGrad(gzPre, blipsOn));
            else
                seq.addBlock(gxPre, mr.scaleGrad(gyPre, blipsOn));
            end
        else
            seq.addBlock(gxPre, mr.scaleGrad(gyPre, blipsOn), mr.scaleGrad(gzPre, p/(nz/2)*blipsOn));
        end

        if isDummyShot
            seq.addBlock(gro, ...
                         mr.scaleGrad(gyBlipUp, blipsOn*kyStep(1)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipUp, blipsOn*kzStep(1)/max(kzStepMax,1)));
        else
            seq.addBlock(gro, adc, ...
                         mr.scaleGrad(gyBlipUp, blipsOn*kyStep(1)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipUp, blipsOn*kzStep(1)/max(kzStepMax,1)));
        end

        for ie = 2:(etl-1)
            % Example: 'gybdu' = Gy blip down up
            gybd = mr.scaleGrad(gyBlipDown, blipsOn*kyStep(ie-1)/max(kyStepMax,1));
            gybu = mr.scaleGrad(gyBlipUp,   blipsOn*kyStep(ie)/max(kyStepMax,1));
            gybdu = mr.addGradients({gybd, gybu}, sys);
            gzbd = mr.scaleGrad(gzBlipDown, blipsOn*kzStep(ie-1)/max(kzStepMax,1));
            gzbu = mr.scaleGrad(gzBlipUp,   blipsOn*kzStep(ie)/max(kzStepMax,1));
            gzbdu = mr.addGradients({gzbd, gzbu}, sys);

            if isDummyShot
                seq.addBlock(mr.scaleGrad(gro, (-1)^(ie-1)), gybdu, gzbdu);
            else
                seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)), gybdu, gzbdu);
            end
        end

        if isDummyShot
            seq.addBlock(mr.scaleGrad(gro, (-1)^(ie)), ...
                         mr.scaleGrad(gyBlipDown, blipsOn*kyStep(ie)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipDown, blipsOn*kzStep(ie)/max(kzStepMax,1)));
        else
            seq.addBlock(adc, ...
                         mr.scaleGrad(gro, (-1)^(ie)), ...
                         mr.scaleGrad(gyBlipDown, blipsOn*kyStep(ie)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipDown, blipsOn*kzStep(ie)/max(kzStepMax,1)));
        end

        % new addition in v3:
        % gz rephaser for 3D
        if strcmp(type, '3D')
            % seq.addBlock(mr.scaleGrad(gzPre, -p/(nz/2)*zBlipsOn*arg.gzPreOn));
            seq.addBlock(mr.scaleGrad(gzPre, -p/(nz/2)*blipsOn*arg.gzPreOn));
        end

        % spoil. Disabling makes the k-space plot neater
        seq.addBlock(gxSpoil, gzSpoil);

        % TR delay
        seq.addBlock(mr.makeDelay(TRdelay));

    end
end
fprintf('\n');


% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end


% Write .seq file
if strcmp(type, 'SMS')
    fnStem = 'sms-epi';
else
    fnStem = '3d-epi';
end

seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'epi');

if nargin < 19
    if use_gre_acs
        ofn = [fnStem, filename, '_use_gre_acs.seq'];
    else
        ofn = [fnStem, filename, '.seq'];
    end
else
    if use_gre_acs
        ofn = [save_path, fnStem, filename, '_use_gre_acs.seq'];
    else
        ofn = [save_path, fnStem, filename, '.seq'];
    end
end

seq.write(ofn);       % Write to pulseq file

% this seems to be necessary to get the ice to work:
seq_acs_ice.setDefinition('FOV', fov);
seq_acs_ice.setDefinition('Name', 'gre_acs_ice');

seq_acs_ice.write([save_path, fnStem, filename, '_only_gre_acs_ice.seq']);

res = [];
res.nz = nz;
res.mb = mb;
res.np = np;

return

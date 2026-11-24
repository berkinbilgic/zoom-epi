%--------------------------------------------------------------------------
%% MATLAB toolboxes
%--------------------------------------------------------------------------

clear
restoredefaultpath

addpath /cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ARES_zoom_fmri/pulseq-1.5.0/matlab/
addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ARES_zoom_fmri

addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/toppe-main
addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/SMS-EPI-develop/sequence/Pulseq/
addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/PulCeq-main/matlab/


caipiPythonPath = '/autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/3DEPI-main/caipi/';


%--------------------------------------------------------------------------
% create 3D-EPI sequence files
%--------------------------------------------------------------------------

% acquisition parameters
voxelSize = [1, 1, 1] * 1.e-3;       % m
% FOV = [220, 220, 130] * 1e-3;       % m
FOV = [220, 66, 64] * 1e-3;       % m

mtx_size = round(FOV ./ voxelSize);
mtx_size = mtx_size + rem(mtx_size,2);

encoding = '3D';

Ry = 1;
% Ry = 2;
% Ry = 3;

pf_ky = 1;
% pf_ky = 7/8;                % partial fourier amount

Rz = 1;
% Rz = 2;
% Rz = 3;
mb = [];


% skyra
gMax = 28;
sMax = 120;


% caipiShiftZ = 2 * (Rz>1);
caipiShiftZ = 0 * (Rz>1);

if Ry == 1
    if pf_ky == 1 
        TE = 47e-3;                         % sec [R=1]
        TR = 7040e-3 / Rz;                     % volume TR (sec) -- not the "inner" TR in 3D encoding [R=1]
    else
        TE = 92e-3;                         % sec [R=1]
        TR = 29000e-3 / Rz;                     % volume TR (sec) -- not the "inner" TR in 3D encoding [R=1]
    end
else 
    if Ry == 2
        if pf_ky == 1 
            TE = 71e-3;                         % sec [R=2] w/o pf
            TR = 19500e-3 / Rz;                     % volume TR (sec) -- not the "inner" TR in 3D encoding [R=2]
        else
            TE = 48e-3;                         % sec [R=2] with pf=7/8
            TR = 15500e-3 / Rz;                     % volume TR (sec) -- not the "inner" TR in 3D encoding [R=2]
        end
    else
        if Ry == 3
            if pf_ky == 1 
                TE = 23e-3;                         % sec [R=3] w/o pf
                TR = 14500e-3 / Rz;                     % volume TR (sec) -- not the "inner" TR in 3D encoding [R=3]
            else
                TE = 33e-3;                         % sec [R=3] with pf=7/8
                TR = 11500e-3 / Rz;                     % volume TR (sec) -- not the "inner" TR in 3D encoding [R=3]
            end
        end
    end
end


% make sure mtx_size(2) is divisible by (Ry*2)
mtx_size_update = ceil(mtx_size(2) / Ry / 2) * Ry * 2;
mtx_size(2) = mtx_size_update;

% update FOV accordingly
FOV = voxelSize .* mtx_size;

disp(['FOV: ',  num2str(FOV), '  mtx_size: ', num2str(mtx_size)])


T1 = 1.25;                     % assume tissue T1

TR_partition = TR / (mtx_size(3) / Rz);
alpha = acosd(exp(-TR_partition/T1));

disp(['TR partition: ', num2str(TR_partition*1e3), ' ms  FA: ', num2str(alpha), '  caipi shift: ', num2str(caipiShiftZ)])


% Design sub-sequence containing 4 frames, to ensure that RF phase cycling (spoiling) 
% is not interrupted across runs (for rf_inc = 117).

% use the following to get all calibration data:
% nDummyFrames = 1*1;       % to get into steady state
% nFrames = 1*2;            % first frame is used for ghost correction if 'doRefScan' is true
% 
% doRefscan = true;         % set to false to remove ghost correction TR 
% use_gre_acs = 0*1;        % set to 0 to disable gre acs refscan

% use the following to get only the EPI fMRI data:
nDummyFrames = 0*1;       % to get into steady state
nFrames = 1*1;            % first frame is used for ghost correction if 'doRefScan' is true

doRefscan = false;        % set to false to remove ghost correction TR 
use_gre_acs = 0;          % set to 0 to disable gre acs refscan



% directory to save seq file to:
seq_path = ['/autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ARES_zoom_fmri/', datestr(datetime('today')), '/' ];

if ~isfolder(seq_path)
    mkdir(seq_path)
end

% filename to use
filename = ['_RyxRz_', num2str(Ry), 'x', num2str(Rz), 'del', num2str(caipiShiftZ), '_pf', num2str(pf_ky)];


%--------------------------------------------------------------------------
% calculate dwell time needed -> Siemens needs it to be multiples of 0.1 usec
%--------------------------------------------------------------------------

gamma_over_2pi = 42.577; % MHz/T
deltak(1) = 1 / FOV(1);  % m-1
dwell = deltak(1) / gMax / gamma_over_2pi * 1e-3;   % sec

disp(['required dwell time to use Gmax during readout: ', num2str(1e6 * dwell), ' usec'])

dwell2use = max(round(10 * dwell * 1e6) / 10 / 1e6, 1e-6);

disp(['used dwell time to use Gmax during readout: ', num2str(1e6 * dwell2use), ' usec'])


%--------------------------------------------------------------------------
% generate both acs and imaging scan, write gre acs with ice to separate seq file
%--------------------------------------------------------------------------

% excitation to match 3d-gre
rf_dur = 8e-3;                      % duration of RF
rf_tbp = 30;                        % time bandwidth product for RF
rf_freqOffset = 1*3e3;              % Hz


rfTB = rf_tbp;
rfDur = rf_dur;

use_trig = 1;                      % set to 1 to add trigger prior to RF           

dur_trig = 1e-3;
delay_after_trig = 5.25e-3;         % sec delay after tigger

% due to rf dead time, there is an extra 0.1 ms of delay time prior to RF
% total duration between trigger's rising edge and start of RF pulse is
% then 1 ms + 5.25 ms + 0.1 ms = 6.35 ms -> needed to program TMS pulse

trig = mr.makeDigitalOutputPulse('ext1', 'duration', dur_trig);    % possible channels: 'osc0', 'osc1', 'ext1'

use_gz_excitation = 0;             % set to 0 to use ares gradient
T_gzReph = 5.25e-3;                % explicitly add this 5.25 ms Gz rephaser after slab select. 
                                   % when Gz is not played, this adds a 5 ms delay. otherwise it is a 5 ms rephasing lobe

nCyclesSpoil_acs = 4;              % try increasing spoiler moment for acs
scaleGrad_acs = 1;                 % set to 0 to remove all gradients during acs. set to 1: normal operation 
scaleGrad_acsSpoiler = 1;          % controls whether to use spoiler during gre acs -> set to 1 to use spoiler in acs

trig_delay_switch = 1;             % set this to 1 to add a final trigger and 500 ms delay after the entire sequence

tic
    [seq, sys, res, gro, seq_acs_ice] = writeEPI_bb_v10_acs_ice_tms_switch(voxelSize, mtx_size, TE, ...
        TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, nDummyFrames, encoding, filename, ...
        use_gre_acs, gMax, sMax, dwell2use, seq_path, rf_freqOffset, rfTB, rfDur, use_trig, delay_after_trig, trig, ...
        use_gz_excitation, T_gzReph, nCyclesSpoil_acs, scaleGrad_acs, scaleGrad_acsSpoiler, trig_delay_switch, ...
        'caipiPythonPath', caipiPythonPath, 'doRefScan', doRefscan);
toc

disp(['ESP: ', num2str(mr.calcDuration(gro)*1e3), ' ms'])



%--------------------------------------------------------------------------
%% seqeyes
%--------------------------------------------------------------------------

addpath /autofs/cluster/berkin/xingwang/share/seqeyes/

seqeyes(seq_acs_ice)

seqeyes(seq)


%--------------------------------------------------------------------------
%% evaluate ice labels
%--------------------------------------------------------------------------

adc_lbl = seq_acs_ice.evalLabels('evolution','adc');
hold on; plot(adc_lbl.LIN);
plot(adc_lbl.PAR) ;
legend('LIN','PAR');
title('evolution of labels/counters');


%--------------------------------------------------------------------------
%% check gradient spectrum
%--------------------------------------------------------------------------

addpath '/autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/utils'

% sys_type = 'prisma';
% sys_type = 'Connectome2';
sys_type = 'skyra';
% sys_type = 'terra';


freq_range4display = 3e3; % Hz

gradSpectrum_bb(seq,sys, sys_type, freq_range4display)

esp = mr.calcDuration(gro);

freq_main = 1 / (2 * esp);  % main epi frequency in Hz

disp(['esp: ', num2str(esp*1e3), ' main epi frequency: ', num2str(freq_main)])


%--------------------------------------------------------------------------
%% check pns
%--------------------------------------------------------------------------

addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/safe_pns_prediction-master/


if strcmp(sys_type,'prisma')
    [pns,tpns] = seq.calcPNS('/autofs/cluster/berkin/xingwang/Syncbb/others_toolboxes/pulseq/matlab/idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'); 
else
    if strcmp(sys_type,'skyra')
        [pns,tpns] = seq.calcPNS('/autofs/cluster/berkin/xingwang/Syncbb/others_toolboxes/pulseq/matlab/idea/asc/MP_GPA_K2309_2250V_793A_GC99.asc'); 
    else
        if strcmp(sys_type,'Connectome2')  
            [pns,tpns] = seq.calcPNS('/autofs/cluster/berkin/xingwang/Syncbb/others_toolboxes/pulseq/matlab/idea/asc/MP_GradSys_P034_c_CX600.asc'); 
        else
            if strcmp(sys_type,'terra')  
                % this is probably pre-upgrade
                % [pns,tpns] = seq.calcPNS('/autofs/cluster/berkin/xingwang/Syncbb/others_toolboxes/pulseq/matlab/idea/asc/MP_GPA_K2298_2250V_793A_SC72CD_EGA.asc'); 
                % is this the correct version?
                [pns,tpns] = seq.calcPNS('/autofs/cluster/berkin/xingwang/Syncbb/others_toolboxes/pulseq/matlab/idea/asc/MP_GradSys_K2298_2250V_1250A_W60_SC72CD.asc');                 
            else
                error('Undefined system')
            end
        end
    end
end


%--------------------------------------------------------------------------
%% generate single TR imaging data scan for loading the trajectory
%--------------------------------------------------------------------------

addpath /cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/

nDummyFrames = 0;
nFrames = 1;            % first frame is used for ghost correction if 'doRefScan' is true
use_gre_acs = 0;        % set to 0 to disable gre acs refscan


filename_traj = ['_traj_1mm_trajectory_RyxRz_', num2str(Ry), 'x', num2str(Rz), 'del', num2str(caipiShiftZ), '_pf', num2str(pf_ky)];

disp(filename_traj)

tic
    seq_traj = writeEPI_bb_v5_acs_ice(voxelSize, mtx_size, TE, TR, alpha, mb, pf_ky, Ry, Rz, caipiShiftZ, nFrames, ...
        nDummyFrames, encoding, filename, use_gre_acs, gMax, sMax, dwell2use, seq_path, rf_freqOffset, rfTB, rfDur, ...
        use_trig, delay_after_trig, trig, use_gz_excitation, T_gzReph, 'caipiPythonPath', caipiPythonPath, 'doRefScan', true);
toc

save([seq_path, 'seq', filename_traj, '.mat'], 'seq_traj')

% h = figure(10); close(h);
% seq_traj.plot();


%--------------------------------------------------------------------------
%% k-space trajectory calculation and plot
%--------------------------------------------------------------------------

[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

h = figure(1); close(h);

figure(1); plot(ktraj(2,:),ktraj(3,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(2,:),ktraj_adc(3,:),'r.'); % plot the sampling points
title([encoding, ' EPI, full k-space trajectory (k_y x k_z)']);


h = figure(2); close(h);
decimate_plot = 15;

figure(2); plot3(ktraj(1,1:decimate_plot:end), ktraj(2,1:decimate_plot:end),ktraj(3,1:decimate_plot:end),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot3(ktraj_adc(1,1:decimate_plot:end), ktraj_adc(2,1:decimate_plot:end),ktraj_adc(3,1:decimate_plot:end),'r.'); % plot the sampling points
title([encoding, ' EPI, full k-space trajectory (k_x x k_y x k_z)']);

%--------------------------------------------------------------------------
%% plot after acs
%--------------------------------------------------------------------------

TR_acs = 25e-3;

Ny_acs = 64;
Nz_acs = 64;

Ndummy = 50;            % dummy shots to enter steady state

time_acs = (Ny_acs * Nz_acs + Ndummy) * TR_acs

figure(100), seq.plot('timeRange', [time_acs,time_acs+.500])
figure(101), seq.plot('timeRange', [time_acs + TR*2,time_acs+TR*2+.500])


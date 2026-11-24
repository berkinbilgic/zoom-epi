%--------------------------------------------------------------------------
%% load data: GRE ACS scan and EPI imaging data combined
%--------------------------------------------------------------------------

% addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/'))      

restoredefaultpath
addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/LIBRARY'))

close all

data_path = '/autofs/space/daisy_001/users/berkin/2026_07_08_bay1_ares_phantom_zoom_rigidcoil/zoom/';
data_path_tsnr = '/autofs/space/daisy_001/users/berkin/2026_07_08_bay1_ares_phantom_zoom_rigidcoil/switch/';

use_tsnr = 0;


if use_tsnr
    % file_name_tsnr = 'meas_MID00712_FID82498_epi_zoom_noTMS_noSwitch_30runs';
    % file_name_tsnr = 'meas_MID00713_FID82499_epi_zoom_noTMS_yesSwitch_30runs';
    % file_name_tsnr = 'meas_MID00716_FID82502_epi_zoom_noTMS_noSwitch_30runs_exp2';
    file_name_tsnr = 'meas_MID00717_FID82503_epi_zoom_yesTMS_yesSwitch_30runs';

    dt_tsnr = mapVBVD2([data_path_tsnr, file_name_tsnr, '.dat']);
    
    data_tsnr = dt_tsnr{end}.image.unsorted();
    data_tsnr = permute(data_tsnr, [1,3,2]);
end

file_name_epi = 'meas_MID00706_FID82492_pulseq_v150_EPI_RF3kHz_6A_opt';


tic
    dt = mapVBVD2([data_path, file_name_epi, '.dat']);
toc

data_epi = dt{end}.image.unsorted();
data_epi = permute(data_epi, [1,3,2]);


Nx = s(data_epi,1);
num_chan = s(data_epi,3);


%--------------------------------------------------------------------------
% parse data: no gre acs in these .seq files -> acquired separately
%--------------------------------------------------------------------------

% num_acs = [48,48];
% 
% data_acs = data_epi(:,1:prod(num_acs),:);
% 
% kspace_acs = reshape(data_acs, [Nx, num_acs, num_chan]);
% 
% img_acs = ifft3call(kspace_acs);
% 
% imagesc3d2( rsos(img_acs,4), s(img_acs)/2, 1, [90,90,-90], [0,1e-4]), setGcf(.5)


% dt_acs = mapVBVD2([data_path, 'meas_MID00133_FID80384_plsq_v1p5_GRE_ACS_zm_RF3kHz_6A_opt.dat']);
% 
% data_acs = dt_acs{end}.image();
% data_acs = permute(data_acs, [1,3,4,2]);
% 
% img_acs = ifft3call(data_acs);
% 
% imagesc3d2( rsos(img_acs,4), s(img_acs)/2, 1, [90,90,-90], [0,1e-4]), setGcf(.5)
% 
% FOV = [220, 66, 64] * 1e-3;        % m -> zoomed FOV R2x2
% 
% genNiigz(flip(flip(flip(rsos(img_acs,4),1),2),3), FOV*1e3./size(img_acs(:,:,:,1)), [data_path, 'img_acs_rsos.nii'])        



%--------------------------------------------------------------------------
% parse data: epi
%--------------------------------------------------------------------------

% acquisition parameters
voxelSize = [1, 1, 1] * 1e-3;       % m
FOV = [220, 66, 64] * 1e-3;        % m -> zoomed FOV R2x2

mtx_size = round(FOV ./ voxelSize);
mtx_size = mtx_size + rem(mtx_size,2);

Ry = 1;
Rz = 1;

pf_ky = 1;
caipiShiftZ = 0 * (Rz>1);

% make sure mtx_size(2) is divisible by (Ry*2)
mtx_size_update = ceil(mtx_size(2) / Ry / 2) * Ry * 2;
mtx_size(2) = mtx_size_update;

if Rz == 3  
    % divisibility by Rz was not taken into account while designing the
    % sequence. address it here
    mtx_size(3) = ceil(mtx_size(3) / Rz) * Rz;
end

disp(num2str(mtx_size))

% data after acs scan
% data_epi_only = data_epi(:,1+prod(num_acs):end,:);  % no acs in this case
data_epi_only = data_epi;

% data without gy/gz blips for phase correction
data_epi_refscan = data_epi_only(:,1:end/2,:);
data_epi_refscan = reshape(data_epi_refscan, [Nx, mtx_size(2:3) ./ [Ry,Rz], num_chan]);

img_epi_refscan = ifft3call(data_epi_refscan);


% imaging data
data_epi_imaging = data_epi_only(:,1+end/2:end,:);
data_epi_imaging = reshape(data_epi_imaging, [Nx, mtx_size(2:3) ./ [Ry,Rz], num_chan]);

img_epi_imaging = ifft3call(data_epi_imaging);

imagesc3d2( rsos(img_epi_refscan,4), s(img_epi_imaging)/2, 3, [90,90,-90], [0,.1e-3]), setGcf(.5)
imagesc3d2( rsos(img_epi_imaging,4), s(img_epi_imaging)/2, 4, [90,90,-90], [0,.1e-3]), setGcf(.5)

if use_tsnr
    % tsnr data
    num_runs = 30;
    data_epi_tsnr = reshape(data_tsnr, [Nx, mtx_size(2:3) ./ [Ry,Rz], num_runs, num_chan]);
    data_epi_tsnr = permute(data_epi_tsnr, [1,2,3,5,4]);
    
    img_epi_tsnr = ifft3call(data_epi_tsnr);
    
    imagesc3d2( mean(rsos(img_epi_tsnr,4),5), s(img_epi_imaging)/2, 5, [90,90,-90], [0,.1e-3]), setGcf(.5)
end

%--------------------------------------------------------------------------
% SVD coil compression 
%--------------------------------------------------------------------------

use_svd = 0;

if use_svd
    num_chan = 20;
    
    [kspace_acs_svd, cmp_mtx] = svd_compress3d(kspace_acs, num_chan, 1);
    
    rmse_acs = rmse(rsos(kspace_acs_svd,4), rsos(kspace_acs,4))
    
    data_epi_imaging_svd = svd_apply3d(data_epi_imaging, cmp_mtx);
    data_epi_refscan_svd = svd_apply3d(data_epi_refscan, cmp_mtx);
    
    size_temp = size(data_epi_tsnr);
    size_temp(4) = num_chan;
    data_epi_tsnr_svd = zeross(size_temp);
    
    for t = 1:num_runs
        disp(t)
        data_epi_tsnr_svd(:,:,:,:,t) = svd_apply3d(data_epi_tsnr(:,:,:,:,t), cmp_mtx);
    end
    
    rmse_epi = rmse(rsos(data_epi_imaging_svd,4), rsos(data_epi_imaging,4))
    rmse_tsnr = rmse(rsos(data_epi_tsnr_svd,4), rsos(data_epi_tsnr,4))
    rmse_ref = rmse(rsos(data_epi_refscan_svd,4), rsos(data_epi_refscan,4))
else
    data_epi_imaging_svd = data_epi_imaging;
    data_epi_refscan_svd = data_epi_refscan;
    
    if use_tsnr
        data_epi_tsnr_svd = data_epi_tsnr;
    end
end


%--------------------------------------------------------------------------
% load trajectory for regridding
%--------------------------------------------------------------------------

addpath /cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ARES_zoom_fmri/pulseq-1.5.0/matlab/

traj_path = '/autofs/cluster/berkin/berkin/Matlab_Code_New/PULSEQ/ARES_zoom_fmri/07-Jul-2026/';

filename_traj = ['_traj_1mm_trajectory_RyxRz_', num2str(Ry), 'x', num2str(Rz), 'del', num2str(caipiShiftZ), '_pf', num2str(pf_ky)];


% sequence saved using a single TR, without ACS data
load([traj_path, 'seq', filename_traj, '.mat'])


traj_recon_delay = 0e-6;

tic
    [ktraj_adc_acq, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq_traj.calculateKspacePP('trajectory_delay', traj_recon_delay);
toc


%--------------------------------------------------------------------------
% apply ramp sampling correction
%--------------------------------------------------------------------------

ro_os = 1;
pF_fe = 1;

ktraj_use = reshape(ktraj_adc_acq, [3, Nx, mtx_size(2:3) ./ [Ry,Rz]]);

kxmin   = min(min(sq(ktraj_use(1,:,:,1))))*(Nx*pF_fe*ro_os/2-1)/(Nx*pF_fe*ro_os/2);
kxmax   = max(max(sq(ktraj_use(1,:,:,1))))/(Nx*pF_fe*ro_os/2-1)*(Nx*pF_fe*ro_os/2);
kxx     = ((-Nx*pF_fe*ro_os/2):(Nx*pF_fe*ro_os/2-1))/(Nx*pF_fe*ro_os/2);

if floor(Nx*pF_fe*ro_os/2) == Nx*pF_fe*ro_os/2
    kxx     = (kxx+1)/2*(kxmax-kxmin)+kxmin;
else
    kxx     = (kxx-min(kxx))/(max(kxx)-min(kxx))*(kxmax-kxmin)+kxmin;
end

h = figure(200), close(h)
figure(200), plot(kxx, sq(ktraj_use(1,:,1,1)))
figure(200), hold on, plot(kxx, ktraj_use(1,:,2,1), 'r')


data_epi_regrid = zeross(size(data_epi_imaging_svd));
data_ref_regrid = zeross(size(data_epi_refscan_svd));

tic
for cc = 1:num_chan
    disp(['channel: ', num2str(cc)])

    for aa = 1:mtx_size(2) / Ry
        for bb = 1:mtx_size(3) / Rz

            if mod(aa,2) == 0
                % correct the flip between even and odd lines
                data_epi_regrid(:,aa,bb,cc) = interp1(ktraj_use(1,:,1,1), flip(data_epi_imaging_svd(:,aa,bb,cc),1), kxx, 'spline', 0);
                data_ref_regrid(:,aa,bb,cc) = interp1(ktraj_use(1,:,1,1), flip(data_epi_refscan_svd(:,aa,bb,cc),1), kxx, 'spline', 0);
            else
                data_epi_regrid(:,aa,bb,cc) = interp1(ktraj_use(1,:,1,1), data_epi_imaging_svd(:,aa,bb,cc), kxx, 'spline', 0);
                data_ref_regrid(:,aa,bb,cc) = interp1(ktraj_use(1,:,1,1), data_epi_refscan_svd(:,aa,bb,cc), kxx, 'spline', 0);
            end   
        end
    end
end
toc

if use_tsnr
    data_epi_regrid_tsnr = zeross(size(data_epi_tsnr_svd));

    tic
    for t = 1:num_runs
        disp(['run: ', num2str(t)])
    
        for cc = 1:num_chan
            for aa = 1:mtx_size(2) / Ry
                for bb = 1:mtx_size(3) / Rz
        
                    if mod(aa,2) == 0
                        % correct the flip between even and odd lines
                        data_epi_regrid_tsnr(:,aa,bb,cc,t) = interp1(ktraj_use(1,:,1,1), flip(data_epi_tsnr_svd(:,aa,bb,cc,t),1), kxx, 'spline', 0);
                    else
                        data_epi_regrid_tsnr(:,aa,bb,cc,t) = interp1(ktraj_use(1,:,1,1), data_epi_tsnr_svd(:,aa,bb,cc,t), kxx, 'spline', 0);
                    end   
                end
            end
        end
    end
    toc

    img_epi_regrid_tsnr = ifft3call(data_epi_regrid_tsnr);

    imagesc3d2( mean(rsos(img_epi_regrid_tsnr,4),5), s(img_epi_imaging)/2, 14, [90,90,-90], [0,1.e-4], [], 'regridded epi tsnr mean'), setGcf(.5)
end

img_epi_regrid = ifft3call(data_epi_regrid);

% zero pad acs to epi's matrix size
% img_acs_pad = ifft3call(padarray(kspace_acs_svd, [0, mtx_size(2:3)-num_acs]/2));

% imagesc3d2( rsos(img_acs_pad,4), s(img_acs_pad)/2, 11, [90,90,-90], [0,1e-4], [], 'gre acs pad'), setGcf(.5)
imagesc3d2( rsos(img_epi_regrid,4), s(img_epi_imaging)/2, 12, [90,90,-90], [0,1.e-4], [], 'regridded epi'), setGcf(.5)
imagesc3d2( rsos(img_epi_imaging,4), s(img_epi_imaging)/2, 13, [90,90,-90], [0,1.e-4], [], 'before regridding'), setGcf(.5)


%--------------------------------------------------------------------------
% k-space based ghost correction
%--------------------------------------------------------------------------

addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/BUDA/fleet_buda
addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/BUDA/fleet_buda/mrir_toolbox

%   dimension 08: Segments
%   dimension 10: Partitions
AccY = 1;

PAT_nav = zeros([Nx, mtx_size(2:3) ./ [Ry,Rz], num_chan, 2]);
PAT_nav(:,1:2:end,:,:,1) = data_ref_regrid(:,1:2:end,:,:); 
PAT_nav(:,2:2:end,:,:,2) = data_ref_regrid(:,2:2:end,:,:); 
PAT_nav = permute(PAT_nav, [1,2,4,6,7,8,9,5,10,3]);

PAT_ref = zeros([size(data_epi_regrid),2]);
PAT_ref(:,1:2:end,:,:,1) = data_epi_regrid(:,1:2:end,:,:); 
PAT_ref(:,2:2:end,:,:,2) = data_epi_regrid(:,2:2:end,:,:); 
PAT_ref = permute(PAT_ref, [1,2,4,6,7,8,9,5,10,3]);

if use_tsnr
    PAT_ref_tsnr = zeros([size(data_epi_regrid),2,num_runs]);
    PAT_ref_tsnr(:,1:2:end,:,:,1,:) = data_epi_regrid_tsnr(:,1:2:end,:,:,:); 
    PAT_ref_tsnr(:,2:2:end,:,:,2,:) = data_epi_regrid_tsnr(:,2:2:end,:,:,:); 
end
% ghost correction function without regridding:
tic
    navCor_kspace = ghost_correct_pat_ref_v1_STD_bb_v1(PAT_ref, PAT_nav, AccY);
toc

navCor_img = sq(ifftc(ifft2call(navCor_kspace), 10));


if use_tsnr
    navCor_kspace_tsnr = zeross([size(navCor_kspace), num_runs]);
    
    for t = 1:num_runs
        disp(t)
        PAT_ref_tsnr2use = permute(PAT_ref_tsnr(:,:,:,:,:,t), [1,2,4,6,7,8,9,5,10,3]);
        navCor_kspace_tsnr2use = ghost_correct_pat_ref_v1_STD_bb_v1(PAT_ref_tsnr2use, PAT_nav, AccY);
        navCor_kspace_tsnr(:,:,:,:,:,:,:,:,:,:,t) = navCor_kspace_tsnr2use;
    end
    
    navCor_img_tsnr = sq(ifftc(ifft2call(navCor_kspace_tsnr), 10));

    imagesc3d2( mean(rsos(navCor_img_tsnr,3),5), s(img_epi_imaging)/2, 15, [90,90,-90], [0,1.5e-4], [], 'tsnr mean ghost corr and regridding'), setGcf(.5)
end

imagesc3d2( rsos(navCor_img,3), s(img_epi_imaging)/2, 14, [90,90,-90], [0,1.5e-4], [], 'ghost corr and regridding'), setGcf(.5)

% save(strcat(data_path, 'img_regrid_lpc_', file_name_epi, '.mat'), 'navCor_img', '-V7.3');


%--------------------------------------------------------------------------
% save as nifti
%--------------------------------------------------------------------------

addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/LIBRARY'))

navCor_img_tsnr_rsos = sq(rsos(navCor_img_tsnr,3));

for t = 1:5:num_runs
    imagesc3d2( navCor_img_tsnr_rsos(:,:,:,t), s(navCor_img_tsnr_rsos)/2, t, [90,90,-90], [0,2.e-4], [], 'tsnr mean ghost corr and regridding'), setGcf(.5)
end

imagesc3d2( rsos(navCor_img,3), s(navCor_img)/2, 11, [90,90,-90], [0,2.e-4], [], 'ghost corr and regridding'), setGcf(.5)


% genNiigz(sq(rsos(navCor_img,3)), [1,1,1], [data_path, 'img_rsos_1TR.nii'])        
% genNiigz(navCor_img_tsnr_rsos, [1,1,1], [data_path_tsnr, 'img_rsos_30TRs_', file_name_tsnr, '.nii'])        



%--------------------------------------------------------------------------
%% espirit coil sensitivity estimation
%--------------------------------------------------------------------------
        
addpath '/autofs/cluster/berkin/berkin/Matlab_Code_New/NEATR/NEATR_Wip_936_Toolbox/neatr_wip_936_sense'

num_acs = 32;
kernel_size = [6,6];
eigen_thresh = 0.7;

sens = zeross(size(img_acs_pad));

delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

total_cores = c.NumWorkers; 
parpool(ceil(total_cores/2))

tic
parfor slc_select = 1:s(img_acs_pad,3)     
    disp(num2str(slc_select))

    [maps, weights] = ecalib_soft( fft2c( sq(img_acs_pad(:,:,slc_select,:)) ), num_acs, kernel_size, eigen_thresh );

    sens(:,:,slc_select,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
end 
toc

delete(gcp('nocreate'))

save(strcat(data_path, 'sens_', file_name_epi, '.mat'), 'sens', '-V7.3');

sens_epi = permute(sens, [1,2,4,3]);

    
%--------------------------------------------------------------------------
%% load coil sensitivity maps
%--------------------------------------------------------------------------

tic
    load(strcat(data_path, 'sens_', file_name_epi, '.mat'));
toc

sens_epi = permute(sens, [1,2,4,3]);


%--------------------------------------------------------------------------
%% zero fill image
%--------------------------------------------------------------------------

kspace_zf = zeros([Nx, mtx_size(2:3), num_chan]);

temp = permute(sq(navCor_kspace), [1,2,4,3]);       % using k-space LPC       

kspace_zf(:,1:2*Ry:end,1:Rz:end,:) = temp(:,1:2:end,:,:);
kspace_zf(:,1+Ry:2*Ry:end,1:Rz:end,:) = temp(:,2:2:end,:,:);

img_zf = ifft3call(kspace_zf);

imagesc3d2( coil_combine(img_zf,sens,4), s(img_zf)/2, 15, [90,90,-90], [0,1e-4], [], 'coil combo zero fill'), setGcf(.5)


kspace_zf_tsnr = zeros([Nx, mtx_size(2:3), num_chan, num_runs]);

temp = permute(sq(navCor_kspace_tsnr), [1,2,4,3,5]);       % using k-space LPC       

kspace_zf_tsnr(:,1:2*Ry:end,1:Rz:end,:,:) = temp(:,1:2:end,:,:,:);
kspace_zf_tsnr(:,1+Ry:2*Ry:end,1:Rz:end,:,:) = temp(:,2:2:end,:,:,:);

img_zf_tsnr = ifft3call(kspace_zf_tsnr);

imagesc3d2( coil_combine(mean(img_zf_tsnr,5),sens,4), s(img_zf)/2, 16, [90,90,-90], [0,1e-4], [], 'tsnr coil combo zero fill'), setGcf(.5)



%--------------------------------------------------------------------------
%% LSQR recon
%--------------------------------------------------------------------------

addpath /autofs/cluster/berkin/berkin/Matlab_Code/TOOLBOXES/SENSE_LSQR_Toolbox

lsqr_iter = 100;
lsqr_tol = 1e-2;

param.m3d = kspace_zf~=0;
param.sens = sens;
param.N = [Nx, mtx_size(2:3)];
param.num_chan = num_chan;
param.lambda = 1e-3;        % L2 reg

tic
    res = lsqr(@apply_sense_tikc3d, cat(1, kspace_zf(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  
toc

Res = reshape(res, param.N);

img_acs_pad_rsos = rsos(img_acs_pad,4);


Res_tsnr = zeross([param.N,num_runs]);

for t = 1:num_runs
    disp(t)

    kspace_zf2use = kspace_zf_tsnr(:,:,:,:,t);
    res = lsqr(@apply_sense_tikc3d, cat(1, kspace_zf2use(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  

    Res_tsnr(:,:,:,t) = reshape(res, param.N);
end


Res_combo = zeross([param.N,num_runs]);

for t = 1:num_runs
    disp(t)

    Res_combo(:,:,:,t) = coil_combine(img_zf_tsnr(:,:,:,:,t), sens, 4);
end


% save(strcat(data_path, 'img_lsqr_lpc_', file_name_epi, '.mat'), 'Res', '-V7.3');
% save(strcat(data_path, 'img_acs_pad_rsos_', file_name_epi, '.mat'), 'img_acs_pad_rsos', '-V7.3');
save(strcat(data_path, 'img_lsqr_lpc_tsnr_30runs_', file_name_tsnr , '.mat'), 'Res_tsnr', '-V7.3');
save(strcat(data_path, 'img_coil_combo_lpc_tsnr_30runs_', file_name_tsnr , '.mat'), 'Res_combo', '-V7.3');


imagesc3d2( mean(Res_combo,4) / sqrt(Rz), [70,16,16], 4, [90,90,-90], [0,1e-4], [], 'coil combo mean tsnr'), setGcf(.5)
imagesc3d2( mean(Res_tsnr,4) / sqrt(Rz), [70,16,16], 5, [90,90,-90], [0,1e-4], [], 'mean tsnr'), setGcf(.5)
imagesc3d2( Res / sqrt(Rz), [70,16,16], 6, [90,90,-90], [0,1e-4], [], 'lsqr'), setGcf(.5)
imagesc3d2( img_acs_pad_rsos, [70,16,16], 7, [90,90,-90], [0,1e-4], [], 'gre acs'), setGcf(.5)


genNiigz(abs(Res_combo), [1,1,1], strcat(data_path, 'img_coil_combo_lpc_tsnr_30runs_', file_name_tsnr, '.nii'))
genNiigz(abs(Res_tsnr), [1,1,1], strcat(data_path, 'img_lsqr_lpc_tsnr_30runs_', file_name_tsnr, '.nii'))
% genNiigz(abs(Res), [1,1,1], strcat(data_path, 'img_lsqr_lpc_', file_name_epi, '.nii'))
% genNiigz(abs(img_acs_pad_rsos), [1,1,1], strcat(data_path, 'img_acs_pad_rsos_', file_name_epi, '.nii'))


 
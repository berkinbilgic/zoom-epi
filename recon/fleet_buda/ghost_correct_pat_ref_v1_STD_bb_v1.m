function [navCor_kspace] = ghost_correct_pat_ref_v1_STD_bb_v1(PAT_ref, PAT_nav, AccY)
% Ghost correction for EPI
% does not perform regridding, thus does not need "prot" header information

disp('Step:  Ghost Correction start')

% estimation
lin_fit_navCor = mrir_artifact_ghost_compute_CorrelationMethod_bb(PAT_nav);

sz_patrefscan = size(PAT_ref);
nPhases = sz_patrefscan(2);
sy = 1; % start line

Ref=zeros(size(PAT_ref));
Ref_navCor=zeros(size(PAT_ref));

for R = 0 : AccY-1
    disp(num2str(R))

    Ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:) = PAT_ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:);

    Ref_navCor(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:)= mrir_artifact_ghost_correct_CorrelationMethod_v3(Ref(:,sy+R:AccY:nPhases,:,1,1,1,1,:,1,:), lin_fit_navCor);
end

data_hybrid = mrir_iDFT_freqencode(Ref_navCor); 

clear Ref_navCor

Ref_navCor_kspace = mrir_fDFT_freqencode(data_hybrid);

navCor_kspace = sum(Ref_navCor_kspace,8);

disp('Step:  Ghost Correction complete')

end
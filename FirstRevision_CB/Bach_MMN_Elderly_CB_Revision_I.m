
%%

%% BACH ELDERLY MEG - COMMUNICATIONS BIOLOGY - REVISION I

%%

%% Defining the conditions - only for old/new paradigm (in the case of MMN the conditions were defined during epoching)

definecond_l = 1; % 1 = defining conditions; 0 = computing descriptive statistics

xlsx_dir_behav = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/MEG_behavioural'; %dir to MEG behavioral results (.xlsx files)
epoch_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/esp*cog*.mat'); %dir to epoched files
clear Block_3
for ii = 1:length(epoch_list) %over epoched data
    if ~strcmp('mmn',epoch_list(ii).name(18:20)) %if the block is not the MMN block..
        %loading SPM object (epoched data)
        D = spm_eeg_load([epoch_list(ii).folder '/' epoch_list(ii).name]);
        dummy = D.fname; %%% OBS!! (SUBJ0041 IS NOT AUD1 AS IN THE NAME BUT IT IS AUD2!!) HERE YOU DID MANUAL ADJUSTMENT OF VARIABLE dummy TO FIX IT..
        %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
        if strcmp(dummy(18:20),'rec') || strcmp(dummy(18:20),'rac')
            dumbloc = 'Block_3.xlsx';
            bl = 3;
        end
        dumls = ['Subj_' dummy(13:16) '_' dumbloc];
        [~,~,raw_recog] = xlsread([xlsx_dir_behav '/' dumls]); %excel files
        %picking the current block
        if bl == 3 %block 3
            if definecond_l == 1
                for k = 1:length(D.trialonset)
                    if raw_recog{(k + 1),3} == 0 %if there was no response
                        if strcmp(raw_recog{(k + 1),2}(8:9),'ol')
                            D = D.conditions(k,'No_response_Old');
                        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1')
                            D = D.conditions(k,'No_response_New_T1');
                        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3')
                            D = D.conditions(k,'No_response_New_T3');
                        end
                    elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
                        D = D.conditions(k,'Old_Correct'); %assign old correct
                    elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 2 %old incorrect
                        D = D.conditions(k,'Old_Incorrect'); %otherwise assign new correct
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
                        if contains(raw_recog{(k + 1),2},{'m1t1e1','m1t1e6','m2t1e1','m2t1e6','m3t1e6'})
                            D = D.conditions(k,'New_T1_Corr_OutK');
                        else
                            D = D.conditions(k,'New_T1_Corr_InK');
                        end
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 1 %new t1 incorrect
                        D = D.conditions(k,'New_T1_Incorrect');
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
                        if contains(raw_recog{(k + 1),2},{'m1t3e1','m2t3e1','m2t3e6','m3t3e1','m3t3e7'})
                            D = D.conditions(k,'New_T3_Corr_OutK');
                        else
                            D = D.conditions(k,'New_T3_Corr_InK');
                        end
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 1 %new t3 incorrect
                        D = D.conditions(k,'New_T3_Incorrect');
                    end
                end
            else %computing descriptive statistics
                %legend
                Block_3{1,1} = 'Subject'; Block_3{1,2} = 'OLD_Cor'; Block_3{1,3} = 'OLD_Cor %'; Block_3{1,4} = 'New_T1_Cor'; Block_3{1,5} = 'New_T1_Cor %'; %1st row
                Block_3{1,6} = 'New_T1_Incor'; Block_3{1,7} = 'New_T1_Incor %'; Block_3{1,8} = 'New_T3_Cor'; Block_3{1,9} = 'New_T3_Cor %'; Block_3{1,10} = 'New_T3_Incor'; Block_3{1,11} = 'New_T3_Incor %'; Block_3{1,12} = 'No response'; Block_3{1,13} = 'No response %'; %1st row
                Block_3{1,14} = 'OLD_RT'; Block_3{1,15} = 'NEWT1_RT'; Block_3{1,16} = 'NEWT2_RT'; Block_3{1,17} = 'NEWT3_RT'; Block_3{1,18} = 'NEWT4_RT'; Block_3{1,19} = 'Old_Incor'; Block_3{1,20} = 'Old_Incor %';
                nr = 0; old = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0; oldi = 0;
                ort = []; n1rt = []; n2rt = []; n3rt = []; n4rt = []; orti = [];
                for k = 1:81
                    if raw_recog{(k + 1),3} == 0 %if there was no response
                        nr = nr + 1;
                    elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
                        old = old + 1;
                        ort = cat(1,ort,raw_recog{(k + 1),4});
                    elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 2 %old incorrect
                        oldi = oldi + 1;
                        orti = cat(1,orti,raw_recog{(k + 1),4});
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
                        n1 = n1 + 1;
                        n1rt = cat(1,n1rt,raw_recog{(k + 1),4});
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 1 %new t1 incorrect
                        n2 = n2 + 1;
                        n2rt = cat(1,n2rt,raw_recog{(k + 1),4});
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
                        n3 = n3 + 1;
                        n3rt = cat(1,n3rt,raw_recog{(k + 1),4});
                    elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 1 %new t3 incorrect
                        n4 = n4 + 1;
                        n4rt = cat(1,n4rt,raw_recog{(k + 1),4});
                    end
                end
                disp(num2str(['Subject ' num2str(ii)]))
                Block_3{ii+1,1} = dummy(13:16); Block_3{ii+1,2} = old; Block_3{ii+1,3} = (old/27)*100; Block_3{ii+1,4} = n1; Block_3{ii+1,5} = (n1/27)*100;
                Block_3{ii+1,6} = n2; Block_3{ii+1,7} = (n2/27)*100; Block_3{ii+1,8} = n3; Block_3{ii+1,9} = (n3/27)*100; Block_3{ii+1,10} = n4; Block_3{ii+1,11} = (n4/27)*100; Block_3{ii+1,12} = nr; Block_3{ii+1,13} = (nr/135)*100;
                Block_3{ii+1,14} = mean(ort); Block_3{ii+1,15} = mean(n1rt); Block_3{ii+1,16} = mean(n2rt); Block_3{ii+1,17} = mean(n3rt); Block_3{ii+1,18} = mean(n4rt); Block_3{ii+1,19} = oldi; Block_3{ii+1,20} = mean(orti);
            end
        end
        %this is for every block
        if ~isempty(D.badtrials) %overwriting badtrials (if any) on condition labels
            BadTrials = D.badtrials;
            for badcount = 1:length(BadTrials) %over bad trials
                D = D.conditions(BadTrials(badcount),'Bad_trial');
            end
        end
        D = D.montage('switch',1);
        D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
        D.save(); %saving data on disk
    end
    disp(num2str(ii))
end

%%

%%

%% *** MEG SOURCE RECONSTRUCTION - BEAMFORMING ***

%%

%% CREATING 8mm PARCELLATION FOR EASIER INSPECTION IN FSLEYES
%OBS!! This section is done only for better handling of some visualization purposes, but it does not affect any of the beamforming algorithm;
% it is just important not to mix up the MNI coordinates, thus I would recommend to use the following lines

%%% OBS!! THIS IS NOT NEEDED SINCE IT IS EXACTLY THE SAME AS IN THE ORIGINAL LEADING SCRIPT %%%

%% CONVERSION T1 - DICOM TO NIFTI

%%% OBS!! THIS IS NOT NEEDED SINCE IT IS EXACTLY THE SAME AS IN THE ORIGINAL LEADING SCRIPT %%%

%% RHINO coregistration

%%% OBS!! THIS IS NOT NEEDED SINCE IT IS EXACTLY THE SAME AS IN THE ORIGINAL LEADING SCRIPT %%%

%% BEAMFORMING

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% FUNCTION FOR SOURCE RECONSTRUCTION

%user settings
clust_l = 1; %1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
timek = 1:1026; %time-points
freqq = []; %frequency range (empty [] for broad band)
% freqq = [0.1 1]; %frequency range (empty [] for broad band)
% freqq = [2 8]; %frequency range (empty [] for broad band)
sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
workingdir2 = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD'; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
invers = 1; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction
absl = 0; % 1 = absolute value of sources; 0 = not

%actual computation
list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*cogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
condss = {'Old_Correct','New_T1_Corr_OutK','New_T1_Corr_InK','New_T3_Corr_OutK','New_T3_Corr_InK'};
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat');
extr = 13:16;
if isempty(freqq)
    workingdir = [workingdir2 '/Block_3/CB_Rev_I_Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers)];
else
    workingdir = [workingdir2 '/Block_3/CB_Rev_I_Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers)];
end
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
if ~exist(workingdir,'dir') %creating working folder if it does not exist
    mkdir(workingdir)
end
for ii = 2:length(list) %over subjects
    extr2 = extr;
    S = [];
    if ~isempty(freqq) %if you want to apply the bandpass filter, you need to provide continuous data
        %and assigning the path to the structure S
        S.norm_megsensors.MEGdata_c = [list(ii).folder '/' list(ii).name(2:end)];
    end    
    S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation
    
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir;
    S.norm_megsensors.MEGdata_e = [list(ii).folder '/' list(ii).name];
    S.norm_megsensors.freq = freqq; %frequency range
    S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')
    
    S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
    S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')
    
    S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
    %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
    %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
    %
    S.inversion.timef = timek; %data-points to be extracted (e.g. 1:300); leave it empty [] for working on the full length of the epoch
    S.inversion.conditions = condss; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
    S.inversion.bc = [1 26]; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
    S.inversion.effects = invers;
    
    S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = [1 2030 3269]; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.
    
    S.nifti = 1; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(ii).name(extr2)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    
    if clust_l ~= 1 %useful  mainly for begugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end

%%

%% *** GETTING NEW ROIs FROM AAL AND FROM FUNCTIONAL ROIs

%% FROM 3D NIFTI IMAGE (OR ALL ROIs OR MNI COORDINATES) TO LBPD COORDINATES
%not needed here but useful for future references

ROIs_l = 2; %1 = AAL ROIs; 2 = systematic variation recognition paper ROIs

if ROIs_l == 1
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/FinalROIs/*.gz');
else
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/FinalROIs/ROIs_SystVarPaper_PlusIFG/*.gz');
end
idx_LBPD = cell(length(list),1);
for ii = 1:length(list)
    S = [];
    S.input = 3; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
    S.coordd = []; %coordinates in MNI space (x,y,z)
    S.AAL_ROIs = []; %AAL ROIs numbers you want to use
    % S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
    S.image = [list(ii).folder '/' list(ii).name];
    %actual function
    idx_LBPD{ii} = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);
end

%% testing that everyting with the previous function was fine after updates

% vect = zeros(3559,1);
% vect(idx_LBPD{1}(:,1)) = 1;
% S.data = vect; %data (voxels x ROIs (time-points))
% S.fname = '//scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/FinalROIs'; %path and name for the image to be saved
% names{1} = 'LAC2';
% S.names = names; %names for the different images
% %actual function
% FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);

%% PREPARING TIME SERIES FOR EACH ROI AND SAVING THEM

ROIs_l = 2; %1 = AAL ROIs; 2 = systematic variation recognition paper ROIs

%list of source reconstructed data (for each subject)
%this is the source reconstruction with the in-key and out-of-key "new" melodies 
list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/CB_Rev_I_Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
%this is from previous source reconstruction ("original" labels) since the computation of this particular piece of information is the same 
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');

timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
% ROIII = {1,2,3,4,5,6}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
load([list(1).folder '/' list(1).name])

% list(36:end) = [];
% list(1:34) = [];
dum2 = zeros(length(idx_LBPD),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3),length(list));
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    load([list(ii).folder '/' list(ii).name])
    dum = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3));
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(OUT.sources_ERFs,1) %over brain voxels
            dum(jj,:,cc) = OUT.sources_ERFs(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            disp(['subject - ' num2str(ii) ' - condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for cc = 1:size(dum,3) %over conditions
        for pp = 1:length(idx_LBPD) %over ROIs (with possibility of averaging together voxels of different ROIs (e.g. left and right auditory cortex (AC))
            dum2(pp,:,cc,ii) = mean(dum(idx_LBPD{pp}(:,1),:,cc),1);
        end
    end
end
condds = OUT.S.inversion.conditions;
if ROIs_l == 1
    save(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8.mat'],'dum2','condds')
else
    save(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/CB_Rev_I_Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_8_SystPaper.mat'],'dum2','condds')
end

%SUBJ0035 HAS NO DATA FOR NEWT1 AND NEWT3!!

%%

%% HERE PLOTTING (AND STATISTICS) FOR DIFFERENT "NEW" MELODIES IN RELATION TO THE MUSICAL KEY

%preparing information for the anova (it counts for both NewT1 and NewT3)
clear gsubj
gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78]; %young
gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76]; %elderly
cnt = 0;
vect = [1 1 2 2];
AGE = cell((length(gsubj{1}) + length(gsubj{2}))*2,1);
COND = cell((length(gsubj{1}) + length(gsubj{2}))*2,1);
for ii = 1:4 %over groups x conditions (2 x 2)
    xx = vect(ii); %getting index for groups
    for pp = 1:length(gsubj{xx}) %over young
        cnt = cnt + 1;
        if ii < 3 %if young group
            AGE{cnt} = 'y';
        else
            AGE{cnt} = 'e';
        end
        if mod(ii,2) == 1 %if ii is odd
            COND{cnt} = 'OutK'; %then you have the out of key condition
        else
            COND{cnt} = 'InK'; %otheriwse the in key condition
        end
    end
end
%statistics (ANOVAs)
%NewT1
P = zeros(size(dum2,1),size(dum2,2),3);
F = zeros(size(dum2,1),size(dum2,2),3);
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %over time-points
        datadum = [];
        for gg = 1:length(gsubj) %over groups
            for cc = 1:2 %over pair of experimental conditions (e.g. in-key and out-of-key NewT1
                datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc+1,gsubj{gg}))); %concatenating data for the anova
            end
        end
        %ANOVA
        %two-way ANOVA (anovan Matlab function)
        [p,f,stats] = anovan(datadum,{AGE,COND},'model','interaction','varnames',{'age','cond'},'display','off');
        P(ii,jj,1) = p(1); %age (main effect) - pvalue
        F(ii,jj,1) = f{2,6}; %age (main effect) - Fvalue
        P(ii,jj,2) = p(2); %cond (main effect) - pvalue
        F(ii,jj,2) = f{3,6}; %cond (main effect) - Fvalue
        P(ii,jj,3) = p(3); %age x cond (interaction effect) - pvalue
        F(ii,jj,3) = f{4,6}; %age x cond (interaction effect) - Fvalue
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
EFF{1} = 'Age'; EFF{2} = 'Condition'; EFF{3} = 'Interaction';
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/CB_Rev_I_Beam_abs_0_sens_1_freq_broadband_invers_1/Stats_Keys_Melodies';
mkdir(outdir)
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:size(F,3) %over effects of the anova
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        F_sel = F(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), F_sel ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        save([outdir '/NewT1' ROIN{ii} '_' EFF{cc} '.mat'],'PDn'); %printing excel file
        writetable(PDn,[outdir '/NewT1' ROIN{ii} '_' EFF{cc} '.xlsx'],'Sheet',1); %printing excel file
    end
end
%NewT3
P = zeros(size(dum2,1),size(dum2,2),3);
F = zeros(size(dum2,1),size(dum2,2),3);
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %over time-points
        datadum = [];
        for gg = 1:length(gsubj) %over groups
            for cc = 1:2 %over pair of experimental conditions (e.g. in-key and out-of-key NewT1
                datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc+3,gsubj{gg}))); %concatenating data for the anova
            end
        end
        %ANOVA
        %two-way ANOVA (anovan Matlab function)
        [p,f,stats] = anovan(datadum,{AGE,COND},'model','interaction','varnames',{'age','cond'},'display','off');
        P(ii,jj,1) = p(1); %age (main effect) - pvalue
        F(ii,jj,1) = f{2,6}; %age (main effect) - Fvalue
        P(ii,jj,2) = p(2); %cond (main effect) - pvalue
        F(ii,jj,2) = f{3,6}; %cond (main effect) - Fvalue
        P(ii,jj,3) = p(3); %age x cond (interaction effect) - pvalue
        F(ii,jj,3) = f{4,6}; %age x cond (interaction effect) - Fvalue
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
EFF{1} = 'Age'; EFF{2} = 'Condition'; EFF{3} = 'Interaction';
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/CB_Rev_I_Beam_abs_0_sens_1_freq_broadband_invers_1/Stats_Keys_Melodies';
mkdir(outdir)
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:size(F,3) %over effects of the anova
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        F_sel = F(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), F_sel ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        save([outdir '/NewT3' ROIN{ii} '_' EFF{cc} '.mat'],'PDn'); %printing excel file
        writetable(PDn,[outdir '/NewT3' ROIN{ii} '_' EFF{cc} '.xlsx'],'Sheet',1); %printing excel file
    end
end

%% plotting (Aarhus' server)

load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/CB_Rev_I_Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_8_SystPaper.mat')

clear conds
gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78]; %young
gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76]; %elderly
GGG{1} = 'yng'; GGG{2} = 'eld';
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
conds{1} = ' Mem '; conds{2} = ' NewT1 OutK '; conds{3} = ' NewT1 InK '; conds{4} = ' NewT3 OutK '; conds{5} = ' NewT3 InK ';
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
for ii = 1:size(dum2,1)%over clusters
    %newT1
    figure
    cnt = 0;
    for cc = 2:3%1:size(dum2,3) %over conditions
        for gg = 1:length(gsubj)
            cnt = cnt + 1;
            plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,cc,gsubj{gg}),4)),'Color',color_line2(cnt,:),'LineWidth',2,'DisplayName',[conds{cc} ' ' GGG{gg}])
            hold on
            plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,cc,gsubj{gg}),4)) + (nanstd(dum2(ii,:,cc,gsubj{gg}),0,4)./sqrt(length(gsubj{gg}))),':','Color',color_line2(cnt,:),'LineWidth',0.5,'HandleVisibility','off')
            hold on
            plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,cc,gsubj{gg}),4)) - (nanstd(dum2(ii,:,cc,gsubj{gg}),0,4)./sqrt(length(gsubj{gg}))),':','Color',color_line2(cnt,:),'LineWidth',0.5,'HandleVisibility','off')
            hold on
        end
    end
    grid minor
    title(ROIN{ii})
    xlim([-0.1 3.4])
    %     ylim(limmy)
    set(gcf,'color','w')
    legend('show')
    %newT3
    figure
    cnt = 0;
    for cc = 4:5%1:size(dum2,3) %over conditions
        for gg = 1:length(gsubj)
            cnt = cnt + 1;
            plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,cc,gsubj{gg}),4)),'Color',color_line2(cnt,:),'LineWidth',2,'DisplayName',[conds{cc} ' ' GGG{gg}])
            hold on
            plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,cc,gsubj{gg}),4)) + (nanstd(dum2(ii,:,cc,gsubj{gg}),0,4)./sqrt(length(gsubj{gg}))),':','Color',color_line2(cnt,:),'LineWidth',0.5,'HandleVisibility','off')
            hold on
            plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,cc,gsubj{gg}),4)) - (nanstd(dum2(ii,:,cc,gsubj{gg}),0,4)./sqrt(length(gsubj{gg}))),':','Color',color_line2(cnt,:),'LineWidth',0.5,'HandleVisibility','off')
            hold on
        end
    end
    grid minor
    title(ROIN{ii})
    xlim([-0.1 3.4])
    %     ylim(limmy)
    set(gcf,'color','w')
    legend('show')
end

%%

%%

%% HERE COMPUTE STATISTICS BETWEEN ORIGINAL CONDITIONS AND FAMILIARITY, ETC. WITH THE MUSICAL PIECE

%SUBJ0035 NO BRAIN DATA!!

%behavioural measures: 4 = familiarity with the Bach's musical piece; 5 = liking the piece; 6 = emotional; 7 = arousal

[fam,~,fam_raw] = xlsread('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Papers/Bach_Elderly_CommunicationsBiology_Revision_I/Familiarity_CB_Revision.xlsx'); %loads working memory data
%this is the file that you have to load and use to run the statistics
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat');
clear ROIN
ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
%t-tests
P = zeros(size(dum2,1),size(dum2,2),size(dum2,3),4); %ROIs x time-points x conditions x 4 correlations
R = zeros(size(dum2,1),size(dum2,2),size(dum2,3),4); %ROIs x time-points x conditions x 4 correlations
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %over time-points
        for cc = 1:size(dum2,3) %over correlations to be performed
            for ee = 1:size(P,4)
                %codes t-test
                a = squeeze(dum2(ii,jj,cc,:)); %neural data
                b = fam(:,ee+3); %behavioral measure (e.g. familiarity with the Bach's musical piece)
                ad = find(isnan(a)); %this is for preventing t-tests to be calculated considering NaNs.. if there are any subject in a or in b corresponding to NaN both a and b are deprivated of that subject
                ab = find(isnan(b));
                abd = [ab ad];
                a(abd) = [];
                b(abd) = [];
                [rho,pval] = corr(a,b);
                P(ii,jj,cc,ee) = pval;
                R(ii,jj,cc,ee) = rho;
            end
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Papers/Bach_Elderly_CommunicationsBiology_Revision_I/MCS_Correlations';
mkdir(outdir)
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:size(dum2,3) %over correlations
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc,1)<p_thresh) = 1; %only testing familiarity (as asked by the reviewer; the other correlations are for my own curiosity)
        rvals = R(ii,:,cc,1);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), rvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        save([outdir '/' ROIN{ii} '_Cond' num2str(cc) '.mat'],'PDn'); %printing excel file
        writetable(PDn,[outdir '/' ROIN{ii} '_Cond' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
    end
end

%% Plotting (correlations)

clear conds
conds{1} = ' Fam '; conds{2} = ' Lik '; conds{3} = ' Emo '; conds{4} = ' Aro ';
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
for ii = 1:size(R,1) %over ROIs
    for cc = 1:size(dum2,3) %over conditions
        figure
        for ee = 1%:size(R,4) %over correlations
            plot(squeeze(time(1:size(R,2))),squeeze(R(ii,:,cc,ee)),'LineWidth',1.5,'DisplayName',[conds{ee}])
            hold on
        end
        grid minor
        title([ROIN{ii} ' Cond ' num2str(cc)])
        xlim([-0.1 3.4])
        %     ylim(limmy)
        set(gcf,'color','w')
        legend('show')
    end
end

%%

%%

%% REDOING THE ANOVA 2-WAY FOR WM AND AGE SO THAT I CAN REPORT INDEPENDENLTY THE TWO MAIN EFFECT AND THE INTERACTION (WITH THE LINE BELOW THE PLOT AS YOU DID IN NATURE COMMUNICATIONS) 

% *** TWO-WAY ANOVA TESTING YOUNG AND ELDERLY AND HIGH AND LOW WM ***

%%

gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68]; %high WM elderly
gsubj{2} = [2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM young
gsubj{3} = [3,4,6,8,18,20,23,26,27,29,31,34,36,38,43,48,49,58,61,62,66]; %low WM elderly
gsubj{4} = [5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM young
%building cell arrays WM and AGE with label corresponding to the different age and WM groups
WM = cell(72,1); %hard coding because I know there are 72 participants (1 excluded for age reasons, 1 for technical reasons and 4 because they did not have the WM)
AGE = cell(size(WM,1),1);
cnt = 0;
for pp = 1:length(gsubj{1}) %over young
    cnt = cnt + 1;
    WM{cnt} = 'hwm';
    AGE{cnt} = 'e';
end
for pp = 1:length(gsubj{2}) %over elderly (60-68)
    cnt = cnt + 1;
    WM{cnt} = 'hwm';
    AGE{cnt} = 'y';
end
for pp = 1:length(gsubj{3}) %over elderly (>68)
    cnt = cnt + 1;
    WM{cnt} = 'lwm';
    AGE{cnt} = 'e';
end
for pp = 1:length(gsubj{4}) %over elderly (>68)
    cnt = cnt + 1;
    WM{cnt} = 'lwm';
    AGE{cnt} = 'y';
end
%this is the file that you have to load and use to run the statistics
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat');
clear ROIN
ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
%t-tests
P = zeros(size(dum2,1),size(dum2,2),size(dum2,3),3); %ROIs x time-points x conditions x ANOVA's effects
F = zeros(size(dum2,1),size(dum2,2),size(dum2,3),3); %ROIs x time-points x conditions x ANOVA's effects
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %ove time-points
        for cc = 1:size(dum2,3) %over experimental conditions
            datadum = [];
            %preparing data for ANOVA
            datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc,gsubj{1}))); %extracting data for group 1
            datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc,gsubj{2}))); %group 2
            datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc,gsubj{3}))); %group 3
            datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc,gsubj{4}))); %group 4
            ab = find(isnan(datadum)); %checking for NaNs
            if ~isempty(ab) %if there are NaNs
                datadum(ab) = []; %removing them
                AGE2 = AGE;
                AGE2(ab) = []; %and removing corresponding labels
                WM2 = WM;
                WM2(ab) = []; %and removing corresponding labels
            else
                AGE2 = AGE;
                WM2 = WM;
            end
            %ANOVA
            %two-way ANOVA (anovan Matlab function)
            [p,f,stats] = anovan(datadum,{AGE2,WM2},'model','interaction','varnames',{'age','wm'},'display','off');
            P(ii,jj,cc,1) = p(1); %age (main effect) - pvalue
            F(ii,jj,cc,1) = f{2,6}; %age (main effect) - Fvalue
            P(ii,jj,cc,2) = p(2); %wm (main effect) - pvalue
            F(ii,jj,cc,2) = f{3,6}; %wm (main effect) - Fvalue
            P(ii,jj,cc,3) = p(3); %age x wm (interaction effect) - pvalue
            F(ii,jj,cc,3) = f{4,6}; %age x wm (interaction effect) - Fvalue
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
EFF{1} = 'Age'; EFF{2} = 'WM'; EFF{3} = 'Interaction';
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/aux/MINDLAB2023_MEG-AuditMemDement/Up_Down/Stats_WM_Updated';
mkdir(outdir)
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:size(F,3) %over experimental conditions
        for ee = 1:size(F,4) %over effects of the anova
            Pbin = zeros(1,size(dum2,2));
            Pbin(P(ii,:,cc,ee)<p_thresh) = 1;
            F_sel = F(ii,:,cc,ee);
            [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), F_sel ) %Monte Carlo simulation function to correct for multiple comparisons
            PDn = cell2table(sign_clust); %table
            save([outdir '/' ROIN{ii} '_Cond' num2str(cc) '_' EFF{ee} '.mat'],'PDn'); %printing excel file
            writetable(PDn,[outdir '/' ROIN{ii} '_Cond' num2str(cc) '_' EFF{ee} '.xlsx'],'Sheet',1); %printing excel file
        end
    end
end

%%

%%

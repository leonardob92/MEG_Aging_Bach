
%% Neural differences between young and older adults in automatic and conscious auditory predictive coding (using musical sequences)

%% Maxfilter

%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script
maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2021_MEG-TempSeqAges';
maxDir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter'; %output path
movement_comp = 1; %1 = yes; 0 = no

path = '/raw/sorted/MINDLAB2021_MEG-TempSeqAges'; %path with all the subjects folders
jj = dir([path '/0*']); %list all the folders starting with '0' in order to avoid hidden files
for ii = 75%:length(jj) %over subjects
    if ii ~= 49
        cart = [jj(ii).folder '/' jj(ii).name]; %create a path that combines the folder name and the file name
        pnana = dir([cart '/2*']); %search for folders starting with '2'
        for pp = 1:length(pnana) %loop to explore ad analyze all the folders inside the path above
            cart2 = [pnana(1).folder '/' pnana(pp).name];
            pr = dir([cart2 '/ME*']); %looks for meg folder
            if ~isempty(pr) %if pr is not empty, proceed with subfolders inside the meg path
                if ii ~= 48 %fixing issue of SUBJ0049 archived as additional 5 blocks of SBUJ0048
                    pnunu = dir([pr(1).folder '/' pr(1).name '/00*']);
                else
                    pnunu = dir([pr(1).folder '/' pr(1).name '/0*']); %only one '0' here since the folder has 10 files and so '009-010'
                    movement_comp = 0; %we cannot do movement compensation for SUBJ0049 since we have movement data of SUBJ0048
                end
                for dd = 1:length(pnunu)
                    fpath = dir([pnunu(1).folder '/' pnunu(dd).name '/files/*.fif']); % looks for .fif file
                    rawName = ([fpath.folder '/' fpath.name]); %assigns the final path of the .fif file to the rawName path used in the maxfilter command
                    if ii ~= 48 %SUBJ0049 was not archived properly (not loaded its preparation..), so we needed a manual operation to deal with it, as reported in the next section
                        maxfName = ['SUBJ' jj(ii).name '_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                    else
                        if dd < 6 %first 5 block are of SUBJ0048
                            maxfName = ['SUBJ' jj(ii).name '_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                        else %blocks from 6 to 10 are of SUBJ0049
                            maxfName = ['SUBJ0049_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                        end
                    end
                    if movement_comp == 1
                        %movement compensation
                        cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    else %no movement compensation (to be used if HPI coils did not work properly)
                        cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    end
                    system(cmd);
                end
            end
        end
    end
end

%% LBPD_startup_D

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path where is the function for submit the jobs to the server

%% Converting the .fif files into SPM objects

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%% conversion to SPM objects

fif_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/*.fif'); %path to SPM objects

for ii = 372:376%length(fif_list) %over the .fif files
    S = []; %structure 'S'                   
    S.dataset = [fif_list(ii).folder '/' fif_list(ii).name];
    D = spm_eeg_convert(S);
%     D = job2cluster(@cluster_spmobject, S); %actual function for conversion
end

%% Removing bad segments using OSLVIEW

%checks data for potential bad segments (periods)
%marking is done by right-clicking in the proximity of the event and click on 'mark event'
%a first click (green dashed label) marks the beginning of a bad period
%a second click indicates the end of a bad period (red)
%this will mean that we are not using about half of the data, but with such bad artefacts this is the best we can do
%we can still obtain good results with what remains
%NB: Push the disk button to save to disk (no prefix will be added, same name is kept)

%OBS! remember to check for bad segments of the signal both at 'megplanar' and 'megmag' channels (you can change the channels in the OSLVIEW interface)
%OBS! remember to mark the trial within the bad segments as 'badtrials' and use the label for removing them from the Averaging (after Epoching) 
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*.mat'); %path to SPM objects

for ii = 31%372:length(spm_list) %over experimental blocks %OBS!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = oslview(D);
    D.save(); %save the selected bad segments and/or channels in OSLVIEW
    disp(ii)
end

%% UPDATE FIDUCIALS FOR SUBJECT0049

%hard-coding solution for a subject saved with the wrong fiducials
%loading empty room with SUBJ0049 fiducials
Dfid = spm_eeg_load('/scratch7/MINDLAB2021_MEG-TempSeqAges/spmeeg_emptyroom.mat');
fid = Dfid.fiducials; %extracting fiducials
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg_SUBJ0049*.mat'); %path to SPM objects of SUBJ0049
%loading files preprocessed for SUBJ0049

for ii = 1:length(spm_list)
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    %pasting fiducials (here you need also to adapt the function "subsasgn.m" (you find comments there to help you..))
    D.fiducials = fid;
    D.save();
    disp(ii)
end

%% AFRICA denoising (part I)

%setting up the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster

%% ICA calculation (part I)

spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*.mat');

for ii = 372:length(spm_list) %OBS!
    S = [];
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    S.D = D;
    jobid = job2cluster(@cluster_africa,S);
%   D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%   D.save();
end

%% AFRICA denoising (part II)

% v = [11 12 19 32];
%visual inspection and removal of artifacted components
%look for EOG and ECG channels (usually the most correlated ones, but check a few more just in case)
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*.mat');

for ii = 372:length(spm_list) %OBS!%38:41
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = osl_africa(D,'do_ident','manual','do_remove',false,'artefact_channels',{'EOG','ECG'});
    %hacking the function to manage to get around the OUT OF MEMORY problem..
    S = [];
    S.D = D;
    jobid = job2cluster(@cluster_rembadcomp,S);
%   D.save();
    disp(ii)
end

%% Epoching: one epoch per old/new excerpt (baseline = (-)100ms)

prefix_tobeadded = 'e'; %adds this prefix to epoched files
spm_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/spmeeg*.mat');

for ii = 372:length(spm_list) %over .mat files
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files
    D = D.montage('switch',0);
    dummy = D.fname; %OBS! D.fname does not work, so we need to use a 'dummy' variable instead
    %if strcmp(dummy(22:26), 'speed') %checks whether characters 22 to 26 are equal to 'speed'; the loop continues if this is true (1) and it stops if this is false (0)
    events = D.events; %look for triggers
    if strcmp('mmn',spm_list(ii).name(17:19)) %if the block is the MMN block..
        pretrig = -100; % epoch start in ms (prestimulus/trigger)
        posttrig = 1200; % epoch end in ms (posttrigger/stimulus)
        %codes for a later (end of this section) detection of (possible) bad trials
        clear trl_sec
        ck = 0;
        for kkk = 1:length(events) %over events
            if strcmp(events(kkk).type,'STI101_up') %with only STI 014_up
                ck = ck + 1;
                trl_sec(ck,1) = events(kkk).time;
                trl_sec(ck,2) = events(kkk).time + posttrig./1000; %divided by 1000 to get time in seconds
            end
        end
        %settings for the epoching
        S2 = [];
        S2.timewin = [pretrig posttrig]; %creating the timewindow of interest
        S2.D = D;
        % event definitions
        S2.trialdef(1).conditionlabel = 'Glob_Stand';     % name of the category
        S2.trialdef(1).eventtype = 'STI101_up';         % the event type based on D for you trial/experiment
        S2.trialdef(1).eventvalue = [11,21,32,42];                  % event value based on D for your trial/experiment
        S2.trialdef(2).conditionlabel = 'Glob_Dev';
        S2.trialdef(2).eventtype = 'STI101_up';
        S2.trialdef(2).eventvalue = [12,22,31,41];
        S2.trialdef(3).conditionlabel = 'Loc_Stand';
        S2.trialdef(3).eventtype = 'STI101_up';
        S2.trialdef(3).eventvalue = [11,21,31,41];
        S2.trialdef(4).conditionlabel = 'Loc_Dev';
        S2.trialdef(4).eventtype = 'STI101_up';
        S2.trialdef(4).eventvalue = [12,22,32,42];
        S2.trialdef(5).conditionlabel = '100';
        S2.trialdef(5).eventtype = 'STI101_up';
        S2.trialdef(5).eventvalue = 100;
        %other settings
        S2.prefix = prefix_tobeadded;
        S2.reviewtrials = 0;
        S2.save = 0;
        S2.epochinfo.padding = 0;
        S2.event = D.events;
        S2.fsample = D.fsample;
        S2.timeonset = D.timeonset;
        %function to define trials specifications
        [epochinfo.trl, epochinfo.conditionlabels, MT] = spm_eeg_definetrial(S2);
        %function for the actual epoching
        S3 = [];
        S3.prefix = S2.prefix;
        S3.epochinfo = epochinfo;
        S3.D = D;
        D = osl_epoch(S3);
        %     D = spm_eeg_load([spm_list(ii).folder '/e' spm_list(ii).name]); %loading epoched data
        D = D.montage('switch',1);
        %from here you go at the end of the if state
    elseif ~strcmp('res',spm_list(ii).name(17:19)) %otherwise, if it is not the resting state..
        %takes the correct triggers sent during the recording
        clear trigcor
        count_evval = 0; %???
        for ieve = 1:length(events) %over triggers
            if strcmp(events(ieve).type,'STI101_up') %only triggers at the beginning of each stimuli
                if strcmp('aud',spm_list(ii).name(17:19)) || strcmp('vis',spm_list(ii).name(17:19))
                    %if events(ieve).value ~= 103 && events(ieve).value ~= 104 && events(ieve).value ~= 128 && events(ieve).value ~= 8 && events(ieve).value ~= 132 && events(ieve).value ~= 48 && events(ieve).value ~= 32 && events(ieve).value ~= 64 %discard 104 and 128 for random triggers
                    if events(ieve).value == 11 || events(ieve).value == 21 %10 and 50 are old and new in recogminor (block 3), while 11 and 21 are old and new in blocks 4 and 5 (aud and vis)
                        count_evval = count_evval + 1;
                        trigcor(count_evval,1) = events(ieve).time; %+ 0.010; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes
                        %variable with all the triggers we need
                    end
                elseif strcmp('rec',spm_list(ii).name(17:19)) || strcmp('rac',spm_list(ii).name(17:19))
                    if events(ieve).value == 10 || events(ieve).value == 50 %10 and 50 are old and new in recogminor (block 3), while 11 and 21 are old and new in blocks 4 and 5 (aud and vis)
                        count_evval = count_evval + 1;
                        trigcor(count_evval,1) = events(ieve).time; %+ 0.010; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes
                        %variable with all the triggers we need
                    end
                end
            end
        end
        trl_sam = zeros(length(trigcor),3); %prepare the samples matrix with 0's in all its cells
        trl_sec = zeros(length(trigcor),3); %prepare the seconds matrix with 0's in all its cells
        %deftrig = zeros(length(trigcor),1); %this is not useful
        for k = 1:length(trigcor) %over selected triggers
            %deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
            trl_sec(k,1) = trigcor(k,1) - 0.1000; %beginning time-window epoch in s (please note that we computed this operation two times, obtaining two slightly different pre-stimulus times.
            %this was done because for some computations was convenient to have a slightly longer pre-stimulus time
            %remove 1000ms of baseline
            trl_sec(k,2) = trigcor(k,1) + 4.4; %end time-window epoch in seconds
            trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in seconds
            trl_sam(k,1) = round(trl_sec(k,1) * 250) + 1; %beginning time-window epoch in samples %250Hz per second
            trl_sam(k,2) = round(trl_sec(k,2) * 250) + 1; %end time-window epoch in samples
            trl_sam(k,3) = -25; %sample before the onset of the stimulus (corresponds to 0.100ms)
        end
        dif = trl_sam(:,2) - trl_sam(:, 1); %difference between the end and the beginning of each sample (just to make sure that everything is fine)
        if ~all(dif == dif(1)) %checking if every element of the vector are the same (i.e. the length of the trials is the same; we may have 1 sample of difference sometimes because of different rounding operations..)
            trl_sam(:,2) = trl_sam(:,1) + dif(1);
        end
        %creates the epochinfo structure that is required for the source reconstruction later
        epochinfo.trl = trl_sam;
        epochinfo.time_continuous = D.time;
        %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
        D = D.montage('switch',0);
        %build structure for spm_eeg_epochs
        S = [];
        S.D = D;
        S.trl = trl_sam;
        S.prefix = prefix_tobeadded;
        D = spm_eeg_epochs(S);
        %store the epochinfo structure inside the D object
        D.epochinfo = epochinfo;
        D.save();
    end
    %take bad segments registered in OSLVIEW and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later
    count = 0;
    Bad_trials = zeros(length(trl_sec),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + events(kkk).duration) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        count = count + 1;
                    end
                end
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    disp(spm_list(ii).name);
    if count == 0
        disp('there are no bad trials marked in oslview');
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        %         D = conditions(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        epochinfo = D.epochinfo;
        xcv = find(Bad_trials == 1);
        %this should be done only later.. in any case.. not a problem..
        for jhk = 1:length(xcv)
            D = D.conditions(xcv(jhk),'Bad');
            epochinfo.conditionlabels(xcv(jhk)) = {'Bad'};
            disp([num2str(ii) ' - ' num2str(jhk) ' / ' num2str(length(xcv))])
        end
        D.epochinfo = epochinfo;
        D.save(); %saving on disk
        disp('bad trials are ')
        length(D.badtrials)
    end
    D.save();
    disp(ii)
end

%% MERGING MMN AND MMN2 OF SUBJECT0073

%fixing an issue that led us to save two files for the same paradigm of SUBJ0073

epoch_list2 = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/espmeeg_SUBJ0073_mmn*mat'); %epoched files

S = []; %empty structure
S.recode = 'same';
S.prefix = 'emerged_mmn_'; %merged
for ss = 1:2 %over blocks for subject ii
    S.D{ss} = spm_eeg_load([epoch_list2(ss).folder '/' epoch_list2(ss).name]);
end
Dout = spm_eeg_merge(S); %actual function
%THEN, I MANUALLY DELETED SUBJ0073_mmn and SUBJ0073_mmn2, the input of the merging

%% Defining the conditions - only for old/new paradigm (in the case of MMN the conditions were defined during epoching)

xlsx_dir_behav = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/MEG_behavioural'; %dir to MEG behavioral results (.xlsx files)
epoch_list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*.mat'); %dir to epoched files
for ii = 297:length(epoch_list) %over epoched data
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
            for k = 1:length(D.trialonset)
                if raw_recog{(k + 1),3} == 0 %if there was no response
                    D = D.conditions(k,'No_response');
                elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
                    D = D.conditions(k,'Old_Correct'); %assign old correct
                elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 2 %old incorrect
                    D = D.conditions(k,'Old_Incorrect'); %otherwise assign new correct
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
                    D = D.conditions(k,'New_T1_Correct');
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 1 %new t1 incorrect
                    D = D.conditions(k,'New_T1_Incorrect');
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
                    D = D.conditions(k,'New_T3_Correct');
                elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 1 %new t3 incorrect
                    D = D.conditions(k,'New_T3_Incorrect');
                end
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

%% *** MEG BEHAVIORAL STATISTICS ***

% see script developed in R

%%

%% *** MASSIVE UNIVARIATE ANALYSIS AND MCS ON MEG SENSORS ***

%% Averaging and Combining planar gradiometers

%settings for cluster (parallel computing)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster'); %set automatically the long run queue
clusterconfig('long_running', 1); %set automatically the long run queue
clusterconfig('slot', 1); %set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)

%% averaging

output_prefix_to_be_set = 'm';

v = [297:312]; %only a selection of files
epoch_list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*mat'); %dir to epoched files (encoding)
for ii = 1:length(v)%1:length(epoch_list) %over epoched files
    input = [];
    input.D = [epoch_list(v(ii)).folder '/' epoch_list(v(ii)).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
    % look the script for more details about the function work
end

%% combining planar gradiometers

average_list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/m*mat'); %dir to epoched files (encoding)
v = [297:312]; %only a selection of files
for ii = 1:length(v)%1:length(average_list) %over files
    input = [];
    input.D = [average_list(v(ii)).folder '/' average_list(v(ii)).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1);
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
end

%% LBPD_startup_D

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path where is the function for submit the jobs to the server

%% Extracting MEG sensor data

block = 4; % 2 = MMN; 3 = recogminor; 4 = aud (numerical sequence); 5 = vis (numerical sequence); 6 = imagination task (auditory and visual together)
% channels_plot = 13; %13;95;;9; %95, 13, 15, 11, 141, 101, 9 % empty for plotting single channels; otherwise number(s) of channels to be averaged and plotted (e.g. [13] or [13 18])
channels_plot = []; % empty for plotting single channels; otherwise number(s) of channels to be averaged and plotted (e.g. [13] or [13 18])
subj = []; %[gsubj{2}(34)]; %empty for all subjects; numbers of subjects with regards to the list indexed here if you want specific subjects

% subj = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78]; %young
% subj = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,35,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76]; %elderly

%1321 1411
waveform_singlechannels_label = 1; %1 to plot individual channels
save_data = 0; %1 to save the data on disk
load_data = 1; %set 1 if you want to load the data instead of extracting it from SPM objects
% v = [1]; %subjects
%bad 8,9 and a bit 6 (recogminor)

S = [];
%computing data
if block == 2
    S.conditions = {'Glob_Stand','Glob_Dev','Loc_Stand','Loc_Dev'};
%     S.conditions = {'Glob_Stand','Glob_Dev','Loc_Stand','Loc_Dev','100'}; 
elseif block == 3
    S.conditions = {'Old_Correct','New_T1_Correct','New_T3_Correct'};
%     S.conditions = {'Old_Correct','New_T1_Correct'};
%     S.conditions = {'Old_Incorrect','New_T1_Incorrect','New_T3_Incorrect'};
end
if block == 2
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/P*mmn*_tsssdsm.mat'); %dir to epoched files (encoding)
    clear listdum
    listdum(1:72,1) = list(2:73); listdum(73,1) = list(1); listdum(74:78,1) = list(74:78);
    list = listdum;
    S.y_lim_ampl_wave = [-180 180]; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
    S.x_lim_temp_wave = []; %limits for time (in secs) (E.g. [-0.1 3.4])
elseif block == 3
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/P*cogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
    S.y_lim_ampl_wave = [-150 150]; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
    S.x_lim_temp_wave = [-0.1 3.4]; %limits for time (in secs) (E.g. [-0.1 3.4])
end
if isempty(subj)
    v = 1:length(list); %subjects
else
    v = subj;
%     S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
end
if ~exist('chanlabels','var')
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter/MEG_sensors/recogminor_all_conditions.mat', 'chanlabels')
end
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors'; %path to write output in
S.outdir = outdir;
S.data = [];
if load_data == 1 %if you already computed and saved on disk the t-tests you can load them here
    load([outdir '/Block_' num2str(block) '.mat'],'time_sel','data_mat','chanlabels');
    S.data = data_mat(:,:,v,:);
    S.chanlabels = chanlabels;
    S.time_real = time_sel;
else %otherwise you can extract the data from SPM MEEG objects (one for each subject)
%     S.spm_list = cell(1,length(list));
% v = 7;
    S.spm_list = cell(1,length(v));
    for ii = 1:length(v)
        S.spm_list(ii) = {[list(v(ii)).folder '/' list(v(ii)).name]};
    end
end

S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = save_data; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = ['Block_' num2str(block)];

%individual waveform plotting
if isempty(channels_plot)
    S.waveform_singlechannels_label = waveform_singlechannels_label; %1 to plot single channel waveforms
else
    S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
end
S.wave_plot_conditions_together = 0; %1 for plotting the average of all
S.mag_lab = 1; %1 for magnetometers; 2 for gradiometers

%averaged waveform plotting
if isempty(channels_plot)
    S.waveform_average_label = 0; %average of some channels
    S.left_mag = 95; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
else
    S.waveform_average_label = 1; %average of some channels
    S.left_mag = channels_plot; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
end
% S.left_mag = [2:2:204];
S.legc = 1; %set 1 for legend
% S.left_mag = 99;
S.signtp = {[]};
% S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'c';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)

%topoplotting
S.topoplot_label = 1;
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 0;
S.topocondsing = [2]; %condition for topoplot
% S.xlim = [0.75 0.85]; %time topolot
% S.xlim = [1.1 1.2]; %time topolot
S.xlim = [3.5 4.3]; 
S.zlimmag = []; %magnetometers amplitude topoplot limits
S.zlimgrad = []; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);

%% MEG sensors - plotting elderly versus young participants (this allows to plot data from other experimental blocks which are not included in the current manuscript)

%%% OBS!! THIS SECTION IS NOT USED FOR CREATING THE IMAGES OF THE MANUSCRIPT %%%

block = 3; %1 = mmn subtracted; 2 = mmn no subtracted; 3 = recogminor; 4 = auditory sequence of numbers; 5 = visual sequence of numbers
condition = [3]; %condition number
young_elderly = 4; %1 = comparing young and elderly; 2 = comparing elderly (60-68 years old) and elderly (>68 years old); 3 = young, high vs low WM; 4 = eldely, high vs low WM
chans = 0; %0 = all channels; or channel index (e.g. 13,31)


if block == 1
    load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/Block_2_MMNsubtracted.mat']); %loading the channels forming the main significant cluster across all deviants
else
    load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/Block_' num2str(block) '.mat']); %loading the channels forming the main significant cluster across all deviants
end

%structure for the function
S = [];
S.data = data_mat;
if young_elderly == 1
    S.groups = {'young','elderly'}; %Group labels
    gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
    gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,35,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76];
elseif young_elderly == 2
    S.groups = {'60-68','>68','young'}; %Group labels
    gsubj{1} = [3,8,10,11,12,20,26,31,34,36,42,45,47,49,51,52,56,58,59,61,62,68]; %60-68
    gsubj{2} = [4,6,7,16,18,23,27,28,29,30,32,35,38,43,48,66,76]; %>68
    gsubj{3} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
elseif young_elderly == 3 %young, high WM vs low WM
    S.groups = {'HWMY','LWMY'}; %Group labels
    gsubj{1} = [2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM
    gsubj{2} = [5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM
elseif young_elderly == 4 %elderly, high WM vs low WM
    S.groups = {'HWME','LWME'}; %Group labels
    gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68]; %high WM
    gsubj{2} = [3,4,6,8,18,20,23,26,27,29,31,34,35,36,38,43,48,49,58,61,62,66]; %low WM
end
S.gsubj = gsubj;
S.signtp = [];
S.time_real = time_sel;
S.legendl = 1;
S.colorline = {'r', 'b', 'k'}; %Select colorline for each group
S.sensors = -1; % -1 for magnetometer, 0 for gradiometers (gradiometers have already been combined)
if block == 1
    S.conditions = {'Glob_mmn','Loc_mmn'};
    S.x_lim = [-0.05 1.2]; % Set x limits
    S.y_lim = [-100 100]; %Set y limits
elseif block == 2
    S.conditions = {'Glob_Stand','Glob_Dev','Loc_Stand','Loc_Dev'};
    %     S.conditions = {'Glob_Stand','Glob_Dev','Loc_Stand','Loc_Dev','100'};
    S.x_lim = [-0.05 1.2]; % Set x limits
    S.y_lim = [-180 180]; %Set y limits
elseif block == 3
    S.conditions = {'Old_Correct','New_T1_Correct','New_T3_Correct'};
    %     S.conditions = {'Old_Correct','New_T1_Correct'};
    S.y_lim = [-200 200]; %Set y limits
    S.x_lim = [-0.1 3.4]; % Set x limits
elseif block > 3
    S.conditions = {'Old_Correct','New_T1_Correct','New_T3_Correct'};
    %     S.conditions = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
    S.y_lim = [-200 200]; %Set y limits
    S.x_lim = [-0.1 3.4]; % Set x limits
end
S.condition_n = condition;
S.chanlabels = chanlabels;
S.STE = 2;
S.chans_index = chans; %plot waveform at single channel or the average across multiple channels (specify channel index - you can find the channel labels in 'chanlabels'); set to 0 to plot all channels

plot_sensors_wavebis2(S) %actual function

%% T-TESTS BETWEEN GROUPS (young vs elderly)

%%% SUBJ0035 HAS BEEN REMOVED SINCE IT DID NOT HAVE DATA FOR NEWT1 AND NEWT3

condition = 3; %1 = Old; 2 = NewT1; 3 = NewT3 (Block 3)

load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/Block_3.mat')
gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78]; %young
gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76]; %elderly
data_mat(1:2:204,:,:,:) = abs(data_mat(1:2:204,:,:,:)); %to avoid the sign ambiguity of the magnetometers and allow easier comparisons between the groups of participants
T = zeros(size(data_mat,1),size(data_mat,2));
Po = zeros(size(data_mat,1),size(data_mat,2));
for jsk = 1:size(data_mat,1) %over channels
    for jsl = 1:size(data_mat,2) %over time-points
        a = squeeze(data_mat(jsk,jsl,gsubj{1},condition)); %young
        b = squeeze(data_mat(jsk,jsl,gsubj{2},condition)); %elderly
        ad = find(isnan(a)); %this is for preventing t-tests to be calculated considering NaNs.. if there are any subject in a or in b corresponding to NaN both a and b are deprivated of that subject
        ab = find(isnan(b));
        abd = [ab ad];
        a(abd) = [];
        b(abd) = [];
        [~,pttest,~,tstats] = ttest2(a,b);
        T(jsk,jsl) = tstats.tstat;
        Po(jsk,jsl) = pttest;
    end
    disp(jsk)
end

% Monte-Carlo simulations (MCS)
%input information
p_thresh = 0.05; %for binarising p-values matrices..

%actual computation
%time-points to be selected
min_time_point = 26;
max_time_point = 526;
clear DATAP2 TSTAT2
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS'; %path where t-test results are stored in
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/Block_3.mat','time_sel','chanlabels')
%here gradiometers and magnetometers are always extracted but in the
%following steps only the requested channels (either magnetometers or
%gradiometers) are used
DATAP2(:,:,1) = Po(1:2:204,min_time_point:max_time_point); %magnetometers p-values
DATAP2(:,:,2) = Po(2:2:204,min_time_point:max_time_point); %gradiometers (combined) p-values
TSTAT2(:,:,1) = T(1:2:204,min_time_point:max_time_point); %mag t-vals
TSTAT2(:,:,2) = T(2:2:204,min_time_point:max_time_point); %grad t-vals
chanlab = chanlabels(1:2:204)';
label = zeros(length(chanlab),1); %channels label
for ii = 1:length(chanlab)
    label(ii,1) = str2double(chanlab{ii,1}(4:end));
end

%individuating positive vs negative t-values
P = DATAP2; %trick for looking only into the positive t-values
P(P < p_thresh) = 1; %binarizing p-values according to threshold
P(P < 1) = 0;
%old (pos) > new (neg)
P(TSTAT2 < 0) = 0; %deleting p-values for negative contrasts
TSTAT_mag_pos = P(:,:,1);
TSTAT_grad = P(:,:,2);
%negative vs positive p-values
P = DATAP2;
P(P < p_thresh) = 1; %binarizing p-values according to threshold
P(P < 1) = 0;
P(TSTAT2 > 0) = 0; %deleting p-values for positive contrasts
TSTAT_mag_neg = P(:,:,1);
%load a 2D approximation in a matrix of the MEG channels location (IT IS CONTAINED IN THE FUNCTIONS FOLDER THAT WE PROVIDED)
[~,~,raw_channels] = xlsread('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MatrixMEGChannelLayout_With0_2.xlsx');
% reshaping data for computational purposes
S = [];
S.label = label;
S.TSTAT_mag_pos = TSTAT_mag_pos;
S.TSTAT_mag_neg = TSTAT_mag_neg;
S.TSTAT_grad = TSTAT_grad;
S.TVal = TSTAT2;
S.raw_channels = raw_channels;

[MAG_data_pos, MAG_data_neg, GRAD_data, MAG_GRAD_tval] = MEG_sensors_MCS_reshapingdata_LBPD_D(S);

%actual Monte Carlo simulations
S = [];
%actual gradiometers data
S.time = time_sel(min_time_point:max_time_point);
S.data(:,:,:,1) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3)); %GRAD_data;
%zeros.. if you do not want to run the function for magnetometers
S.data(:,:,:,2) = MAG_data_pos;
S.data(:,:,:,3) = MAG_data_neg;
S.sensortype = [];
S.MAG_GRAD_tval = MAG_GRAD_tval;
S.MEGlayout = cell2mat(raw_channels);
S.permut = 1000;
S.clustmax = 1;
S.clust_sort = 1; %1 to sort according to size; 2 according to max-tval (in absolute terms)
S.permthresh = 0.001;

[MAG_clust_pos, MAG_clust_neg, GRAD_clust] = MEG_sensors_MonteCarlosim_LBPD_D(S);

save([outdir '/MAG_Cond' num2str(condition) '_YoungvsEld_tvalp_05.mat'],'MAG_clust_pos','MAG_clust_neg');

%% exctracting information for MEG sensor clusters

list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS/*.mat');
for ii = 1:3 %over experimental conditions
    load([list(ii).folder '/' list(ii).name]) %loading clusters outputted from MCS on MEG sensors
    [ PDn, PDnchan ] = Extract_MEGSensor_Information_LBPD_D( MAG_clust_pos );
    writetable(PDn,['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS/Images/' list(ii).name(1:end-4) 'pos_overall.xlsx'],'Sheet',1); %printing excel file
    writetable(PDnchan,['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS/Images/' list(ii).name(1:end-4) 'pos_channels.xlsx'],'Sheet',1); %printing excel file
    [ PDn, PDnchan ] = Extract_MEGSensor_Information_LBPD_D( MAG_clust_neg );
    writetable(PDn,['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS/Images/' list(ii).name(1:end-4) 'neg_overall.xlsx'],'Sheet',1); %printing excel file
    writetable(PDnchan,['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS/Images/' list(ii).name(1:end-4) 'neg_channels.xlsx'],'Sheet',1); %printing excel file
end

%%

%% *** MEG SOURCE RECONSTRUCTION - BEAMFORMING ***

%%

%% CREATING 8mm PARCELLATION FOR EASIER INSPECTION IN FSLEYES
%OBS!! This section is done only for better handling of some visualization purposes, but it does not affect any of the beamforming algorithm;
% it is just important not to mix up the MNI coordinates, thus I would recommend to use the following lines

%1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
imag_8mm = load_nii('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_8mm_brain.nii.gz');
Minfo = size(imag_8mm.img); %get info about the size of the original image
M8 = zeros(Minfo(1), Minfo(2), Minfo(3)); %Initialize an empty matrix with the same dimensions as the original .nii image
cc = 0; %set a counter
M1 = imag_8mm.img;
for ii = 1:Minfo(1) %loop across each voxel of every dimension
    for jj = 1:Minfo(2)
        for zz = 1:Minfo(3)
            if M1(ii,jj,zz) ~= 0 %if we have an actual brain voxel
                cc = cc+1;
                M8(ii,jj,zz) = cc;
            end
        end
    end
end
%2) PUT YOUR MATRIX IN THE FIELD ".img"
imag_8mm.img = M8; %assign values to new matrix 
%3) SAVE NIFTI IMAGE USING save_nii
save_nii(imag_8mm,'/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.gz');
%4) USE FSLEYES TO LOOK AT THE FIGURE
%Create parcellation on the 8mm template
for ii = 1:3559 %for each 8mm voxel
    cmd = ['fslmaths /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_80mm_3559ROIs/' num2str(ii) '.nii.gz'];
    system(cmd)
    disp(ii)
end
%5) GET MNI COORDINATES OF THE NEW FIGURE AND SAVE THEM ON DISK
MNI8 = zeros(3559,3);
for mm = 1:3559 %over brain voxel
    path_8mm = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/' num2str(mm) '.nii.gz']; %path for each of the 3559 parcels
    [mni_coord,pkfo] = osl_mnimask2mnicoords(path_8mm);  %getting MNI coordinates
    MNI8(mm,:) = mni_coord; %storing MNI coordinates
end
%saving on disk
save('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat', 'MNI8');

%% CONVERSION T1 - DICOM TO NIFTI

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/dicm2nii'); %adds path to the dcm2nii folder in osl
MRIsubj = dir('/projects/MINDLAB2021_MEG-TempSeqAges/raw/0*');
MRIoutput = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/MRI_nifti';
MRIout_block{1} = 'Block_1'; MRIout_block{2} = 'Block_2'; MRIout_block{3} = 'Block_3'; MRIout_block{4} = 'Block_4'; MRIout_block{5} = 'Block_5'; MRIout_block{6} = 'Block_6';

for bb = 1:length(MRIout_block) %over experimental blocks
    for ii = 1:length(MRIsubj) %over subjects
        asd = [MRIoutput '/Block_' num2str(bb) '/' MRIsubj(ii).name];
        if ~exist(asd,'dir') %checking whether the directory exists
            mkdir(asd); %if not, creating it
        end
        if strcmp(MRIsubj(ii).name,'0072') %we previously obtained subject 0072's MRI so we convert it from a different folder
            dcmSource = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0011/20201208_173601/MR/005.t1_mprage_3D_sag_fatsat/files/';
            niiFolder = asd;
            dicm2nii(dcmSource, niiFolder, '.nii');
        else %actual subjects (only) from the new data collection
            if isempty(dir([asd '/*.nii'])) %if there are no nifti images.. I need to convert them
                flagg = 0;
                MRIMEGdate = dir([MRIsubj(ii).folder '/' MRIsubj(ii).name '/20*']);
                niiFolder = [MRIoutput '/' MRIout_block{bb} '/' MRIsubj(ii).name];
                for jj = 1:length(MRIMEGdate) %over dates of recording
                    %                     if ~isempty(dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR*'])) %if we get an MRI recording
                    MRI2 = dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR/*fatsat']); %looking for T1
                    if length(MRI2) == 1 %if we have one T1
                        flagg = 1; %determining that I could convert MRI T1
                        dcmSource = [MRI2(1).folder '/' MRI2(1).name '/files/'];
                        dicm2nii(dcmSource, niiFolder, '.nii');
                    end
                    if length(MRI2) ~= 1 && jj == length(MRIMEGdate)
                        warning(['subject ' MRIsubj(ii).name ' has no MRI T1 or has more than 1 MRI T1']);
                        warning('copying brain template..')
                        copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
                    end
                end
                if ~isempty(dir([niiFolder '/*.txt'])) %if something goes wrong with the conversion of the MRI file, I copy-paste a template
                    copyfile('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_2mm.nii',[niiFolder '/MNI152_T1_2mm.nii'])
                end
            end
        end
        disp(ii)
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
% clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% RHINO coregistration

%block to be run RHINO coregistrartion on
block = 2; % 2 = MMN; 3 = recogminor; 4 = aud; 5 = vis; 6 = imagery/replay (IN BLOCK 6 THERE ARE 5 SUBJECTS WHICH ARE MISSING AND THAT CAN BE RETRIEVED BY MODIFYING THE SPM OBJECTS THAT HAVE BEEN MERGED..!!)

if block == 2
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*mmn*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16; extr1 = 25:28;
elseif block == 3
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*cogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16;
elseif block == 4
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/es*aud*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16;
elseif block == 5
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/es*vis*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16;
elseif block == 6
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/me_r*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 22:25;
end

%running rhino
%OBS! check that all MEG data are in the same order and number as MRI nifti files!
a = ['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/MRI_nifti/Block_' num2str(block)]; %set path to MRI subjects' folders
for ii = 1%:length(list) %OBS! change this depending on atonal vs. major
    if block == 2 && ii == 1 %subject with merged files (and so different name..)
       extr2 = extr1; %getting name from different characters (25:28 instead of 13:16) 
    else
        extr2 = extr;
    end
    S = [];
    S.ii = ii;
    S.D = [list(ii).folder '/' list(ii).name]; %path to major files
    D = spm_eeg_load(S.D);
    if ~isfield(D,'inv') %checking if the coregistration was already run
        dummyname = D.fname;
        if 7 == exist([a '/' dummyname(extr2)],'dir') %if you have the MRI folder
            dummymri = dir([a '/' dummyname(extr2) '/*.nii']); %path to nifti files (ending with .nii)
            if ~isempty(dummymri)
                S.mri = [dummymri(1).folder '/' dummymri(1).name];
                %standard parameters
                S.useheadshape = 1;
                S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
                %         S.forward_meg = 'MEG Local Spheres';
                S.forward_meg = 'Single Shell'; %CHECK WHY IT SEEMS TO WORK ONLY WITH SINGLE SHELL!!
                S.fid.label.nasion = 'Nasion';
                S.fid.label.lpa = 'LPA';
                S.fid.label.rpa = 'RPA';
                jobid = job2cluster(@coregfunc,S); %running with parallel computing
            else
                warning(['subject ' dummyname(extr2) ' does not have the MRI'])
            end
        end
    else
        if isempty(D.inv{1}) %checking whether the coregistration was run but now it is empty..
            warning(['subject ' D.fname ' has an empty rhino..']);
        end
    end
    disp(ii)
end

%% checking (or copying) RHINO

copy_label = 1; % 1 = pasting inv RHINO from epoched data (where it was computed) to continuous data (not supported for block 6); 0 = simply showing RHINO coregistration
block = 5; % 2 = MMN; 3 = recogminor; 4 = aud; 5 = vis; 6 = imagery/replay (IN BLOCK 6 THERE ARE 5 SUBJECTS WHICH ARE MISSING AND THAT CAN BE RETRIEVED BY MODIFYING THE SPM OBJECTS THAT HAVE BEEN MERGED..!!)

if block == 2
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*mmn*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16; extr1 = 25:28;
elseif block == 3
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*cogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16;
elseif block == 4
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/es*aud*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16;
elseif block == 5
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/es*vis*_tsssdsm.mat'); %dir to epoched files (encoding)
    extr = 13:16;
elseif block == 6
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/me_r*_tsssdsm.mat'); %dir to epoched files (encoding)
end

for ii = 1:length(list)
    D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
    if isfield(D,'inv')
        if copy_label == 0 %simply displaying RHINO coregistration
            if isfield(D,'inv') %checking if the coregistration was already run
                rhino_display(D)
            end
        else %pasting inv RHINO from epoched data (where it was computed) to continuous data
            if block ~= 6
                inv_rhino = D.inv;
                if block == 2 && ii == 1
                    warning('you cannot do that for subject 1 (which would be 0073) of block 2, since it is a merged file')
                else
                    D2 = spm_eeg_load([list(ii).folder '/' list(ii).name(2:end)]);
                    D2.inv = inv_rhino;
                    D2.save();
                end
            else
                error('you cannot copy-paste rhino in the continuous file for hte imagery/replay data since these files have been obtained by merging epoched data in resting state, auditory and visual encoding of numbers')
            end
        end
    end
    disp(['Block ' num2str(block) ' - Subject ' num2str(ii)])
end

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
block = 6; % 2 = MMN; 3 = recogminor; 4 = aud; 5 = vis; 6 = imagery/replay (IN BLOCK 6 THERE ARE 5 SUBJECTS WHICH ARE MISSING AND THAT CAN BE RETRIEVED BY MODIFYING THE SPM OBJECTS THAT HAVE BEEN MERGED..!!)
invers = 1; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction
absl = 0; % 1 = absolute value of sources; 0 = not

%actual computation
%list of subjects with coregistration (RHINO - OSL/FSL) - epoched
if block == 2
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*mmn*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Glob_Stand','Loc_Stand','Glob_Dev','Loc_Dev','100'};
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_MMN.mat');
    extr = 13:16; extr1 = 25:28;
    timek = 1:326;
elseif block == 3
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/e*cogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Correct','New_T1_Correct','New_T3_Correct'};
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat');
    extr = 13:16;
elseif block == 4
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/es*aud*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat');
    extr = 13:16;
elseif block == 5
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/es*vis*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat');
    extr = 13:16;
elseif block == 6
    list = dir ('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/me_r*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Replay_Aud','Replay_Vis','Replay_Res'};
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_replay.mat');
    extr = 22:25;
end
if isempty(freqq)
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers)];
else
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers)];
end
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
if ~exist(workingdir,'dir') %creating working folder if it does not exist
    mkdir(workingdir)
end
for ii = 1:length(list) %over subjects
    if block == 2 && ii == 1 %subject with merged files (and so different name..)
       extr2 = extr1; %getting name from different characters (25:28 instead of 13:16) 
    else
        extr2 = extr;
    end
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

%% BRAIN SOURCES OF THE SIGNIFICANT CLUSTERS AT MEG SENSOR LEVEL (MAIN EFFECTS AND CONTRAST BETWEEN YOUNG VS ELDERLY)

clust_num = 10; %number of clusters you want (in decrescent size)

% 1) extracting and preparing data in proper time-windows
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS'; %path where t-test results are stored in
load([outdir '/MAG_Cond' num2str(condition) '_YoungvsEld_tvalp_05.mat']) %loading clusters outputted from MCS on MEG sensors
list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/SUBJ*.mat');
POS = zeros(3559,clust_num,length(list)); %brain voxels x clusters x subjects (young > elderly)
NEG = zeros(3559,clust_num,length(list)); %brain voxels x clusters x subjects (elderly > young)
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    load([list(ii).folder '/' list(ii).name])
    for cc = 1:clust_num %over clusters that you want
        POS(:,cc,ii) = mean(OUT.sources_ERFs(:,round(MAG_clust_pos{cc+1,6}*250)+26:round(MAG_clust_pos{cc+1,7}*250)+26,condition),2); %elaborated but proper solution to get time-points (in time-samples) for cluster cc        
        NEG(:,cc,ii) = mean(OUT.sources_ERFs(:,round(MAG_clust_neg{cc+1,6}*250)+26:round(MAG_clust_neg{cc+1,7}*250)+26,condition),2); %elaborated but proper solution to get time-points (in time-samples) for cluster cc
    end
end

% 2) %making contrasts (young vs elderly)
%POS (sources of clusters where young > elderly)
TPOS = zeros(size(POS,1),size(POS,2));
PoPOS = zeros(size(POS,1),size(POS,2));
for jsk = 1:size(POS,1) %over brain voxels
    for jsl = 1:size(POS,2) %over time-windows (clusters)
        a = squeeze(POS(jsk,jsl,gsubj{1})); %young
        b = squeeze(POS(jsk,jsl,gsubj{2})); %elderly
        ad = find(isnan(a)); %this is for preventing t-tests to be calculated considering NaNs.. if there are any subject in a or in b corresponding to NaN both a and b are deprivated of that subject
        ab = find(isnan(b));
        abd = [ab ad];
        a(abd) = [];
        b(abd) = [];
        [~,pttest,~,tstats] = ttest2(a,b);
        TPOS(jsk,jsl) = tstats.tstat;
        PoPOS(jsk,jsl) = 1 - pttest;
    end
    disp(jsk)
end
%NEG (sources of clusters where elderly > young)
TNEG = zeros(size(NEG,1),size(NEG,2));
PoNEG = zeros(size(NEG,1),size(NEG,2));
for jsk = 1:size(NEG,1) %over brain voxels
    for jsl = 1:size(NEG,2) %over time-windows (clusters)
        a = squeeze(NEG(jsk,jsl,gsubj{1})); %young
        b = squeeze(NEG(jsk,jsl,gsubj{2})); %elderly
        ad = find(isnan(a)); %this is for preventing t-tests to be calculated considering NaNs.. if there are any subject in a or in b corresponding to NaN both a and b are deprivated of that subject
        ab = find(isnan(b));
        abd = [ab ad];
        a(abd) = [];
        b(abd) = [];
        [~,pttest,~,tstats] = ttest2(a,b);
        TNEG(jsk,jsl) = tstats.tstat;
        PoNEG(jsk,jsl) = 1 - pttest;
    end
    disp(jsk)
end

% 3) plotting main effects and contrasts in brain templates
clear DAT
DAT{1} = PoPOS; DAT{2} = TPOS; DAT{3} = mean(POS(:,:,gsubj{1}),3); DAT{4} = mean(POS(:,:,gsubj{2}),3); DAT{5} = mean(NEG(:,:,gsubj{1}),3); DAT{6} = mean(NEG(:,:,gsubj{2}),3); DAT{7} = PoNEG; DAT{8} = TNEG;
names{1} = ['PvalPOS_Cond' num2str(condition) '_YoungVsEld_10Clusts']; names{2} = ['TvalPOS_Cond' num2str(condition) '_YoungVsEld_10Clusts']; names{3} = ['MEPOS_Cond' num2str(condition) '_Young_10Clusts']; names{4} = ['MEPOS_Cond' num2str(condition) '_Elderly_10Clusts']; names{5} = ['MENEG_Cond' num2str(condition) '_Young_10Clusts']; names{6} = ['MENEG_Cond' num2str(condition) '_Elderly_10Clusts']; names{7} = ['PvalNEG_Cond' num2str(condition) '_YoungVsEld_10Clusts']; names{8} = ['TvalNEG_Cond' num2str(condition) '_YoungVsEld_10Clusts'];
%FROM LBPD COORDINATES TO 3D NIFTI IMAGE
for ii = 1:length(DAT)
    S = [];
    S.singleimage = 1;
    S.data = DAT{ii}; %data (voxels x ROIs (time-points))
    S.fname = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources'; %path and name for the image to be saved
    S.names = names(ii); %names for the different images
    %actual function
    FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);
end

%% PREPARING MEG SOURCE IMAGES FOR ALL PARTICIPANTS TOGETHER (NO CONTRAST BETWEEN GROUPS - THIS IS FOR PROVIDING READERS WITH AN UNDERSTANDING OF THE GENERAL BRAIN ACTIVITY OVER ALL PARTICIPANTS)

condition = 3; %1 = Old; 2 = NewT1; 3 = NewT3 (Block 3)
clust_num = 10; %number of clusters you want (in decrescent size)

% 1) extracting and preparing data in proper time-windows
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS'; %path where t-test results are stored in
load([outdir '/MAG_Cond' num2str(condition) '_YoungvsEld_tvalp_05.mat']) %loading clusters outputted from MCS on MEG sensors
list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/SUBJ*.mat');
POS = zeros(3559,clust_num,length(list)); %brain voxels x clusters x subjects (young > elderly)
NEG = zeros(3559,clust_num,length(list)); %brain voxels x clusters x subjects (elderly > young)
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    load([list(ii).folder '/' list(ii).name])
    for cc = 1:clust_num %over clusters that you want
        POS(:,cc,ii) = mean(OUT.sources_ERFs(:,round(MAG_clust_pos{cc+1,6}*250)+26:round(MAG_clust_pos{cc+1,7}*250)+26,condition),2); %elaborated but proper solution to get time-points (in time-samples) for cluster cc        
        NEG(:,cc,ii) = mean(OUT.sources_ERFs(:,round(MAG_clust_neg{cc+1,6}*250)+26:round(MAG_clust_neg{cc+1,7}*250)+26,condition),2); %elaborated but proper solution to get time-points (in time-samples) for cluster cc
    end
end
POS2 = mean(POS,3); %mean over subjects
NEG2 = mean(NEG,3); %mean over subjects
if condition == 1
    maxp = abs(max(max(POS2(:,1:5)))); %maximum value (in absolute terms) for first 5 clusters (the one reported in the manuscript
    maxn = abs(max(max(NEG2(:,1)))); %same concept (only 1 cluster here reported in the manuscript)
    POS2 = POS2./maxp;
    NEG2 = NEG2./maxn;
elseif condition == 2
    maxp = abs(max(max(POS2(:,1)))); %maximum value (in absolute terms) for first 5 clusters (the one reported in the manuscript
    maxn = abs(max(max(NEG2(:,1)))); %same concept (only 1 cluster here reported in the manuscript)
    POS2 = POS2./maxp;
    NEG2 = NEG2./maxn;
elseif condition == 3
    maxp = abs(max(max(POS2(:,1:3)))); %maximum value (in absolute terms) for first 5 clusters (the one reported in the manuscript
    maxn = abs(max(max(NEG2(:,1)))); %same concept (only 1 cluster here reported in the manuscript)
    POS2 = POS2./maxp;
    NEG2 = NEG2./maxn;
end
    
% 2) plotting main effects and contrasts in brain templates
%POS (young > elderly)
%FROM LBPD COORDINATES TO 3D NIFTI IMAGE
for ii = 1:clust_num
    S = [];
    S.singleimage = 1;
    S.data = POS2(:,ii); %data (voxels x ROIs (time-points))
    S.fname = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/MainEffect_AllParticipants'; %path and name for the image to be saved
    S.names = {['Cond' num2str(condition) '_clust' num2str(ii) '_POS']}; %names for the different images
    %actual function
    FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);
end
%NEG (elderly > young)
%FROM LBPD COORDINATES TO 3D NIFTI IMAGE
for ii = 1:clust_num
    S = [];
    S.singleimage = 1;
    S.data = NEG2(:,ii); %data (voxels x ROIs (time-points))
    S.fname = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/MainEffect_AllParticipants'; %path and name for the image to be saved
    S.names = {['Cond' num2str(condition) '_clust' num2str(ii) '_NEG']}; %names for the different images
    %actual function
    FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);
end

%% Extracting information about the clusters at source level and reporting it in xlsx files

%here we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

path = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/MCS/ForPaper/MainEffect_AllParticipants';
list = dir([path '/*gz']);
for ii = 1:length(list)
    fname = [path '/' list(ii).name]; %tone x cluster 1
    S = [];
    S.fname = fname;
    %actual function
    PDn = Extract_BrainCluster_Information_3D_LBPD_D(S);
    writetable(PDn,[path '/' list(ii).name(1:end-7) '.xlsx'],'Sheet',1); %printing excel file
end

%% MCS on MEG source results

t_val_thresh = 3; %t-value threshold for binarizing the data (referred to contrast young vs elderly)
% = 2 for YOUNG vs ELDERLY
% = 3 for ELDERLY vs YOUNG

%cluster-based MCS to clean things up (either by removing random voxels automatically or by selecting only the biggest cluster)
%Cluster permutations test
% list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/*gz');
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
MASK = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
clear LIST
LIST{1} = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/TvalPOS_Cond1_YoungVsEld_10Clusts.nii.gz';
LIST{2} = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/TvalNEG_Cond1_YoungVsEld_10Clusts.nii.gz';
LIST{3} = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/TvalPOS_Cond2_YoungVsEld_10Clusts.nii.gz';
LIST{4} = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/TvalNEG_Cond2_YoungVsEld_10Clusts.nii.gz';
LIST{5} = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/TvalPOS_Cond3_YoungVsEld_10Clusts.nii.gz';
LIST{6} = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/TvalNEG_Cond3_YoungVsEld_10Clusts.nii.gz';
NAMES{1} = 'Old_YoungVsElderly_POS'; NAMES{2} = 'Old_YoungVsElderly_NEG'; NAMES{3} = 'NewT1_YoungVsElderly_POS'; NAMES{4} = 'NewT1_YoungVsElderly_NEG'; NAMES{5} = 'NewT3_YoungVsElderly_POS'; NAMES{6} = 'NewT3_YoungVsElderly_NEG';
for ff = 1%:length(LIST)
    T = load_nii(LIST{ff});
    %preparing function for MCS
    for ii = 1%:size(T.img,4) %over time-windows
        %getting MNI coordinates
        %extracting matrix with statistics
        T2 = T.img(:,:,:,ii); %extracting time-window ii
        %mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
        mask = zeros(size(T2,1),size(T2,2),size(T2,3));
        data = mask;
        data(abs(T2)>t_val_thresh) = 1; %binarizing data according to t-value threshold
        mask(MASK.img~=0) = 1; %assigning 1 when you have real brain voxels
        %preparation of information and data for the actual function
        S = [];
        S.T = T2;
        S.outdir = ['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/MCS']; %output path
        S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
        S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
        S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
        S.data = data;
        S.mask = mask; %mask of the brain layout you have your results in
        S.permut = 1000; %number of permutations for Monte Carlo simulation
        S.clustmax = 0; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict).
        S.permthresh = 0.001; %threshold for MCS
        %UNTIL HERE
        S.anal_name = [NAMES{ff} '_MEGSensorClust_' num2str(ii) '_MCS']; %name for the analysis (used to identify and save image and results)
        
        %actual function
        PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
        disp(ii)
    end
end
%after this, we created the final images using Workbench

%% Extracting information about the clusters at source level and reporting it in xlsx files

%here we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

path = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/MEGChannelsMCS_Sources/MCS/ForPaper';
list = dir([path '/*gz']);
for ii = 1:length(list)
    fname = [path '/' list(ii).name]; %tone x cluster 1
    S = [];
    S.fname = fname;
    %actual function
    PDn = Extract_BrainCluster_Information_3D_LBPD_D(S);
    writetable(PDn,[path '/' list(ii).name(1:end-7) '.xlsx'],'Sheet',1); %printing excel file
end

%% PLOTTING MEG SENSOR RESULTS (TIME SERIES)

%% MEG sensors - plotting elderly versus young participants (this allows to plot data from other experimental blocks which are not included in the current manuscript)

condition = 3; %1 = Old; 2 = NewT1; 3 = NewT3
clust_num = 1; %numbers of clusters you want
POS_l = 2; %1 = young > elderly; 2 = elderly > young

%loading results (MCS on MEG sensors)
load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS/MAG_Cond' num2str(condition) '_YoungvsEld_tvalp_05.mat']);
if POS_l == 1
    clustplotd = MAG_clust_pos(clust_num + 1,:);
else
    clustplotd = MAG_clust_neg(clust_num + 1,:);
end
clustplot = clustplotd{:,3};
block = 3; %1 = mmn subtracted; 2 = mmn no subtracted; 3 = recogminor; 4 = auditory sequence of numbers; 5 = visual sequence of numbers
young_elderly = 1; %1 = comparing young and elderly; 2 = comparing elderly (60-68 years old) and elderly (>68 years old); 3 = young, high vs low WM; 4 = eldely, high vs low WM
chans = [];
for ii = 1:size(clustplot,1)
    chans(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
end
S = [];
if block == 1
    load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/Block_2_MMNsubtracted.mat']); %loading the channels forming the main significant cluster across all deviants
else
    load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/Block_' num2str(block) '.mat']); %loading the channels forming the main significant cluster across all deviants
end
data_mat = abs(data_mat); %here, doing absolute value of the data since the contrasts were computed on absolute values
%structure for the function
S.data = data_mat;
if young_elderly == 1
    S.groups = {'young','elderly'}; %Group labels
    gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
    gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,35,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76];
elseif young_elderly == 2
    S.groups = {'60-68','>68','young'}; %Group labels
    gsubj{1} = [3,8,10,11,12,20,26,31,34,36,42,45,47,49,51,52,56,58,59,61,62,68]; %60-68
    gsubj{2} = [4,6,7,16,18,23,27,28,29,30,32,35,38,43,48,66,76]; %>68
elseif young_elderly == 3 %young, high WM vs low WM
    S.groups = {'HWMY','LWMY'}; %Group labels
    gsubj{1} = [2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM
    gsubj{2} = [5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM
elseif young_elderly == 4 %elderly, high WM vs low WM
    S.groups = {'HWME','LWME'}; %Group labels
    gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68]; %high WM
    gsubj{2} = [3,4,6,8,18,20,23,26,27,29,31,34,35,36,38,43,48,49,58,61,62,66]; %low WM
end
S.gsubj = gsubj;
S.time_real = time_sel;
S.legendl = 0;
S.colorline = {'r', 'b', 'k'}; %Select colorline for each group
S.sensors = -1; % -1 for magnetometer, 0 for gradiometers (gradiometers have already been combined)
if block == 1
    S.conditions = {'Glob_mmn','Loc_mmn'};
    S.x_lim = [-0.05 1.2]; % Set x limits
    S.y_lim = [-100 100]; %Set y limits
elseif block == 2
    S.conditions = {'Glob_Stand','Glob_Dev','Loc_Stand','Loc_Dev'};
    %     S.conditions = {'Glob_Stand','Glob_Dev','Loc_Stand','Loc_Dev','100'};
    S.x_lim = [-0.05 1.2]; % Set x limits
    S.y_lim = [-180 180]; %Set y limits
elseif block == 3
    S.conditions = {'Old_Correct','New_T1_Correct','New_T3_Correct'};
    %     S.conditions = {'Old_Correct','New_T1_Correct'};
    S.y_lim = [-200 200]; %Set y limits
    S.x_lim = [-0.1 3.4]; % Set x limits
elseif block > 3
    S.conditions = {'Old_Correct','New_T1_Correct','New_T3_Correct'};
    %     S.conditions = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
    S.y_lim = [-200 200]; %Set y limits
    S.x_lim = [-0.1 3.4]; % Set x limits
end
S.y_lim = []; %Set y limits
S.x_lim = [-0.1 3.4]; % Set x limits
S.condition_n = condition;
S.chanlabels = chanlabels;
S.STE = 2;
S.signtp = [clustplotd{1,6} clustplotd{1,7}];
S.chans_index = chans; %plot waveform at single channel or the average across multiple channels (specify channel index - you can find the channel labels in 'chanlabels'); set to 0 to plot all channels

plot_sensors_wavebis2(S) %actual function

%% EXTRACTING DATA TO BE EXPORTED AND PLOTTED LOCALLY

POS_l = 2; %1 = young > elderly; 2 = elderly > young
condition = 3; %1 = Old; 2 = NewT1; 3 = NewT3
clust_num = 1; %numbers of clusters you want

load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/Block_3.mat']); %loading the channels forming the main significant cluster across all deviants
data_mat = abs(data_mat); %here, doing absolute value of the data since the contrasts were computed on absolute values
%loading results (MCS on MEG sensors)
load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/MEG_sensors/MCS/MAG_Cond' num2str(condition) '_YoungvsEld_tvalp_05.mat']);
if POS_l == 1 %young > elderly
    clustplotd = MAG_clust_pos(clust_num + 1,:);
else %elderly > young
    clustplotd = MAG_clust_neg(clust_num + 1,:);
end
%extracting significant channels for requested cluster 
clustplot = clustplotd{:,3};
chans = [];
for ii = 1:size(clustplot,1)
    chans(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
end
%averaging over significant channels
data2 = mean(data_mat(chans,:,:,condition),1);
%extracting significant time-points
signn = [clustplotd{1,6} clustplotd{1,7}];
%saving on disk data prepared for plotting
save(['/aux/MINDLAB2021_MEG-TempSeqAges/06_04_2023/Cond_' num2str(condition) '_clust_' num2str(clust_num) '_POS_l_' num2str(POS_l)],'data2','time_sel','chanlabels','signn')

%%

%% *** GETTING NEW ROIs FROM AAL (BASED ON THE KEY BRAIN REGIONS EMERGED FROM THE PREVIOUS CONTRASTS) ***

%%% OBS! This was done as an exploratory procedure. The results were good
%%% but I decided to use functional ROIs derived from one of our previous
%%% studies. Both solutions returned similar results and in this way I was
%%% more coherent with our previous works. To allow full disclosure, I
%%% decided not to remove the few sections below, even if their results are
%%% not reported in the paper.

% The ROIs are 8 (numbers refers to AAL ROIs (non-symmetric)):
% 1)L auditory cortex: 79,81,83,85
% 2)R auditory cortex: 80,82,84,86
% 3)L hippocampal area and inferior temporal cortex: 87,89,37,39,41
% 4)R hippocampal area and inferior temporal cortex: 88,90,38,40,42
% 5)L inferior frontal gyrus (and operculum and insula): 11,13,15,17,29
% 6)R inferior frontal gyrus (and operculum and insula): 12,14,16,18,30
% 7)Ventro-medial prefrontal cortex (L and R): 31,32,27,28
% 8)Medial cingulate (L and R): 33,34,35,36,67,68


%% CREATING INDEPENDENT FILES FOR EACH AAL REGION

%%% THIS SHOULD BE RUN ONLY ONCE AND THEN THE MASKS CAN BE USED FOR DIFERENT PROJECTS %%%

%creating indepdent image files with AAL ROIs
%command to split one ROI (80) from AAL
% cmd = 'fslmaths aal_1mm.nii.gz -thr 80 -uthr 80 -bin 80.nii.gz';
% system(cmd)

for ii = 1:90
    cmd = ['fslmaths /scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/aal_8mm_try5.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/' num2str(ii) '.nii.gz'];
    system(cmd)
end

%% JOINING THE ROIs IN AAL SPACE TO GET THE FINAL ROIs

path = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm';
%1)LAC
name1 = '79'; name2 = '81'; name3 = '83'; name4 = '85';
output = 'LAC';
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' -add ' path '/' name4 ' ' path  '/FinalROIs/' output]; %OBS! be careful with the spacing
system(cmd)
%2)RAC
name1 = '80'; name2 = '82'; name3 = '84'; name4 = '86';
output = 'RAC';
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' -add ' path '/' name4 ' ' path  '/FinalROIs/' output]; %OBS! be careful with the spacing
system(cmd)
%3)LHIT
name1 = '87'; name2 = '89'; name3 = '37'; name4 = '39'; name5 = '41';
output = 'LHIT';
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' -add ' path '/' name4 ' -add ' path '/' name5 ' ' path  '/FinalROIs/' output]; %OBS! be careful with the spacing
system(cmd)
%4)RHIT
name1 = '88'; name2 = '90'; name3 = '38'; name4 = '40'; name5 = '42';
output = 'RHIT';
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' -add ' path '/' name4 ' -add ' path '/' name5 ' ' path  '/FinalROIs/' output]; %OBS! be careful with the spacing
system(cmd)
%5)RIFG is taken from here:
copyfile /scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband_invers_1/IFG_young_vs_elderly_Old_CorrectTval.nii.gz /scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/FinalROIs/RIFG.nii.gz
%6)LIFG (mirror image to RIFG)
list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/FinalROIs/RIF*gz');
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat'); %loading new cooridnate system..
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
ROIs = zeros(3559,length(list)); %initializing ROIs
for mm = 1:length(list) %over ROIs
    roi = load_nii([list(mm).folder '/' list(mm).name]); %loading parcel (ROI) mm
    for ii = 1:size(maskk.img,1) %over dimension 1 of the 3D image
        for jj = 1:size(maskk.img,2) %over dimension 2 of the 3D image
            for zz = 1:size(maskk.img,3) %over dimension 3 of the 3D image
                if roi.img(ii,jj,zz) ~= 0 %if there is a voxel forming the ROI
                    ROIs(maskk.img(ii,jj,zz),mm) = 1; %assigning one to that voxel in ROIs (by taking number of voxel from the maskk.img)
                end
            end
        end
    end
end
%LIFG
ab = find(ROIs(:,1)==1);
LIFG = zeros(3559,1);
for ii = 1:length(ab) %over voxels of HITR
    %[MIN,IDX] = min(abs(sum(MNI8-[MNI8(ab(ii),1)*(-1),MNI8(ab(ii),2),MNI8(ab(ii),3)],2)))
    dummy = zeros(3559,3);
    dummy((MNI8(ab(ii),1)*(-1)+4)==MNI8(:,1),1) = 1; %getting values in the other hemisphere (multiplying by -1 o reverse the sign and adding 4 since there is a very small (negligible) imperfection of the coordinate system (symmetric across hemispheres voxels are slightly misalligned.. nothing that really matters here in MEG)
    dummy((MNI8(ab(ii),2)==MNI8(:,2)),2) = 1;
    dummy((MNI8(ab(ii),3)==MNI8(:,3)),3) = 1;
    LIFG((sum(dummy,2)==3),1) = 1;
end
%from LBPD coordinates to nifti image
S.data = LIFG; %data (voxels x ROIs (time-points))
S.fname = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/FinalROIs'; %path and name for the image to be saved
names{1} = 'LIFG';
S.names = names; %names for the different images
%actual function
FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);
%7)VMPFC
name1 = '31'; name2 = '32'; name3 = '27'; name4 = '28';
output = 'VMPFC';
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' -add ' path '/' name4 ' ' path  '/FinalROIs/' output]; %OBS! be careful with the spacing
system(cmd)
%8)MC
name1 = '33'; name2 = '34'; name3 = '35'; name4 = '36'; name5 = '67'; name6 = '68';
output = 'MC';
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' -add ' path '/' name3 ' -add ' path '/' name4 ' -add ' path '/' name5 ' -add ' path '/' name6 ' ' path  '/FinalROIs/' output]; %OBS! be careful with the spacing
system(cmd)


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
list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
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
    save(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat'],'dum2','condds')
end

%SUBJ0035 HAS NO DATA FOR NEWT1 AND NEWT3!!

%%

%% FROM HERE THE CODES REFLECT AGAIN THE PROCEDURES REPORTED IN THE PAPER

%%

%% *** EXTRACTING MNI COORDINATES FOR THE SELECTED FUNCTIONAL ROIs - TO REPORT THEM IN A SUPPLEMENTARY TABLE ***

%%

list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/SingleImages_8mm/FinalROIs/ROIs_SystVarPaper_PlusIFG/*gz');
outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/09_03_2023';

for ii = 1:length(list) %over ROIs
    CELL = cell(size(mni_coords,1)+2,4);
    %getting MNI coordinates of significant voxels within the provided image
    [ mni_coords, ~ ] = osl_mnimask2mnicoords([list(ii).folder '/' list(ii).name]);
    CELL{1,1} = [list(ii).name(1:end-7) ' voxel #'];
    CELL{1,2} = 'MNI coordinates';
    CELL{2,2} = 'x'; CELL{2,3} = 'y'; CELL{2,4} = 'z';
    numm = num2cell(1:size(mni_coords,1)); %progressive number of the brain voxels
    CELL(3:size(mni_coords,1)+2,1) = numm;
    CELL(3:size(mni_coords,1)+2,2) = num2cell(mni_coords(:,1)); %x
    CELL(3:size(mni_coords,1)+2,3) = num2cell(mni_coords(:,2)); %y
    CELL(3:size(mni_coords,1)+2,4) = num2cell(mni_coords(:,3)); %z
    PDn = cell2table(CELL); %table
    writetable(PDn,[outdir '/ROIs_MNI_Coordinates.xlsx'],'Sheet',ii); %printing excel file
end

%%

%% *** COMPUTING THE STATISTICS (ELDERLY VS YOUNG) FOR THE TIME SERIES OF THOSE ROIs (INDPENDENTLY FOR THE THREE EXPERIMENTAL CONDITIONS) ***

%%

%SUBJ0035 HAS NO DATA FOR NEWT1 AND NEWT3!! MUST BE EXCLUDED BY THIS ANALYSIS

ROIs_l = 2; %1 = AAL ROIs; 2 = systematic variation recognition paper ROIs

gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78]; %young
gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76]; %elderly
%this is the file that you have to load and use to run the statistics
if ROIs_l == 1
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8.mat');
else
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat');
end
clear ROIN
ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
%t-tests
P = zeros(size(dum2,1),size(dum2,2),size(dum2,3)); %ROIs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dum2,1),size(dum2,2),size(dum2,3)); %ROIs x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %ove time-points
        for cc = 1:size(dum2,3) %over experimental conditions
            %codes t-test
            a = squeeze(dum2(ii,jj,cc,gsubj{1})); %young
            b = squeeze(dum2(ii,jj,cc,gsubj{2})); %elderly
            ad = find(isnan(a)); %this is for preventing t-tests to be calculated considering NaNs.. if there are any subject in a or in b corresponding to NaN both a and b are deprivated of that subject
            ab = find(isnan(b));
            abd = [ab ad];
            a(abd) = [];
            b(abd) = [];
            [~,pttest,~,tstats] = ttest2(a,b);
            P(ii,jj,cc) = pttest;
            T(ii,jj,cc) = tstats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
if ROIs_l == 1
    outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics';
else
    outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics/SystPaperROIs';
    mkdir(outdir)
end
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:size(dum2,3) %over experimental conditions
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        save([outdir '/' ROIN{ii} '_Cond' num2str(cc) '.mat'],'PDn'); %printing excel file
        writetable(PDn,[outdir '/' ROIN{ii} '_Cond' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
    end
end

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER

for qq = 1:8
    ROIs_l = 2; %1 = AAL ROIs; 2 = systematic variation recognition paper ROIs
    ROI = [qq]; %1 = 'LAC'; 2 = 'LHIT'; 3 = 'LIFG'; 4 = 'MC'; 5 = 'RAC'; 6 = 'RHIT'; 7 = 'RIFG'; 8 = 'VMPFC'; % ROIs in dum2; (UNLESS IT'S APR2020 OR SEQUENCES OF NUMBERS)
    condition = [4]; %1 = Old; 2 = NewT1; 3 = NewT3; (UNLESS IT'S APR2020 OR SEQUENCES OF NUMBERS)
    if ROI == 1 || ROI == 5 %LAC anord RAC
        ylimm = [-100 40]; %amplitude limits; leave empty [] for automatic adjustment
    else %the other ROIs
        ylimm = [-60 30]; %amplitude limits; leave empty [] for automatic adjustment
    end
    eld_l = 7; %1 = comparing young and elderly; 2 = comparing elderly (60-68 years old) and elderly (>68 years old);
               %3 = young, high vs low WM; 4 = eldely, high vs low WM; 5 = high vs low WM (all participants together);
               %6 all participants toegheter high and low WM - data from APR2020; 7 = young/elderly - sequences of numbers
               %-1 = high WM elderly, high WM young, low WM elderly, low WM young 
    export_l = 0; %1 = export images; 0 = not
    
               
    clear ROIN
    S = [];
    if eld_l == 6 %all participants high and low WM - APR2020
        load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/ROIs_6.mat');
        ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';
        S.conds = {'Old Sys','NewT1 Sys','NewT2 Sys','NewT3 Sys','NewT4 Sys'};
    elseif eld_l == 7 %sequences of numbers
        load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
        S.conds = {'Old Sys','NewT1 Sys','NewT2 Sys','NewT3 Sys','NewT4 Sys'};
        ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';
    elseif eld_l < 6
        if ROIs_l == 1 %1 = AAL ROIs; 2 = systematic variation recognition paper ROIs
            S.conds = {'Old AAL','NewT1 AAL','NewT3 AAL'};
        else
            S.conds = {'Old Sys','NewT1 Sys','NewT3 Sys'};
        end
        ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
        if ROIs_l == 1 %1 = AAL ROIs; 2 = systematic variation recognition paper ROIs
            load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8.mat');
        else
            load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat');
        end
    end
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
    if eld_l < 6
        signtp = cell(length(ROI),length(condition));
        for ii = 1:length(ROI)
            for cc = 1:length(condition)
                if ROIs_l == 1
                    load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics/' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.mat'])
                else
                    load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics/SystPaperROIs/' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.mat'])
                end
                signtp{ii,cc} = table2array(PDn(:,3)); %4D matrix (clusters x time (begin and end) x ROIs x conditions)
            end
        end
    end
    dum2 = permute(dum2,[1 2 4 3]);
    dum3 = dum2(:,:,:,:);
    
    %structure for the function
    S.data = dum3;
    if eld_l == 1
        S.groups = {'young','elderly'}; %Group labels
        gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
        gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76];
    elseif eld_l == 2
        S.groups = {'60-68','>68','young'}; %Group labels
        gsubj{1} = [3,8,10,11,12,20,26,31,34,36,42,45,47,49,51,52,56,58,59,61,62,68]; %60-68
        gsubj{2} = [4,6,7,16,18,23,27,28,29,30,32,38,43,48,66,76]; %>68
        gsubj{3} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
    elseif eld_l == 3 %young, high WM vs low WM
        S.groups = {'HWMY','LWMY'}; %Group labels
        gsubj{1} = [2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM
        gsubj{2} = [5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM
    elseif eld_l == 4 %elderly, high WM vs low WM
        S.groups = {'HWME','LWME'}; %Group labels
        gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68]; %high WM
        gsubj{2} = [3,4,6,8,18,20,23,26,27,29,31,34,36,38,43,48,49,58,61,62,66]; %low WM
    elseif eld_l == 5 %all participants (high and low WM)
        S.groups = {'HWMEY','LWMEY'}; %Group labels
        gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68,2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM
        gsubj{2} = [3,4,6,8,18,20,23,26,27,29,31,34,36,38,43,48,49,58,61,62,66,5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM
    elseif eld_l == 6 %all participants (high and low WM) - APR2020
        [numWM,~,WM] = xlsread('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/AuditoryPatternRecognition/behavioral_data/WM_APR2020_Edited_SR.xlsx'); %loads working memory data
        iop = nanmean(numWM(:,2)); %mean WM
        gsubj{1} = find(numWM(:,2)>iop); %high WM
        gsubj{2} = find(numWM(:,2)<iop); %low WM
        S.groups = {'HWM APR','LWM APR'}; %Group labels
    elseif eld_l == 7 %sequences of numbers
        S.groups = {'young','elderly'}; %Group labels
        gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
        gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76];
    elseif eld_l == -1
        S.groups = {'WMHE','WMHY','WMLE','WMLY'}; %Group labels
        gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68]; %high WM elderly
        gsubj{2} = [2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM young
        gsubj{3} = [3,4,6,8,18,20,23,26,27,29,31,34,36,38,43,48,49,58,61,62,66]; %low WM elderly
        gsubj{4} = [5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM young
    end
    S.gsubj = gsubj;
    S.time_real = time(1:1026);
    S.legendl = 1;
    dumcol = colormap(lines(length(condition) * length(S.groups)* length(ROI))); %extracting some colours from a colormap
    if eld_l == 1 %if young vs elderly we want young = red and elderly = blue
        S.colorline = zeros(size(dumcol,1),3); S.colorline(1,:) = dumcol(2,:); S.colorline(2,:) = dumcol(1,:); %swapping colors for young and elderly
    else
        S.colorline = dumcol;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = ROI;
    S.condition_n = condition;
    S.ROIs_labels = ROIN(ROI);
    if eld_l > 1 || eld_l == -1
        S.signtp = [];
    else
        S.signtp = signtp;
    end
    
    waveplot_groups_server_LBPD(S) %actual function
    if export_l == 1
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/03_03_2023/ROIs_Images/' ROIN{ROI} '_Cond' num2str(condition(cc)) '.png'])
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/03_03_2023/ROIs_Images/' ROIN{ROI} '_Cond' num2str(condition(cc)) '.eps'])
    end
end
    
%% EXTRACTING DATA FOR RUNNING LINEAR MODELS IN R

[num,~,raw] = xlsread('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/TempSeqAges2021_WM_Education_MusicalTraining.xlsx');
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/03_03_2023';
timeold = [0.3 0.7 1.1 1.45 1.8; 0.2 0.6 0.95 1.35 1.7; 0.28 0.7 1.1 1.45 1.8; 0.1 0.55 0.9 1.28 1.65; 0.1 0.55 0.9 1.28 1.65; 0.3 0.65 1.1 1.45 1.8; 0.28 0.75 1.15 1.45 1.85; 0.1 0.55 0.9 1.28 1.65];
timenewt1 = [0.62; 0.62; 0.62; 0.52; 0.52; 0.62; 0.62; 0.52];
timenewt3 = [1.37; 1.37; 1.37; 1.28; 1.28; 1.37; 1.37; 1.28];
times{1} = timeold; times{2} = timenewt1; times{3} = timenewt3;
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat');
datt = cell(79,13);
datt{1,1} = 'Subject';
datt{1,2} = 'LAC'; datt{1,3} = 'LHIT'; datt{1,4} = 'LIFG'; datt{1,5} = 'MC'; datt{1,6} = 'RAC'; datt{1,7} = 'RHIT'; datt{1,8} = 'RIFG'; datt{1,9} = 'VMPFC'; % ROIs in dum2
datt{1,10} = 'age'; datt{1,11} = 'years education'; datt{1,12} = 'WM'; datt{1,13} = 'years musical training'; datt{1,14} = 'sex';
for cc = 1:length(times) %over conditions
    timedum = times{cc}; %extracting time-points for condition cc
    for ii = 1:size(dum2,4) %over subjects
        datt{ii+1,1} = ii;
        if ismember(ii,[3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76])
            datt{ii+1,10} = 2; %2 = elderly
        elseif ismember(ii,[2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78])
            datt{ii+1,10} = 1; %1 = young
        else
            datt{ii+1,10} = 'NaN'; %neither young nor elderly
        end
        datt(ii+1,11) = raw(ii+1,5); %years of education
        datt(ii+1,12) = raw(ii+1,7); %WM
        datt(ii+1,13) = raw(ii+1,6); %years of musical training
        if strcmp(raw{ii+1,4},'M') %if male
            datt{ii+1,14} = 1; %male = 1
        elseif strcmp(raw{ii+1,4},'F') %if female
            datt{ii+1,14} = 2; %female = 1
        end
        for ss = 1:size(timedum,1) %over ROIs
            dumt = zeros(1,length(timedum));
            for tt = 1:size(timedum,2) %over selected time-points
                [a,ai] = min(abs(time-timedum(ss,tt)));
                dumt(1,tt) = squeeze(mean(dum2(ss,ai-3:ai+3,cc,ii),2));
            end
            datt{ii+1,ss+1} = mean(dumt);
        end
    end
    PDn = cell2table(datt); %table
    writetable(PDn,[outdir '/R_Old_NewT1_NewT3.xlsx'],'Sheet',cc); %printing excel file
    disp(cc)
end

% please, see the R script for additional details in these analysis

%%

%% *** ANOVA FOR THREE DIFFERENT GROUPS (YOUNG, ELDERLY 60-68 AND ELDERLY OLDER THAN 68) ***

%%

%SUBJ0035 HAS NO DATA FOR NEWT1 AND NEWT3!! MUST BE EXCLUDED BY THIS ANALYSIS

groups = {'young','60-68','>68'}; %Group labels
gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78]; %young
gsubj{2} = [3,8,10,11,12,20,26,31,33,34,36,42,45,47,49,51,52,56,58,59,61,62,68]; %60-68
gsubj{3} = [4,6,7,16,18,23,27,28,29,30,32,38,43,48,66,76]; %>68
%building cell array groupp with label corresponding to the different age groups
groupp = cell(76,1); %hard coding because I know there are 76 participants
cnt = 0;
for pp = 1:length(gsubj{1}) %over young
    cnt = cnt + 1;
    groupp{cnt} = 'y';
end
for pp = 1:length(gsubj{2}) %over elderly (60-68)
    cnt = cnt + 1;
    groupp{cnt} = 'e';
end
for pp = 1:length(gsubj{3}) %over elderly (>68)
    cnt = cnt + 1;
    groupp{cnt} = 'ee';
end
%this is the file that you have to load and use to run the statistics
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat');
clear ROIN
ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
%t-tests
P = zeros(size(dum2,1),size(dum2,2),size(dum2,3)); %ROIs x time-points x conditions
F = zeros(size(dum2,1),size(dum2,2),size(dum2,3)); %ROIs x time-points x conditions
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %ove time-points
        for cc = 1:size(dum2,3) %over experimental conditions
            datadum = [];
            %preparing data for ANOVA
            datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc,gsubj{1}))); %extracting data for group 1
            datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc,gsubj{2}))); %group 2
            datadum = cat(1,datadum,squeeze(dum2(ii,jj,cc,gsubj{3}))); %group 3
            ab = find(isnan(datadum)); %checking for NaNs
            if ~isempty(ab) %if there are NaNs
                datadum(ab) = []; %removing them
                groupp2 = groupp;
                groupp2(ab) = []; %and removing corresponding labels
            else
                groupp2 = groupp;
            end
            %ANOVA
            [p,f,~] = anova1(datadum,groupp,'off'); %'off' for not showing the plot
            P(ii,jj,cc) = p;
            F(ii,jj,cc) = f{2,5};
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics/ThreeGroupsAge';
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:size(dum2,3) %over experimental conditions
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        F_sel = F(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), F_sel ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        save([outdir '/' ROIN{ii} '_Cond' num2str(cc) '.mat'],'PDn'); %printing excel file
        writetable(PDn,[outdir '/' ROIN{ii} '_Cond' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
    end
end

%%

%% *** TWO-WAY ANOVA TESTING YOUNG AND ELDERLY AND HIGH AND LOW WM ***

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
P = zeros(size(dum2,1),size(dum2,2),size(dum2,3)); %ROIs x time-points x conditions
F = zeros(size(dum2,1),size(dum2,2),size(dum2,3)); %ROIs x time-points x conditions
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
            P(ii,jj,cc) = min(p); %getting smallest p-value of the three (the idea is that any significance (either WM, age or the interaction) is relevant for this analysis 
            F(ii,jj,cc) = f{2,6};
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics/WM_Age';
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:size(dum2,3) %over experimental conditions
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        F_sel = F(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), F_sel ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        save([outdir '/' ROIN{ii} '_Cond' num2str(cc) '.mat'],'PDn'); %printing excel file
        writetable(PDn,[outdir '/' ROIN{ii} '_Cond' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
    end
end

%%

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 3); %slot in the queu

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER (THREE AGE-RELATED GROUPS OR WM-AGE TOGETHER)

for qq = 1
    ROI = [qq]; %1 = 'LAC'; 2 = 'LHIT'; 3 = 'LIFG'; 4 = 'MC'; 5 = 'RAC'; 6 = 'RHIT'; 7 = 'RIFG'; 8 = 'VMPFC'; % ROIs in dum2; (UNLESS IT'S APR2020 OR SEQUENCES OF NUMBERS)
    condition = [1]; %1 = Old; 2 = NewT1; 3 = NewT3; (UNLESS IT'S APR2020 OR SEQUENCES OF NUMBERS)
    if ROI == 1 || ROI == 5 %LAC anord RAC
        ylimm = [-100 40]; %amplitude limits; leave empty [] for automatic adjustment
    else %the other ROIs
        ylimm = [-60 30]; %amplitude limits; leave empty [] for automatic adjustment
    end
    eld_l = 1; %1 = comparing elderly (60-68 years old) and elderly (>68 years old);
               %2 = high WM elderly, high WM young, low WM elderly, low WM young 
    export_l = 0; %1 = export images; 0 = not
    
    clear ROIN
    S = [];
    S.conds = {'Old Sys','NewT1 Sys','NewT3 Sys'};
    ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/ROIs_8_SystPaper.mat');
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
    %getting significant time-points
    signtp = cell(length(ROI),length(condition));
    for ii = 1:length(ROI)
        for cc = 1:length(condition)
            if eld_l == 1
                load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics/ThreeGroupsAge/' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.mat'])
            else
                load(['/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/ROIs_LearningBach_ElderlyvsYoung/AAL_ROIs/Statistics/WM_Age/' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.mat'])
            end
            signtp{ii,cc} = table2array(PDn(:,3)); %4D matrix (clusters x time (begin and end) x ROIs x conditions)
        end
    end
    dum2 = permute(dum2,[1 2 4 3]); %reshaping order of the dimensions in the data
    dum3 = dum2(:,:,:,:);
    %structure for the function
    S.data = dum3;
    if eld_l == 1
        S.groups = {'60-68','>68','young'}; %Group labels
        gsubj{1} = [3,8,10,11,12,20,26,31,33,34,36,42,45,47,49,51,52,56,58,59,61,62,68]; %60-68
        gsubj{2} = [4,6,7,16,18,23,27,28,29,30,32,38,43,48,66,76]; %>68
        gsubj{3} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78]; %young
    elseif eld_l == 2 %young, high WM vs low WM
        S.groups = {'WMHE','WMHY','WMLE','WMLY'}; %Group labels
        gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68]; %high WM elderly
        gsubj{2} = [2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM young
        gsubj{3} = [3,4,6,8,18,20,23,26,27,29,31,34,36,38,43,48,49,58,61,62,66]; %low WM elderly
        gsubj{4} = [5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM young
    end
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.gsubj = gsubj;
    S.time_real = time(1:1026);

    if eld_l == 1
        S.colorline = [0 0 0.502; 0 0.55 1; 0.852 0 0];
    else
        S.colorline = [0 0 0.502; 0.702 0 0; 0 0.55 1; 1 0.2 0.2];
    end
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = ROI;
    S.condition_n = condition;
    S.ROIs_labels = ROIN(ROI);
    S.signtp = signtp;
    
%     jobid = job2cluster(@waveplot_groups_server_LBPD,S);

    waveplot_groups_server_LBPD(S) %actual function
    
    if eld_l == 1
        barb = 'age3';
    else
        barb = 'WM_age';
    end
    if export_l == 1
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/15_03_2023/' barb '_' ROIN{ROI} '_Cond' num2str(condition(cc)) '.png'])
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/15_03_2023/' barb '_' ROIN{ROI} '_Cond' num2str(condition(cc)) '.eps'])
    end
end
   
%%

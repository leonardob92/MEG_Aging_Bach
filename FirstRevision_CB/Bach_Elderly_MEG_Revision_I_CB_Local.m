%%

%% CHALLENGING AGE-RELATED DECLINE IN BRAIN FUNCTION: EVIDENCE FROM FAST NEUROIMAGING OF MUSICAL SEQUENCE RECOGNITION - TEMPSEQAGES 2021 - LEONARDO BONETTI

%%

%% LOCAL COMPUTER - REVISION I - COMMUNICATIONS BIOLOGY

%This code was used locally to prepare plots with higher resolution compared to the solutions offered by the Aarhus' cluster of computers

%%

%% IN-KEY AND OUT-OF-KEY "NEW" MELODIES

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
addpath('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Codes')

clear
%loading data, labels, significant time-windows, etc.
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/ROIs_8_SystPaper_Rev_I.mat') %data
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/data_corr.mat','ROIN') %ROIs label
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/time_normal.mat'); %loading time

%% actual code for plotting

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
% ylimm = []; %amplitude limits; leave empty [] for automatic adjustment

ROI = 8; %ROIs (they are 6)
condition = [2:3]; %look in "condds"
export_l = 1; %1 = export images; 0 = not
if ROI == 1 || ROI == 5
    ylimm = [-100 40]; %amplitude limits; leave empty [] for automatic adjustment
else
    ylimm = [-60 30]; %amplitude limits; leave empty [] for automatic adjustment
end
cnt = 0;
close all
for ii = 1%:size(ROIs_to_AAL,1) %over original parcels (+1 which is the voxels which did not belong to any parcel)
    signtp_col = [];
    signtp = [];
    if sum(condition) == 5 %barbaric trick to know what conditions user requested (this is almost hard-coding and does not generalise well to other contexts)
        bim = dir(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/Stats_Keys_Melodies/NewT1*' ROIN{ROI(ii)} '_*.mat']);
    elseif sum(condition) == 9
        bim = dir(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/Stats_Keys_Melodies/NewT3*' ROIN{ROI(ii)} '_*.mat']);
    end
    for ss = 1:length(bim) %over main/interaction effects of ANOVA
        load([bim(ss).folder '/' bim(ss).name]) %loading significant clusters
        sbum = table2cell(PDn); %converting the table into a cell
        sbom = sbum(:,3); %extracting time-windows of significant clusters
        if sum(double(cellfun(@min,sbom)<0.0)+double(cellfun(@max,sbom)>3)) > 0
            sbom(find(double(cellfun(@min,sbom)<0.35)+double(cellfun(@max,sbom)>3))) = []; %removing extremes for plotting purposes
        end
        signtp = cat(1,signtp,sbom); %concatenating significant clusters over the main/interaction effects of ANOVA
        dumss = repmat(ss,length(sbom),1); %number of clusters repeated with effect ID (ss)
        signtp_col = cat(1,signtp_col,dumss); %concatenating them
    end
    lineplot = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    cnt = cnt + 1;
    S = [];
    S.ii = ii;
    S.conds = condds;
    %structure for the function
    data2 = dum2(:,:,:,:);
    S.data = permute(data2,[1 2 4 3]);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.groups = {'young','elderly'}; %Group labels
    gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
    gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76];
    S.gsubj = gsubj;
    S.time_real = time(1:1026);
    
%     S.colorline = [0 0 0.502; 0.702 0 0; 0 0.55 1; 1 0.2 0.2]; %colors for the waveforms
    
    S.colorline = [0.702 0 0; 0 0 0.502; 1 0.2 0.2; 0 0.55 1]; %colors for the waveforms

    
%     S.colorsign = {'-';'--';':'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
    S.colorsign = {'o';'+';'*'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
    
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = ROI;
    S.condition_n = condition;
    S.ROIs_labels = ROIN(ROI(ii));
    S.lineplot = lineplot(cnt,:);
    %         S.subplot = [];
    
    S.signtp = signtp;
    if col_l == 1
        S.signtp_col = signtp_col;
    else
        S.signtp_col = [];
    end
    
    waveplot_groups_local_v2(S) %actual function
    %         title(lab(ROIs_to_AAL{ii,1}(pp),:))
    
    if export_l == 1
        if sum(condition) == 5 %barbaric trick to know what conditions user requested (this is almost hard-coding and does not generalise well to other contexts)
            sbrim = 'NewT1';
        elseif sum(condition) == 9
            sbrim = 'NewT3';
        end
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Francesco_Key/' sbrim '_ROI_' ROIN{ROI} '.pdf'],'Resolution',300)
    end
end

%%

%% WM AND AGE - RECOMPUTED WITH BETTER PLOTTING AND ANALYSISES (MAIN/INTERACTION EFFECTS FROM ANOVAs KEPT SEPARATE)

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
addpath('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Codes')

clear
%loading data, labels, significant time-windows, etc.
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/ROIs_8_SystPaper.mat'); %loading data
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/time_normal.mat'); %loading time
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/data_corr.mat','ROIN') %ROIs label

%% actual code for plotting

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = [-70 30]; %amplitude limits; leave empty [] for automatic adjustment

ROI = 2; %ROIs (they are 6)
condition = 3; %1 = Old; 2 = NewT1; 3 = NewT3;
export_l = 1; %1 = export images; 0 = not
cnt = 0;
close all
for ii = 1%:size(ROIs_to_AAL,1) %over original parcels (+1 which is the voxels which did not belong to any parcel)
    signtp_col = [];
    signtp = [];
    bim = dir(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/Stats_WM_Updated/NewT3*' ROIN{ROI(ii)} '*Cond' num2str(condition) '*.mat']); %note that the beginning of the name "NewT3) is wrong and does not relate to "NewT3". Refer to CondX for the actual experimental condition
    for ss = 1:length(bim) %over main/interaction effects of ANOVA
        load([bim(ss).folder '/' bim(ss).name]) %loading significant clusters
        sbum = table2cell(PDn); %converting the table into a cell
        sbom = sbum(:,3); %extracting time-windows of significant clusters
        if sum(double(cellfun(@min,sbom)<0.35)+double(cellfun(@max,sbom)>3)) > 0
            sbom(find(double(cellfun(@min,sbom)<0)+double(cellfun(@max,sbom)>3))) = []; %removing extremes for plotting purposes
        end
        signtp = cat(1,signtp,sbom); %concatenating significant clusters over the main/interaction effects of ANOVA
        dumss = repmat(ss,length(sbom),1); %number of clusters repeated with effect ID (ss)
        signtp_col = cat(1,signtp_col,dumss); %concatenating them
    end
    lineplot = [27 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    cnt = cnt + 1;
    S = [];
    S.ii = ii;
    S.conds = condds;
    %structure for the function
    data2 = dum2(:,:,:,:);
    S.data = permute(data2,[1 2 4 3]);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.groups = {'WMHE','WMHY','WMLE','WMLY'}; %Group labels
    gsubj{1} = [7,10,11,12,16,28,30,32,33,42,45,47,51,52,56,59,68]; %high WM elderly
    gsubj{2} = [2,17,19,22,24,39,41,46,50,53,54,57,60,64,69,72,78]; %high WM young
    gsubj{3} = [3,4,6,8,18,20,23,26,27,29,31,34,36,38,43,48,49,58,61,62,66]; %low WM elderly
    gsubj{4} = [5,9,13,14,15,21,25,37,40,44,55,63,65,70,71,73,74]; %low WM young
    S.gsubj = gsubj;
    S.time_real = time(1:1026);
    
    S.colorline = [0 0 0.502; 0.702 0 0; 0 0.55 1; 1 0.2 0.2]; %colors for the waveforms
%     S.colorsign = {'-';'--';':'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
    S.colorsign = {'o';'+';'*'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
    
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = ROI;
    S.condition_n = condition;
    S.ROIs_labels = ROIN(ROI(ii));
    S.lineplot = lineplot(cnt,:);
    %         S.subplot = [];
    
    S.signtp = signtp;
    if col_l == 1
        S.signtp_col = signtp_col;
    else
        S.signtp_col = [];
    end
    
    waveplot_groups_local_v2(S) %actual function
    %         title(lab(ROIs_to_AAL{ii,1}(pp),:))
    
    if export_l == 1
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Francesco_WM/WM_ROI_' ROIN{ROI} '_Cond' num2str(condition) '.pdf'],'Resolution',300)
    end
    title(condds{condition})
end

%%

%% CORRELATION BETWEEN FAMILIARITY WITH THE BACH'S MUSICAL PIECE AND THE NEURAL DATA

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
addpath('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Codes')

clear
%loading data, labels, significant time-windows, etc.
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/time_normal.mat'); %loading time
load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/data_corr.mat') %ROIs label

%% actual code for plotting

%%% OBS!! HERE I USED A FUNCTION NOT ORIGINALLY THOUGHT FOR THIS PURPOSE BUT VERY HANDY IN THIS CONTEXT
%%% HOWEVER, FOR FUTURE REFERENCES, PERHAPS BETTER NOT TO USE IT AGAIN IN THIS WAY

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = [-0.5 0.5]; %amplitude limits; leave empty [] for automatic adjustment

ROI = 1; %ROIs (they are 6)
condition = 1:3; %1 = Old; 2 = NewT1; 3 = NewT3;
export_l = 0; %1 = export images; 0 = not
cnt = 0;
close all
for ii = 1%:size(ROIN,1) %over original parcels (+1 which is the voxels which did not belong to any parcel)
    signtp_col = [];
    signtp = [];
    cnt2 = 0;
    for cc = 1:length(condition) %over conditions
        %reading significant clusters for each condition
        bim = dir(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Up_Down/MCS_Correlations/' ROIN{ROI(ii)} '*Cond' num2str(condition(cc)) '*.mat']); %note that the beginning of the name "NewT3) is wrong and does not relate to "NewT3". Refer to CondX for the actual experimental condition
        for ss = 1:length(bim) %over main/interaction effects of ANOVA
            cnt2 = cnt2 + 1;
            load([bim(ss).folder '/' bim(ss).name]) %loading significant clusters
            sbum = table2cell(PDn); %converting the table into a cell
            sbom = sbum(:,3); %extracting time-windows of significant clusters
            if sum(double(cellfun(@min,sbom)<0.35)+double(cellfun(@max,sbom)>3)) > 0
                sbom(find(double(cellfun(@min,sbom)<0)+double(cellfun(@max,sbom)>3))) = []; %removing extremes for plotting purposes (still reporting this information in the supplementary tables)
            end
            signtp = cat(1,signtp,sbom); %concatenating significant clusters over the main/interaction effects of ANOVA
            dumss = repmat(cnt2,length(sbom),1); %number of clusters repeated with effect ID (ss)
            signtp_col = cat(1,signtp_col,dumss); %concatenating them
        end
    end
    
    lineplot = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    cnt = cnt + 1;
    S = [];
    S.ii = ii;
    S.conds = condds;
    %structure for the function
    data2 = R(:,:,:,1);
    S.data = permute(data2,[1 2 4 3]);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.groups = {'R'}; %Group labels
    gsubj{1} = [1]; %trick to directly plot the correlation coefficient

    S.gsubj = gsubj;
    S.time_real = time(1:1026);
    S.colorline = [0.702 0 0; 0 0.55 1; 0 0 0.502];%; 0.702 0 0; 0 0.55 1]; %colors for the waveforms
%     bumba = S.colorline(1:length(condition),:); %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
%     S.colorline = bumba; S.colorline(1,:) = bumba(2,:); S.colorline(2,:) = bumba(1,:); %barbaric trick to get the colors that I want
%     S.colorsign = {'-';'--';':'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
%     S.colorsign = {'*'}; %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
    S.colorsign = S.colorline(1:length(condition),:); %colors/different symbols for significant clusters (different colors for different effects of the ANOVA)
    
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = ROI;
    S.condition_n = condition;
    S.ROIs_labels = ROIN(ROI(ii));
    S.lineplot = lineplot(cnt,:);
    %         S.subplot = [];
    
    S.signtp = signtp;
    if col_l == 1
        S.signtp_col = signtp_col;
    else
        S.signtp_col = [];
    end
    
    waveplot_groups_local_v2(S) %actual function
    %         title(lab(ROIs_to_AAL{ii,1}(pp),:))
    
    if export_l == 1
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Manuscript/CommunicationsBiology/FirstRevision/Francesco_CorrFam/WM_ROI_' ROIN{ROI} '_Cond' num2str(condition) '.pdf'],'Resolution',300)
    end
end

%%

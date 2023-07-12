
%% Neural differences between young and older adults in automatic and conscious auditory predictive coding

%%% waveform plotting - computed on a local computer %%%

%%

%% MEG SOURCES

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER (ALL SUBJECTS AT THE SAME TIME)

addpath('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots');

for qq = 1
    ROI = [qq]; %1 = 'LAC'; 2 = 'LHIT'; 3 = 'LIFG'; 4 = 'MC'; 5 = 'RAC'; 6 = 'RHIT'; 7 = 'RIFG'; 8 = 'VMPFC'; % ROIs in dum2;
    condition = [3]; %1 = Old; 2 = NewT1; 3 = NewT3;
    if ROI == 1 || ROI == 5 %LAC anord RAC
        ylimm = [-100 40]; %amplitude limits; leave empty [] for automatic adjustment
    else %the other ROIs
        ylimm = [-60 30]; %amplitude limits; leave empty [] for automatic adjustment
    end    
    export_l = 0; %1 = export images; 0 = not
    
    clear ROIN
    S = [];
    S.conds = {'Old Sys','NewT1 Sys','NewT3 Sys'};
    ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
    load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/ROIs_8_SystPaper.mat'); %loading data
    load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/time_normal.mat'); %loading time
    signtp = cell(length(ROI),length(condition));
    for ii = 1:length(ROI)
        for cc = 1:length(condition)
            load(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/SystPaperROIs/' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.mat'])
            signtp{ii,cc} = table2array(PDn(:,3)); %4D matrix (clusters x time (begin and end) x ROIs x conditions)
        end
    end
    dum2 = permute(dum2,[1 2 4 3]);
    dum3 = dum2(:,:,:,:);
    S.data = dum3;
    S.groups = {'young','elderly'}; %Group labels
    gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
    gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76];
    S.gsubj = gsubj;
    S.time_real = time(1:1026);
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    dumcol = colormap(lines(length(condition) * length(S.groups)* length(ROI))); %extracting some colours from a colormap
    S.colorline = zeros(size(dumcol,1),3); S.colorline(1,:) = dumcol(2,:); S.colorline(2,:) = dumcol(1,:); %swapping colors for young and elderly
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = ROI;
    S.condition_n = condition;
    S.ROIs_labels = ROIN(ROI);
    S.signtp = signtp;
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.qq = qq;
    
    waveplot_groups_local(S) %actual function
    
    if export_l == 1
%         exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/YounVsEld_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.png'],'Resolution',300)
%         exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/YounVsEld_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.eps'],'Resolution',300)
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/YounVsEld_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.pdf'],'Resolution',300)
%             saveas(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/YounVsEld_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.svg'])
%             saveas(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/YounVsEld_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.bmp'])
    end
    

end

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER (ALL SUBJECTS AT THE SAME TIME)

addpath('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots');

for qq = 1:8
    ROI = [qq]; %1 = 'LAC'; 2 = 'LHIT'; 3 = 'LIFG'; 4 = 'MC'; 5 = 'RAC'; 6 = 'RHIT'; 7 = 'RIFG'; 8 = 'VMPFC'; % ROIs in dum2; (UNLESS IT'S APR2020 OR SEQUENCES OF NUMBERS)
    condition = [3]; %1 = Old; 2 = NewT1; 3 = NewT3; (UNLESS IT'S APR2020 OR SEQUENCES OF NUMBERS)
    eld_l = 2; %1 = comparing elderly (60-68 years old) and elderly (>68 years old);
    %2 = high WM elderly, high WM young, low WM elderly, low WM young
    if ROI == 1 || ROI == 5 %LAC anord RAC
        ylimm = [-100 40]; %amplitude limits; leave empty [] for automatic adjustment
    else %the other ROIs
        if eld_l == 1
            ylimm = [-60 30]; %amplitude limits; leave empty [] for automatic adjustment
        else
            ylimm = [-75 35]; %amplitude limits; leave empty [] for automatic adjustment
        end
    end
    export_l = 0; %1 = export images; 0 = not
    
    clear ROIN
    S = [];
    S.conds = {'Old Sys','NewT1 Sys','NewT3 Sys'};
    ROIN{1} = 'LAC'; ROIN{2} = 'LHIT'; ROIN{3} = 'LIFG'; ROIN{4} = 'MC'; ROIN{5} = 'RAC'; ROIN{6} = 'RHIT'; ROIN{7} = 'RIFG'; ROIN{8} = 'VMPFC'; % ROIs in dum2
    load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/ROIs_8_SystPaper.mat'); %loading data
    load('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/time_normal.mat'); %loading time
    %getting significant time-points
    signtp = cell(length(ROI),length(condition));
    for ii = 1:length(ROI)
        for cc = 1:length(condition)
            if eld_l == 1
                load(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/ThreeGroupsAge/' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.mat'])
            else
                load(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/WM_Age/' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.mat'])
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
    S.qq = qq;
    
    waveplot_groups_local(S) %actual function
    
    if eld_l == 1
        barb = 'age3';
    else
        barb = 'WM_age';
    end
    if export_l == 1
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/Ages_WM_' num2str(eld_l) '_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.pdf'],'Resolution',300)
%         exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/Ages_WM_' num2str(eld_l) '_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.png'],'Resolution',300)
%         exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/Ages_WM_' num2str(eld_l) '_' ROIN{ROI(ii)} '_Cond' num2str(condition(cc)) '.eps'],'Resolution',300)
    end
end
   

%%

%% MEG SENSORS

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER (ALL SUBJECTS AT THE SAME TIME)

POS_l = 2; %1 = young > elderly; 2 = elderly > young
condition = 3; %1 = Old; 2 = NewT1; 3 = NewT3
clust_num = 1; %numbers of clusters you want
export_l = 1; %1 = export images; 0 = not


ylimm = [15 140]; %amplitude limits; leave empty [] for automatic adjustment
addpath('/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots');
ROI = 1; %1 = 'LAC'; 2 = 'LHIT'; 3 = 'LIFG'; 4 = 'MC'; 5 = 'RAC'; 6 = 'RHIT'; 7 = 'RIFG'; 8 = 'VMPFC'; % ROIs in dum2;
clear ROIN
S = [];
if condition == 1
    S.conds = {'Old Sys'};
elseif condition == 2
    S.conds = {'NewT1 Sys'};
elseif condition == 3
    S.conds = {'NewT3 Sys'};
end
ROIN{1} = 'MEG sensors';
load(['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/MEG_sensors_significant_clusters/Cond_' num2str(condition) '_clust_' num2str(clust_num) '_POS_l_' num2str(POS_l) '.mat']); %loading data
signtp = {signn};
S.data = data2(:,1:1026,:);
S.groups = {'young','elderly'}; %Group labels
gsubj{1} = [2,5,9,13,14,15,17,19,21,22,24,25,37,39,40,41,44,46,50,53,54,55,57,60,63,64,65,67,69,70,71,72,73,74,75,77,78];
gsubj{2} = [3,4,6,7,8,10,11,12,16,18,20,23,26,27,28,29,30,31,32,33,34,36,38,42,43,45,47,48,49,51,52,56,58,59,61,62,66,68,76];
S.gsubj = gsubj;
S.time_real = time_sel(1:1026);
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 1;
end
dumcol = colormap(lines(length(condition) * length(S.groups)* length(ROI))); %extracting some colours from a colormap
S.colorline = zeros(size(dumcol,1),3); S.colorline(1,:) = dumcol(2,:); S.colorline(2,:) = dumcol(1,:); %swapping colors for young and elderly
S.x_lim = [-0.1 3.4]; % Set x limits
S.y_lim = ylimm; %Set y limits
S.ROI_n = ROI;
S.condition_n = 1;
S.ROIs_labels = ROIN(ROI);
S.signtp = signtp;
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow

waveplot_groups_local(S) %actual function

if export_l == 1
    exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/SEMPRE_TempSeqAges2021/Papers/Paper_Elderly_Bach_MEG/Codes_Data_Plots/Images/MEG_sensors/Cond_' num2str(condition) '_clust_' num2str(clust_num) '_POS_l_' num2str(POS_l) '.pdf'],'Resolution',300)
end



%%





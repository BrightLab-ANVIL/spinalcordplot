function [heatmap,voxelDir,tstats]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,phys_prefix,phys_loc,moco_loc)
% DESCRIPTION
% Use masks created by the x.carpetPlots and produce the carpet plot
% figures in time and frequency domains, plotted with regressors, and
% alongside t-stats.
% 
% USAGE 
% SCheatmap(input_folder,write_out,bySlice,useLevels,TR,phys_prefix,phys_loc,moco_loc)
%
% 
% MANDATORY ARGUMENTS
% input_folder ------->  full path to input folder (x.carpetPlots output folder)
% write_out----------->  1 or 0: 1 will write out files to current folder                        
% bySlice ------------>  1 or 0: 1 will sort the carpetplot by slice, 0 will sort by tissue type
% useLevels ---------->  1 or 0: 1 will indicate vertebral levels on plot, 0 will not. 
%                        WARNING: only use this if CSF is not included in masks
% TR ----------------->  TR in seconds
% phys_prefix -------->  phys regressor prefix  -  e.g. 'sub-03_ses-BH'
% phys_loc ----------->  full path to phys regressors folder
% GLM ---------------->  1 or 0: to use the GLM                             (WIP)
% task --------------->  specify which regressor is the task (full path?)         (WIP)
% moco_loc ----------->  OPTIONAL: full path to 6DOF motion traces (i.e. mc.txt)
%                        enter empty brackets otherwise: []
% 
% 
% OUTPUTS
% heatmap ------------>  heat map matrix (spatially normalized)
% tstats ------------->  matrix of t-statistics  
% voxelDir ----------->  directory of voxels that corresponds to heat map
%                        
%
% Kimberly Hemmerling 2020
% Concepts inspired by Power et al. 2017

%% Do checks and add paths
if nargin ~= 8
    help SCheatmap
    error('Incorrect number of input arguments!')
end
close all
addpath(input_folder)
addpath(phys_loc)
fprintf('\nBeginning... \n \n')
%% Load data
maskdir=dir([input_folder '/mask*ts.txt']);
% Loop through tissue type masks to load data
maskts={};
for i=1:size(maskdir,1)
    % Access each matrix in cell: maskts{i}
    maskts{i}=readmatrix(maskdir(i).name); 
end
% % % GLM=1; task=''; % WIP
%% Delete NaN column 
% The following lines of code delete the NaN column. This column is 
% likely a result of the data import.
for i=1:size(maskdir,1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    [~,nan_col]=find(isnan(maskts{i}(:,:)));
    if nan_col
        maskts{i}(:,nan_col(1))=[];
    end
end
fprintf('NaN column deleted \n \n')
%% Delete first and last slice
firstSlice=min(maskts{1}(3,:));
lastSlice=max(maskts{1}(3,:));
for i=size(maskts,2):-1:1
    current_ts=maskts{i}; 
    for j=size(current_ts,2):-1:1 % Start at end to avoid going > array bounds
        if current_ts(3,j)==firstSlice 
            current_ts(:,j)=[]; %This gets rid of inferior most slice
        elseif current_ts(3,j)==lastSlice
            current_ts(:,j)=[]; %This gets rid of superior most slice
        end
    end
    maskts{i}=current_ts;
end
fprintf('Ignoring first and last slices. \n \n')
%% Rearrange data ordered by tissue type ------------------------ OPTION 1A
if bySlice==0
    fprintf('Sorting by tissue type. \n \n')
    % Rearrange data to one matrix to create a heatmap later
    % Inner SC structures are at the top of the map (i.e. starting with last mask)
    nCol=size(maskts{1},1)-3; % Horizontal axis dim
    nRow=0;
    for i=1:size(maskts,2)
        nRow=nRow+size(maskts{i},2); % Set vertical axis dim
    end
    heatmap=zeros(nRow,nCol);
    voxelDir=zeros(nRow,3); % Directory of voxels to locate in 'heatmap' (x,y,z)
    mean_ts=zeros(nRow,1); % Initialize for spatial normalization
    row=0;
    for i=size(maskts,2):-1:1
        current_ts=maskts{i}; 
        for j=1:size(current_ts,2)
            row=row+1;
            % Starts at row 4 because there is voxel information contained in 1:3
            voxelDir(row,:)=current_ts(1:3,j); % Save voxel coord info to the directory
            heatmap(row,:)=current_ts(4:end,j); % Save timeseries to plot
            mean_ts(row)=mean(current_ts(4:end,j)); % Save means for spatial normalization

        end
    end
end
%% Rearrange data ordered by slice ------------------------------ OPTION 1B
if bySlice==1
    fprintf('Sorting longitudinally by vertebral level/slice. \n \n')
    nCol=size(maskts{1},1); % Set horizontal axis dim
    nRow=0;
    for i=1:size(maskts,2)
        nRow=nRow+size(maskts{i},2); % Set vertical axis dim
    end
    tempPlot=zeros(nRow,nCol);
%     voxelDir=zeros(nRow,3); % Directory of voxels to locate in 'heatmap' (x,y,z)
    mean_ts=zeros(nRow,1); % Initialize for spatial normalization
    num_masks=size(maskts,2);
    a=1; b=0; % a & b are bounds to fill in heatmap
    for i=size(maskts,2):-1:1
        a=b+1;
        b=b+size(maskts{i},2);
        current_ts=maskts{i}; 
        tempPlot(a:b,:)=current_ts';
    end
    
    % Sort by vertebre levels if indicated
    if useLevels==1
        load('vertebralLevels.txt') 
        levels=vertebralLevels';
        tempLevels=[tempPlot(:,1:3) zeros(length(tempPlot),1) tempPlot(:,4:end)];
        for i=1:length(levels)
            x=levels(i,1);
            y=levels(i,2);
            z=levels(i,3);
            level=levels(i,4);
            loc=find( tempPlot(:,1)==x & tempPlot(:,2)==y & tempPlot(:,3)==z );
            if loc>0
                tempLevels(loc,4)=level;
            end
        end
        % Loop through any empty voxels and interpolate if missing (zeros)
        lev_idx=(find(tempLevels(:,4)==0));
        for i=1:length(lev_idx)
            idx=lev_idx(i);
            % Check to avoid exceeding array bounds
            if idx==length(tempLevels)
                tempLevels(idx,4)=tempLevels(idx-1,4);
            else
                tempLevels(idx,4)=round((tempLevels(idx-1,4) + tempLevels(idx+1,4))/2);
            end
        end
        % Using 'descend' puts superior slices at top and inferior at bottom
        heatmap=sortrows(tempLevels,[3 1],'descend');
        % Additionally sort by the vertebral levels
        heatmap=sortrows(heatmap,4);
        % Move voxel info to directory and delete from plot
        voxelDir=heatmap(:,1:4);
        heatmap(:,1:4)=[]; 
    else
        % This section is without vertebral level information
        % Sorting primarily by slice (3), secondarily by x (1)
        heatmap=sortrows(tempPlot,[3 1],'descend');
        % Move voxel info to directory and delete from plot
        voxelDir=heatmap(:,1:3);
        heatmap(:,1:3)=[]; 
    end    
    % Calculate mean of each voxel
    for i=1:size(heatmap,1)
        mean_ts(i)=mean(heatmap(i,:));
    end
end
%% Demean and normalize heatmap
heatmap_preNorm=heatmap;
heatmap=zeros(size(heatmap_preNorm));
for r=1:nRow
    heatmap(r,:)=(heatmap_preNorm(r,:)-mean_ts(r))./range(mean_ts);
end
%% Load physiological and motion data
% Load convolved physiological regressors
physHR_conv=load(strcat(phys_prefix,'_HR_CRFconv.txt'));
physCO2_conv=load(strcat(phys_prefix,'_CO2_HRFconv.txt'));
physRVT_conv=load(strcat(phys_prefix,'_RVT_RRFconv.txt'));
physO2_conv=load(strcat(phys_prefix,'_O2_HRFconv.txt'));
% Load nonconvolved physiological traces
physHR=load(strcat(phys_prefix,'_HR.txt'));
physCO2=load(strcat(phys_prefix,'_CO2.txt'));
physRVT=load(strcat(phys_prefix,'_RVT.txt'));
physO2=load(strcat(phys_prefix,'_O2.txt'));
% Load motion and demean
if moco_loc
    motion=load(moco_loc);
    for i=1:size(motion,2)
        motion(:,i)=motion(:,i)-mean(motion(:,i));
    end
    if size(motion,2) ~=6
        moco_loc=[];
        warning('Motion file does not have 6 columns. Will ignore.')
    end
end
% Demean convolved physiological regressors
physHR_conv=physHR_conv-mean(physHR_conv);
physCO2_conv=physCO2_conv-mean(physCO2_conv);
physRVT_conv=physRVT_conv-mean(physRVT_conv);
physO2_conv=physO2_conv-mean(physO2_conv);
% Define caxis bounds for heatmap
c1=-0.4; c2=-c1;
% Do more checks (equal lengths)
if ~isequal(size(heatmap,2),length(physHR), length(physCO2), length(physRVT), length(physO2), size(motion,1))
    error('Length of traces does not match number of TRs (%d). Check fMRI data, physiological traces, and motion traces.', size(heatmap,2))
end
%% Plot data ordered by tissue
if bySlice==0
    %%%%%% Set phys data to plot:
    phys=[physCO2 physHR];
    %%%%%%
    gmEnds=size(maskts{size(maskts,2)},2); % Size of GM graph portion
    figure('Name','By Tissue','Renderer', 'painters', 'Position', [50 1000 630 700])
    subplot(4,1,[3,4])
    imagesc(heatmap)
    set(gca,'YTickLabel',[]); pbaspect([2 1 1])
    xlabel('{\bfTRs}')
    colormap gray
    caxis([c1 c2])
    % Draw white line b/t tissue types
    line([0 nRow], [gmEnds+0.5 gmEnds+0.5],'Color','white','LineWidth',0.7) 

    % Add physio
    subplot(411)
    plot(phys(:,1),'c','LineWidth',1.5); xlim([0 length(phys)])
    ylabel({'{\bfP_{ET}CO_{2}}','[mmHg]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    subplot(412)
    plot(phys(:,2),'g','LineWidth',1.5); xlim([0 length(phys)])
    ylabel({'{\bfHR}','[bpm]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')

    % Semi-manually define FSL's greengray colormap (using polynomials)
    x=1:256; p1 = -7.3246e-20; p2 = 7.6618e-17; p3 = -3.2357e-14; p4 = 7.0892e-12; p5 = -8.6807e-10;
    p6 = 5.9635e-08; p7 = -1.9699e-06; p8 = 1.1824e-05; p9 = 0.0010542; p10 = -0.0020195;
    d1 = p1.*x.^9 + p2.*x.^8 + p3.*x.^7 + p4.*x.^6 + p5.*x.^5 + p6.*x.^4 + p7.*x.^3 + p8.*x.^2 + p9.*x + p10 ;
    d1(1)=0; d1(256)=1;
    p1 = -4.1809e-14; p2 = 3.5582e-11; p3 = -1.1274e-08; p4 = 1.6419e-06;
    p5 = -0.00011092; p6 = 0.0069924; p7 = -0.0011795;
    d2 = p1.*x.^6 + p2.*x.^5 + p3.*x.^4 + p4.*x.^3 + p5.*x.^2 + p6.*x + p7 ;
    p1 = -4.3403e-20; p2 = 4.5503e-17; p3 = -1.9218e-14; p4 = 4.2255e-12; p5 = -5.3269e-10;
    p6 = 4.1003e-08; p7 = -1.8802e-06; p8 = 4.8024e-05; p9 = -0.00010991; p10 = 0.0027377;
    d3 = p1.*x.^9 + p2.*x.^8 + p3.*x.^7 + p4.*x.^6 + p5.*x.^5 + p6.*x.^4 + p7.*x.^3 + p8.*x.^2 + p9.*x + p10 ;
    d3(256)=1;
    greengrayMap=[d1' d2' d3'];
%     Load FSL's colormap:
%     greengrayMap=load('/usr/local/fsl/fslpython/envs/fslpython/lib/python3.7/site-packages/fsleyes/assets/colourmaps/brain_colours/greengray.cmap');
    tissueTypes=ones(size(maskts,2),3);
    tissueTypes(:,1)=size(maskts,2):-1:1;
    idx=1;
    tissueTypes(idx,3)=gmEnds;
    for t=size(maskts,2)-1:-1:1
        idx=idx+1;
        tissueTypes(idx,3)=size(maskts{t},2)+tissueTypes(idx-1,3);
    end
    tissueTypes(2:end,2)=tissueTypes(1:end-1,3)+1;
    tissueColorbar=zeros(size(heatmap,1),1);
    for t=1:size(tissueTypes,1)
        tissueColorbar(tissueTypes(t,2):tissueTypes(t,3))=tissueTypes(t,1);
    end
    % Plot tissue next to heatmap [x0 y0 width height]
    subplot('Position',[0.11 0.125 0.019 0.348])
    imagesc(tissueColorbar); colormap(gca,greengrayMap)
    caxis([0 max(tissueColorbar(:,1))]);
    set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
end

%% Plot data ordered by slice (and vertebral level if opted for)
if bySlice==1
    %%%%%% Set phys data to plot:
    phys=[physCO2 physHR];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','By Slice','Renderer', 'painters', 'Position', [50 1000 630 700])
    subplot(4,1,[3,4])
    imagesc(heatmap)
    set(gca,'YTickLabel',[]); pbaspect([2 1 1])
    xlabel('{\bfTRs}'); set(gca,'FontSize',12)
%     ylabel('\leftarrow Inferior                     Superior \rightarrow')
    colormap gray
    caxis([c1 c2])
    % Add physio
    subplot(411)
    plot(phys(:,1),'c','LineWidth',1.5); xlim([0 length(phys)]); set(gca,'XTickLabel',[],'FontSize',12)
    ylabel({'{\bfP_{ET}CO_{2}}','[mmHg]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    subplot(412)
    plot(phys(:,2),'g','LineWidth',1.5); xlim([0 length(phys)]); set(gca,'XTickLabel',[],'FontSize',12)
    ylabel({'{\bfHR}','[bpm]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    if useLevels==1
        % Add small indicators of where the vertebral levels are
        levChange=(min(voxelDir(:,4)):max(voxelDir(:,4))-1); l=1;
        levChange=[levChange' zeros(1,range(voxelDir(:,4)))'];
        for i=2:length(voxelDir)
            if voxelDir(i-1,4)~=voxelDir(i,4)
                levChange(l,2)=i;
                l=l+1;
            end
        end
        % The levChange vector has the row locations of the first row of
        % the new vertebral level
        subplot(4,1,[3,4])
        for level=levChange(:,2)
            line([nCol-10 nCol], [level-0.5 level-0.5], 'Color','y','LineWidth',1)
        end
        % Create a colorbar of the vertebral levels and add to figure
        vertebralLevels=ones(size(levChange,1)+1,3);
        vertebralLevels(1:end,1)=[levChange(:,1); levChange(end,1)+1];
        vertebralLevels(1:end,3)=[levChange(:,2); size(heatmap,1)];
        vertebralLevels(2:end,2)=vertebralLevels(1:end-1,3)+1;
        blueLightBlueMap = [zeros(256,1), linspace(0,1,256)', ones(256,1)]; % from FSL
        vertebralColorbar=zeros(size(heatmap,1),1);
        for v=1:size(vertebralLevels,1)
            vertebralColorbar(vertebralLevels(v,2):vertebralLevels(v,3))=vertebralLevels(v,1);
        end
        % Plot vertebral level next to heatmap
        subplot('Position',[0.11 0.125 0.019 0.348])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 max(vertebralLevels(:,1))]); 
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    end
end

%% PHASE / FREQUENCY / POWER
%%%%%% Set phys data to plot:
phys=[physCO2 physHR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform & calculations for phys data
for ph=1:size(phys,2)
    fs=1/TR; % 0.5 for TR=2
    phys_FT(:,ph)=fft(phys(:,ph)); % Fourier transformed phys regressor
    n=length(phys_FT(:,ph));
    phys_power(:,ph)=abs(phys_FT(:,ph)).^2/n;                % non-zero centered power
    phys_powershift(:,ph)=abs(fftshift(phys_FT(:,ph))).^2/n; % zero-centered power
end
freq = (0:n-1)*(fs/n);      % non-zero centered frequency range
fshift=(-n/2:n/2-1)*(fs/n); % zero-centered frequency range

% Fourier transform & calculations for fMRI timeseries
heatmap_FT=zeros(size(heatmap));
powerPlot=zeros(size(heatmap));
powershiftPlot=zeros(size(heatmap));
for i=1:size(heatmap,1)
    voxel_ts_FT=fft(heatmap(i,:));
    power=abs(voxel_ts_FT).^2/n;                % non-zero centered power
    powershift=abs(fftshift(voxel_ts_FT)).^2/n; % zero-centered power
    heatmap_FT(i,:)=voxel_ts_FT; % Fourier transformed (1D) carpet plot
    powerPlot(i,:)=power;
    powershiftPlot(i,:)=powershift;
end
% REMOVING ZERO FREQ., (for ease of visualization) %
powershiftPlot(:,136)=0; phys_powershift(136,:)=0; 
powerPlot(:,1)=0; phys_power(1,:)=0;
% Plot carpetplot frequency
figure('Name','Power','Renderer', 'painters', 'Position', [50 1000 630 700])
subplot(4,1,[3,4])
imagesc(powerPlot(:,1:round(n/2))); set(gca,'xtick',[]); colormap pink
caxis([0 0.8]); 
set(gca,'YTickLabel',[]); pbaspect([2 1 1])
% % if bySlice==1
% %     ylabel('\leftarrow Inferior                 Superior \rightarrow')
% % elseif bySlice==0
% %     ylabel('\leftarrow Outer/WM                 Inner/GM \rightarrow')
% % end
if bySlice==1
    % Add small indicators of where the vertebral levels are
    if useLevels==1
        levChange=(min(voxelDir(:,4)):max(voxelDir(:,4))-1); l=1;
        levChange=[levChange' zeros(1,range(voxelDir(:,4)))'];
        for i=2:length(voxelDir)
            if voxelDir(i-1,4)~=voxelDir(i,4)
                levChange(l,2)=i;
                l=l+1;
            end
        end
        subplot(4,1,[3,4])
        for level=levChange(:,2)
            line([nCol/2-5 nCol/2], [level-0.5 level-0.5], 'Color','white')
        end
    end
    % Plot vertebral level next to heatmap
    subplot('Position',[0.11 0.125 0.019 0.348])
    imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
    caxis([0 max(vertebralLevels(:,1))]); 
    set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
elseif bySlice==0
    % Plot tissue next to heatmap [x0 y0 width height]
    subplot('Position',[0.11 0.125 0.019 0.348])
    imagesc(tissueColorbar); colormap(gca,greengrayMap)
    caxis([0 max(tissueColorbar(:,1))]);
    set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
end
% Add physio to plot
subplot(411)
plot(freq,phys_power(:,1),'c','LineWidth',1.5)
title({'{\bfP_{ET}CO_{2}}'})
xlim([0 0.25]) % limit plot by Nyquist frequency
subplot(412)
plot(freq,phys_power(:,2),'g','LineWidth',1.5)
title({'{\bfHR}'})
xlim([0 0.25])

%% Run GLM
if moco_loc
    X=[ones(size(heatmap,2),1) physCO2_conv physHR_conv motion];
    % Demean the design matrix
    for i=2:size(X,2)
        X(:,i)=(X(:,i)-mean(X(:,i)))./range(X(:,i));
    end
    % Loop through each voxel Y individually when calculating beta=pinv(X)*Y
    % the length of beta will be the number of regressors. Then calculate t-stat
    Bfit = zeros(nRow,9);%
    tstats=zeros(nRow,8);
    pin_mult=pinv(X)*pinv(X)';
    pinv_X=pinv(X);
    % Define contrast vectors
    contrast=[0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0; 0 0 0 0 1 0 0 0 0;...
    0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1]; 
    for row=1:nRow
        % Calculate beta coefficients
        Y=heatmap(row,:)';
        beta=pinv_X*Y;
        Bfit(row,:)=beta;
        % Calculate t-statistics
        r = Y-X*beta;
        sumsq=r'*r; 
        DOF=size(X,1)-size(X,2); % Number of timepoints minus regressors
        rvar = sumsq/DOF;
        for i=1:size(contrast,1)
            c=contrast(i,:);
            tstats(row,i)=(c*beta)/sqrt(rvar*c*(pin_mult)*c');
        end
    end
else
    tstats=[];
end
%% Motion, phys, GLM plot
if moco_loc
    % Define colormaps
    greenMap = [zeros(256,1), linspace(0,1,256)', zeros(256,1)];
    cyanMap = [zeros(256,1), linspace(0,1,256)', linspace(0,1,256)'];
    magMap = [linspace(0,1,256)', zeros(256,1), linspace(0,1,256)'];
    rotationsMap = [linspace(0,1,256)', linspace(0,0.5686,256)', zeros(256,1)];
    translationsMap = [linspace(0,1,256)', zeros(256,1), zeros(256,1)];
    figure('Name','GLM and Regressors','Renderer', 'painters','Position',[50 1 693 804])%[50 1000 830 944])
    % Physio ( Position: [x0 y0 width height] )
    subplot('Position',[0.13 0.8472 0.3347 0.0604])
    plot(physCO2,'c','LineWidth',1.5); xlim([0 length(physCO2)]); set(gca,'xtick',[],'FontSize',12)
    ylabel({'{\bfP_{ET}CO_{2}}','[mmHg]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    subplot(10,2,3)
    subplot('Position',[0.13 0.7808 0.3347 0.0604])
    plot(physHR,'g','LineWidth',1.5); xlim([0 length(physHR)]); set(gca,'xtick',[],'FontSize',12)
    ylabel({'{\bfHR}','[bpm]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    % Carpet plot
    two_sd=2*std2(heatmap);
    subplot('Position', [0.13 0.6131 0.3347 0.1442])
    imagesc(heatmap)
    set(gca,'YTickLabel',[],'FontSize',11); pbaspect([2 1 1])
    caxis([c1 c2])
    xlabel('{\bfTRs}')
    colormap gray
    % Add small indicators of where the vertebral levels are
    if bySlice==1
        if useLevels==1
            levChange=(min(voxelDir(:,4)):max(voxelDir(:,4))-1); l=1;
            levChange=[levChange' zeros(1,range(voxelDir(:,4)))'];
            for i=2:length(voxelDir)
                if voxelDir(i-1,4)~=voxelDir(i,4)
                    levChange(l,2)=i;
                    l=l+1;
                end
            end
            for level=levChange(:,2)
                line([nCol-10 nCol], [level-0.5 level-0.5], 'Color','white')
            end
        end
        % Plot vertebral level next to heatmap
        subplot('Position',[0.12 0.6131 0.01 0.1442])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 max(vertebralLevels(:,1))]); 
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    elseif bySlice==0
        % Plot tissue next to heatmap [x0 y0 width height]
        subplot('Position',[0.12 0.6131 0.01 0.1442])
        imagesc(tissueColorbar); colormap(gca,greengrayMap)
        caxis([0 max(tissueColorbar(:,1))]);
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    end    
    % Motion regressors
    subplot('Position',[0.13 0.5292-0.01 0.3347 0.0525]) % Rx 
    plot(motion(:,1),'Color',[1 0.5686 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
    ylabel('R_x','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    subplot('Position',[0.13 0.4707-0.01 0.3347 0.0525]) % Ry 
    plot(motion(:,2),'Color',[1 0.5686 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
    ylabel('R_y','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    subplot('Position',[0.13 0.4122-0.01 0.3347 0.0525]) % Rz 
    plot(motion(:,3),'Color',[1 0.5686 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
    ylabel('R_z','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    subplot('Position',[0.13 0.3537-0.01 0.3347 0.0525]) % Tx 
    plot(motion(:,4),'Color',[1 0 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
    ylabel('T_x','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    subplot('Position',[0.13 0.2952-0.01 0.3347 0.0525]) % Ty 
    plot(motion(:,5),'Color',[1 0 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
    ylabel('T_y','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    subplot('Position',[0.13 0.2367-0.01 0.3347 0.0525]) % Tz 
    plot(motion(:,6),'Color',[1 0 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
    ylabel('T_z','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    % GLM tstats ( Position: [x0 y0 width height] )
    subplot('Position',[0.4680 0.6131 0.0419 0.1442]) % CO2
    imagesc(abs(tstats(:,1))); set(gca,'xtick',[],'ytick',[]); colormap(gca,cyanMap); caxis([0 5])
    subplot('Position',[0.5099 0.6131 0.0419 0.1442]) % HR
    imagesc(abs(tstats(:,2))); set(gca,'xtick',[],'ytick',[]); colormap(gca,greenMap); caxis([0 5])
    subplot('Position',[0.5518 0.6131 0.1257 0.1442]) % MOTION x6
    imagesc(abs(tstats(:,3:5))); set(gca,'xtick',[],'ytick',[]);  colormap(gca,rotationsMap); caxis([0 5])
    title('t-statistics','HorizontalAlignment','right','FontSize',12)
    subplot('Position',[0.6775 0.6131 0.1257 0.1442]) % MOTION x6
    imagesc(abs(tstats(:,6:8))); set(gca,'xtick',[],'ytick',[],'FontSize',12);  colormap(gca,translationsMap); caxis([0 5])
end
%% Reorganize data and plot by t-statistic magnitude (CO2)
if moco_loc
    % CO2 heatmap
    temp_sorter=[abs(tstats(:,1)) heatmap];
    heatmap_CO2_sort=sortrows(temp_sorter,1,'descend');
    heatmap_CO2_sort=heatmap_CO2_sort(:,2:end);
    % CO2 tstats
    temp_sorter=[abs(tstats(:,1)) tstats];
    tstats_CO2_sort=sortrows(temp_sorter,1,'descend');
    tstats_CO2_sort=tstats_CO2_sort(:,2:end);
    % Physio
    figure('Name','Plot data by HR t-statistic magnitude','Renderer', 'painters', 'Position', [50 1000 1398 621])%1457 648])
    subplot(4,2,1)
    plot(physCO2,'c','LineWidth',1.5); xlim([0 length(physCO2)]); set(gca,'xtick',[],'FontSize',12)
    ylabel({'{\bfP_{ET}CO_{2}}','[mmHg]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    subplot(4,2,3)
    plot(physHR,'g','LineWidth',1.5); xlim([0 length(physHR)]); set(gca,'xtick',[],'FontSize',12)
    ylabel({'{\bfHR}','[bpm]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    % Heatmap
    subplot(4,2,[5,7])
    imagesc(heatmap_CO2_sort)
    set(gca,'YTickLabel',[],'FontSize',12); pbaspect([2 1 1])
    caxis([c1 c2])
    xlabel('{\bfTRs}')
    colormap gray
    % GLM tstats ( Position: [x0 y0 width height] )
    subplot('Position',[0.4680 0.1100 0.0419 0.3768]) % CO2
    imagesc(abs(tstats_CO2_sort(:,1))); set(gca,'xtick',[],'ytick',[]); colormap(gca,cyanMap); caxis([0 5]); title('CO2')
    subplot('Position',[0.5099 0.1100 0.0419 0.3768]) % HR
    imagesc(abs(tstats_CO2_sort(:,2))); set(gca,'xtick',[],'ytick',[]); colormap(gca,greenMap); caxis([0 5]); title('HR')
end

%% Reorganize data and plot by t-statistic magnitude (HR)
if moco_loc
    % HR heatmap
    temp_sorter=[abs(tstats(:,2)) heatmap];
    heatmap_HR_sort=sortrows(temp_sorter,1,'descend');
    heatmap_HR_sort=heatmap_HR_sort(:,2:end);
    % HR tstats
    temp_sorter=[abs(tstats(:,2)) tstats];
    tstats_HR_sort=sortrows(temp_sorter,1,'descend');
    tstats_HR_sort=tstats_HR_sort(:,2:end);
    clear temp_sorter
    % Physio
    figure('Name','Plot data by HR t-statistic magnitude','Renderer', 'painters', 'Position', [50 1000 1398 621])%1457 648])
    subplot(4,2,1)
    plot(physCO2,'c','LineWidth',1.5); xlim([0 length(physCO2)]); set(gca,'xtick',[],'FontSize',12)
    ylabel({'{\bfP_{ET}CO_{2}}','[mmHg]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    subplot(4,2,3)
    plot(physHR,'g','LineWidth',1.5); xlim([0 length(physHR)]); set(gca,'xtick',[],'FontSize',12)
    ylabel({'{\bfHR}','[bpm]'},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    % Heatmap
    subplot(4,2,[5,7])
    imagesc(heatmap_HR_sort)
    set(gca,'YTickLabel',[],'FontSize',12); pbaspect([2 1 1])
    caxis([c1 c2])
    xlabel('{\bfTRs}')
    colormap gray
    % GLM tstats ( Position: [x0 y0 width height] )
    subplot('Position',[0.4680 0.1100 0.0419 0.3768]) % CO2
    imagesc(abs(tstats_HR_sort(:,1))); set(gca,'xtick',[],'ytick',[]); colormap(gca,cyanMap); caxis([0 5]); title('CO2')
    subplot('Position',[0.5099 0.1100 0.0419 0.3768]) % HR
    imagesc(abs(tstats_HR_sort(:,2))); set(gca,'xtick',[],'ytick',[]); colormap(gca,greenMap); caxis([0 5]); title('HR')
end

%% Condensed motion regressor plot
if moco_loc
    figure('Name','Translations and Rotations','Renderer', 'painters', 'Position', [50 1000 683 700])
    subplot(4,1,[1,2])
    imagesc(heatmap)
%     if bySlice==1
%         ylabel('\leftarrow Inferior                    Superior \rightarrow')
%     elseif bySlice==0
%         
%         ylabel('\leftarrow Outer/WM                           Inner/GM \rightarrow')
%     end
    set(gca,'YTickLabel',[],'XTickLabel',[]); pbaspect([2 1 1])
    caxis([c1 c2]); colormap gray
    subplot(413)
    hold on
    plot(motion(:,1),'Color',[1 0.5686 0],'LineWidth',1)
    plot(motion(:,2),'Color',[1 0.5686 0],'LineWidth',1)
    plot(motion(:,3),'Color',[1 0.5686 0],'LineWidth',1)
    ylim([-0.06 0.06]); xlim([1 length(motion)]); set(gca,'FontSize',12,'XTickLabel',[])
    ylabel({'{\bfRotations}','[rads]'},'rotation',90)
    hold off
    subplot(414)
    hold on
    plot(motion(:,4),'Color',[1 0 0],'LineWidth',1.5)
    plot(motion(:,5),'Color',[1 0 0],'LineWidth',1.5)
    plot(motion(:,6),'Color',[1 0 0],'LineWidth',1.5)
    ylim([-6 6]); xlim([1 length(motion)])
    ylabel({'{\bfTranslations}','[mm]'},'rotation',90)
    xlabel('{\bfTRs}'); set(gca,'FontSize',12)
    hold off
    if bySlice==1
        % Plot vertebral level next to heatmap
        subplot('Position',[0.111 0.5482 0.019 0.3768])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 max(vertebralLevels(:,1))]); 
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    elseif bySlice==0
        if useLevels==1
            % Draw line b/t tissue types
            subplot(4,1,[1,2])
            line([0 nRow], [gmEnds+0.5 gmEnds+0.5],'Color','white','LineWidth',0.7) 
        end
        % Plot tissue next to heatmap [x0 y0 width height]
        subplot('Position',[0.111 0.5482 0.019 0.3768])
        imagesc(tissueColorbar); colormap(gca,greengrayMap)
        caxis([0 max(tissueColorbar(:,1))]);
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    end
end

%% DVARS plotting (and FD?)
dvars=load('dvars.txt');
dvars(1)=NaN;
figure('Name','DVARS', 'Renderer', 'painters','Position', [50 1000 630 500])
subplot(311)
plot(1:n,dvars,'b','LineWidth',1.5); xlim([0 length(dvars)]); ylim([0 max(dvars)])
ylabel('DVARS','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
subplot(3,1,[2,3])
imagesc(heatmap)
set(gca,'YTickLabel',[]); pbaspect([2 1 1])
xlabel('{\bfTRs}'); caxis([c1 c2]); colormap gray
if bySlice==1
    % Plot vertebral level next to heatmap
    subplot('Position',[0.11 0.124 0.019 0.488])
    imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
    caxis([0 max(vertebralLevels(:,1))]); 
    set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
elseif bySlice==0
    % Plot tissue next to heatmap [x0 y0 width height]
    subplot('Position',[0.11 0.124 0.019 0.488])
    imagesc(tissueColorbar); colormap(gca,greengrayMap)
    caxis([0 max(tissueColorbar(:,1))]);
    set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
end
%% Write out files if requested
if write_out==1
    if bySlice==1
        order='bySlice';
    elseif bySlice==0
        order='byTissue';
    end
    writematrix(heatmap,[phys_prefix order '_heatmap.txt'],'Delimiter','tab')
    writematrix(voxelDir,[phys_prefix order 'voxelDirectory.txt'],'Delimiter','tab')
    if tstats
        writematrix(voxelDir,[phys_prefix order '_tstats.txt'],'Delimiter','tab')
    end
end


fprintf('\n...done!\n')
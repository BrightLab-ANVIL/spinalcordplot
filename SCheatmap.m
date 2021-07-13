function [heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc,options)
% SCheatmap DESCRIPTION
% Use masks created by the x.carpetPlots and produce the carpet plot
% figures in time and frequency domains, plotted with regressors, and
% alongside t-stats.
% 
% USAGE (minimum arguments)
% SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc)
%
% 
% MANDATORY ARGUMENTS
% input_folder ------->  full path to input folder (x.heatmapPrep output folder)
% write_out----------->  1 or 0: 1 will write out files to current folder                        
% bySlice ------------>  1 or 0: 1 will sort the heatmap by slice, 0 will sort by tissue type
% useLevels ---------->  1 or 0: 1 will indicate vertebral levels on plot, 0 will not. 
%                        WARNING: only use this if CSF is not included in masks
% TR ----------------->  TR in seconds
% prefix ------------->  regressor/trace file prefix  -  e.g. 'sub-03'
% trace_loc ---------->  full path to folder with physiological, task or other traces
%                        required format: 'prefix_Trace' - e.g. sub-03_CO2.txt'
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS
% "moco", "/path/file.txt" OR ["/path/mocoX.txt" "/path/mocoY.txt"] --> path or paths to motion data
%                                                 Either: path to one motion file (6DOF or X&Y) or 
%                                                 two paths to X and Y slicewise motion
% "mocoLabel", ["Tx" "Ty" "Tz" "Rx" "Ry" "Rz"] -> order of columns
%
% "cBound", 0.3 --------------------------------> caxis abs value of limit around zero (default: 0.4)
% "PlotSmoothData", 1 --------------------------> 1 or 0: 1 will plot smoothed data (default: 0)
% "Traces", ["RVT" "HR"] -----------------------> choose traces to focus on (same length as num TRs)
%                                                 1x2 string array (default: ["CO2" "HR"])
% "GLMtask", ["Task" "HR"] ---------------------> indicate that you would like to perform a GLM and 
%                                                 identify two traces to be used as GLM regressors
% "plots", 1 -----------------------------------> 0: output no plots; 1: output basic plot only; 
%                                                 2: output all available plots (default)
% "stim", "/path/file.txt" ---------------------> Add binary stimulus timings to plot of first trace 
% "slices", [1 10] -----------------------------> Limit slices to include in output heatmap (inclusive)
%                                                 default excludes top & bottom slices
%            
%
% Example using name-value pair arguments (assuming minimum arguments as ...)
% SCheatmap(..., "mocoLabel", ["Tx" "Ty" "Tz" "Rx" "Ry" "Rz"], "PlotSmoothData", 1, "GLMtask", ["CO2" "HR"])
% 
% OUTPUTS
% heatmap ------------>  heat map matrix (spatially normalized)
% tstats ------------->  matrix of t-statistics  
% voxelDir ----------->  directory of voxels that corresponds to heat map
%                        
%
% NOTE: This script assumes a certain format for the prefix, followed by
% the trace name, as follows: 'prefix_Trace.txt'
%   E.g., Normal trace: 'sub-03_ses-BH_HR.txt'
%         Convolved trace: 'sub-03_ses-BH_HR_CRFconv.txt'
%         (similarly, _CO2_HRFconv.txt, _RVT_RRFconv.txt, _O2_HRFconv.txt)
% 
% 
% Kimberly Hemmerling 2020
% Concepts inspired by Power et al. 2017

%% Do checks and add paths
arguments
    input_folder (1,:) char
    write_out (1,1) {mustBeMember(write_out,[0,1])}
    bySlice (1,1) {mustBeMember(bySlice,[0,1])}
    useLevels (1,1) {mustBeMember(useLevels,[0,1])}
    TR (1,1) {mustBeNumeric,mustBePositive}
    prefix (1,:) char
    trace_loc (1,:) char
    options.moco (1,:) string = "-"
    options.mocoLabel (1,:) string = ["Rx" "Ry" "Rz" "Tx" "Ty" "Tz"]
    options.cBound (1,1) double {mustBePositive} = 0.4
    options.PlotSmoothData (1,1) {mustBeMember(options.PlotSmoothData,[0,1])} = 0
    options.Traces (1,2) string = ["CO2" "HR"] %{mustBeMember(options.Traces,["CO2","HR","O2","RVT"])}
    options.GLMtask (1,2) string = ["-" "-"] %{mustBeMember(options.GLMtask,["CO2","HR","O2","RVT","-"])} = ["-" "-"]
    options.plots (1,1) {mustBeMember(options.plots,[0,1,2])} = 2
    options.stim (1,1) string = "-"
    options.slices (1,2) double {mustBeNumeric,mustBeInteger} = [-1 -1]
end
close all
addpath(input_folder)
addpath(trace_loc)
fprintf('\nBeginning... \n \n')
%% Load data
if options.PlotSmoothData==0
    maskdir=dir([input_folder '/mask*ts.txt']);
elseif options.PlotSmoothData==1
    fprintf('\nUsing data smoothed within masks.\n')
    maskdir=dir([input_folder '/blur_mask*ts.txt']);
else
    error('PlotSmoothData input argument not binary.')
end
% Loop through tissue type masks to load data
maskts={};
for i=1:size(maskdir,1)
    % Access each matrix in cell: maskts{i}
    maskts{i}=readmatrix(maskdir(i).name); 
end
% % GLM=1; task=''; % WIP
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
fprintf('\nNaN column deleted \n \n')
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
        % NEW: to make within each tissue mask sorted from top to bottom
        temp_maskts=maskts{i}';
        temp_maskts=sortrows(temp_maskts,[3 1],'descend');
        maskts{i}=temp_maskts';
        % % % % % % % % % % % % % % % % % % % % 
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
%% Limit slices to include in heatmap ----- new, 04/14/21
if useLevels==1
    % Set for vertebral level colorbar prior to deleting rows to correctly
    % correspond with output maps of vertebral level masks
    maxLevel=max(voxelDir(:,4));
end
if all(options.slices >= 0)
    bottomSlice=options.slices(1);
    topSlice=options.slices(2);
    idxs=[];
    for v=1:size(voxelDir,1)
        if (voxelDir(v,3)>topSlice) || (voxelDir(v,3)<bottomSlice)
            idxs=[idxs v];
        end
        
    end
    voxelDir(idxs,:)=[];
    heatmap(idxs,:)=[];
    nRow=size(voxelDir,1);
    % Recalculate mean timeseries of each voxel after deletions
    mean_ts=zeros(nRow,1);
    for i=1:size(heatmap,1)
        mean_ts(i)=mean(heatmap(i,:));
    end
% else
%     bottomSlice=max(voxelDir(:,3))-1; % Exclude bottom-most slice
%     topSlice=min(voxelDir(:,3))+1; % Exclude top-most slice
end
%% Demean and standardize/normalize heatmap (mean normalization)
heatmap_preNorm=heatmap;
heatmap=zeros(size(heatmap_preNorm));
for r=1:nRow
%     heatmap(r,:)=(heatmap_preNorm(r,:)-mean_ts(r))./range(mean_ts);
% Realized that I may have been dividing by the range of the means of the
% ts... which doesn't make sense. Changing to dividing by range of each ts
    heatmap(r,:)=(heatmap_preNorm(r,:)-mean_ts(r))./range(heatmap_preNorm(r,:));
end
%% Define colormaps
% Tissue colormap:
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
% Load FSL's colormap:
%greengrayMap=load('/usr/local/fsl/fslpython/envs/fslpython/lib/python3.7/site-packages/fsleyes/assets/colourmaps/brain_colours/greengray.cmap');
% Slice colormap:
blueLightBlueMap = [zeros(256,1), linspace(0,1,256)', ones(256,1)]; % from FSL
%% Define caxis bounds for heatmap (user input or default 0.4)
c1=-options.cBound; c2=options.cBound;
%% Plot basic plot then exit function
if (options.plots==1) && (bySlice==0)
    figure('Name','Basic Plot: By Tissue Type','Renderer', 'painters', 'Position', [50 1000 887 538])
    imagesc(heatmap)
    set(gca,'YTickLabel',[]); set(gca,'FontSize',20); pbaspect([2 1 1])
    xlabel('{\bfTRs}')
    colormap gray
    caxis([c1 c2])
    tissueTypes=ones(size(maskts,2),3);
    tissueTypes(:,1)=size(maskts,2):-1:1;
    idx=1;
    gmEnds=size(maskts{size(maskts,2)},2);
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
    if all(options.slices >= 0)
        tissueColorbar(idxs,:)=[];
    end
    % Plot tissue next to heatmap [x0 y0 width height]
    subplot('Position',[0.11 0.197 0.019 0.64])
    imagesc(tissueColorbar); colormap(gca,greengrayMap)
    caxis([0 max(tissueColorbar(:,1))]);
    set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[]) 
    freqmap=0;
    fprintf('\nPlotted basic plot... done!\n')
    return
elseif (options.plots==1) && (bySlice==1)
    figure('Name','Basic Plot: By Vertebral Level','Renderer', 'painters', 'Position', [50 1000 887 538])
    imagesc(heatmap)
    set(gca,'YTickLabel',[]); pbaspect([2 1 1])
    xlabel('{\bfTRs}'); set(gca,'FontSize',20)
    colormap gray
    caxis([c1 c2])
    if useLevels==1
        levChange=(min(voxelDir(:,4)):max(voxelDir(:,4))-1); l=1;
        levChange=[levChange' zeros(1,range(voxelDir(:,4)))'];
        for i=2:length(voxelDir)
            if voxelDir(i-1,4)~=voxelDir(i,4)
                levChange(l,2)=i;
                l=l+1;
            end
        end
        % Create a colorbar of the vertebral levels and add to figure
        vertebralLevels=ones(size(levChange,1)+1,3);
        vertebralLevels(1:end,1)=[levChange(:,1); levChange(end,1)+1];
        vertebralLevels(1:end,3)=[levChange(:,2); size(heatmap,1)];
        vertebralLevels(2:end,2)=vertebralLevels(1:end-1,3)+1;
        vertebralColorbar=zeros(size(heatmap,1),1);
        for v=1:size(vertebralLevels,1)
        vertebralColorbar(vertebralLevels(v,2):vertebralLevels(v,3))=vertebralLevels(v,1);
        end
        % Plot vertebral level next to heatmap
        subplot('Position',[0.11 0.197 0.019 0.64])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 maxLevel]); maxLevel
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
        freqmap=0;
    end
    fprintf('\nPlotted basic plot... done!\n')
    return
end
%% Load physiological, motion, or trace data
% Load convolved physiological regressors if user chooses GLMtask
if all(options.GLMtask ~= "-")
    suffix.(options.GLMtask(1))='';
    suffix.(options.GLMtask(2))='';
end
% To do: increase robustness here for no available convolved regressors
suffix.HR='_CRFconv';
suffix.CO2='_HRFconv';
suffix.RVT='_RRFconv';
suffix.O2='_HRFconv';
if all(options.GLMtask ~= "-")
    physconv.(options.GLMtask(1))=load(strcat(prefix,'_',options.GLMtask(1),suffix.(options.GLMtask(1)),'.txt'));
    physconv.(options.GLMtask(2))=load(strcat(prefix,'_',options.GLMtask(2),suffix.(options.GLMtask(2)),'.txt'));
    % Demean convolved physiological regressors
    physconv.(options.GLMtask(1))=physconv.(options.GLMtask(1))-mean(physconv.(options.GLMtask(1)));
    physconv.(options.GLMtask(2))=physconv.(options.GLMtask(2))-mean(physconv.(options.GLMtask(2)));
    % Load non-convolved versions too for visualization
    phys.(options.GLMtask(1))=load(strcat(prefix,'_',options.GLMtask(1),'.txt'));
    phys.(options.GLMtask(2))=load(strcat(prefix,'_',options.GLMtask(2),'.txt'));
    % Loading statement 
    fprintf(strcat("Loading: ",prefix,"_",options.GLMtask(1),suffix.(options.GLMtask(1)),".txt",...
    " and ",prefix,"_",options.GLMtask(2),suffix.(options.GLMtask(2)),".txt\n"))
end
% Load nonconvolved traces
phys.(options.Traces(1))=load(strcat(prefix,'_',options.Traces(1),'.txt'));
phys.(options.Traces(2))=load(strcat(prefix,'_',options.Traces(2),'.txt'));
% Loading statement
fprintf(strcat("Loading: ",prefix,"_",options.Traces(1),".txt",...
    " and ",prefix,'_',options.Traces(2),".txt\n"))
if all(options.moco ~= "-") && (length(options.moco) == 1)
    motion=load(options.moco);
    for i=1:size(motion,2)
        motion(:,i)=motion(:,i)-mean(motion(:,i));
    end
    if (size(motion,2) ~=6) && (size(motion,2) ~=2)
        warning('Motion file should have 2 (X and Y) or 6 (6DOF motion) columns. Will ignore.')
    end
    if ~isequal(size(motion,2),length(options.mocoLabel))
        error('Number of motion parameters (%d) does not equal number of mocoLabel entries (%d). Exiting.',...
            size(motion,2),length(options.mocoLabel))
    end
    if ~isequal(size(heatmap,2), size(motion,1))
        error('Length of motion traces does not match number of TRs (%d).', size(heatmap,2))
    end
else
    % Define dummy vector for passing if statement in GLM section
    motion=zeros(size(heatmap,2),2);
end
% Do more checks (equal lengths)
if options.GLMtask ~= "-"
    if ~isequal(size(heatmap,2), length(phys.(options.Traces(1))), length(phys.(options.Traces(2))), ...
            length(physconv.(options.GLMtask(1))), length(physconv.(options.GLMtask(2))))
        error('Length of traces does not match number of TRs (%d). Check fMRI data, and physiological traces.', size(heatmap,2))
    end
else
    if ~isequal(size(heatmap,2), length(phys.(options.Traces(1))), length(phys.(options.Traces(2))))
        error('Length of traces does not match number of TRs (%d). Check fMRI data, and physiological traces.', size(heatmap,2))
    end
end
% Create labels for phys data                       
% (can use this format for labels: [label.(options.Traces(1)) label.(option.Traces(2))]
label.(options.Traces(1))={['{\bf' convertStringsToChars(options.Traces(1))  '}']}; % Add units here later as optional input?
label.(options.Traces(2))={['{\bf' convertStringsToChars(options.Traces(2))  '}']};  
if all(options.GLMtask ~= "-")
    label.(options.GLMtask(1))={['{\bf' convertStringsToChars(options.GLMtask(1))  '}']};  
    label.(options.GLMtask(2))={['{\bf' convertStringsToChars(options.GLMtask(2))  '}']};  
end
label.CO2={'{\bfP_{ET}CO_{2}}','[mmHg]'};
label.HR={'{\bfHR}','[bpm]'};
label.RVT={'{\bfRVT}'};
label.O2={'{\bfO2}'};   
% Load stimulus (task timing vector)
if options.stim ~= "-"
    stim=load(options.stim);
    % Check that the stimulus vector is the same length as heatmap timeseries
    if ~isequal(length(stim), size(heatmap,2))
        error('Length of stimulus timing vector does not match number of TRs (%d).', size(heatmap,2))
    end
    stim_percent=0.1; % Percent of y-axis for stimulus vector to fill
    stim_gray=[150/255 150/255 150/255];
end
%% Load and organize slicewise motion correction parameters 
if all(options.moco ~= "-") && (length(options.moco) == 2)
    mocoX=load(options.moco(1));
    mocoY=load(options.moco(2));
    slicewise_moco_params_X=zeros(size(heatmap));
    slicewise_moco_params_Y=zeros(size(heatmap));
    for s=1:size(mocoX,2)
        slice=mocoX(3,s);
        temp_voxel_loc=find(voxelDir(:,3)==slice);
        for v=1:length(temp_voxel_loc)
            % Find where voxel directory is slice, and add moco timeseries
            idx=temp_voxel_loc(v);
            slicewise_moco_params_X(idx,:)=mocoX(4:end,slice);
            slicewise_moco_params_Y(idx,:)=mocoY(4:end,slice);
        end
    end
end
%% Plot data ordered by tissue
if (bySlice==0) && (options.plots==2)
    gmEnds=size(maskts{size(maskts,2)},2); % Size of GM graph portion
    figure('Name','By Tissue Type','Renderer', 'painters', 'Position', [50 1000 630 700])
    subplot(4,1,[3,4])
    imagesc(heatmap)
    set(gca,'YTickLabel',[],'FontSize',12); pbaspect([2 1 1])
    xlabel('{\bfTRs}')
    colormap gray
    caxis([c1 c2])
    % Draw white line b/t tissue types
%     line([0 nRow], [gmEnds+0.5 gmEnds+0.5],'Color','white','LineWidth',0.7) 
    % Add physio
    subplot(411); hold on
    plot(phys.(options.Traces(1)),'c','LineWidth',1.5); xlim([0 length(phys.(options.Traces(1)))])
    if options.stim ~= "-"
        ybound=ylim; 
        stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
        plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
    end
    ylabel(label.(options.Traces(1)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    set(gca,'XTickLabel',[],'FontSize',12); hold off
    subplot(412); hold on
    plot(phys.(options.Traces(2)),'g','LineWidth',1.5); xlim([0 length(phys.(options.Traces(1)))])
    if options.stim ~= "-"
        ybound=ylim;
        stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
        plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
    end
    ylabel(label.(options.Traces(2)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    set(gca,'XTickLabel',[],'FontSize',12); hold off
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
    if all(options.slices >= 0)
        tissueColorbar(idxs,:)=[];
    end
    % Plot tissue next to heatmap [x0 y0 width height]
    subplot('Position',[0.11 0.125 0.019 0.348])
    imagesc(tissueColorbar); colormap(gca,greengrayMap)
    caxis([0 max(tissueColorbar(:,1))]);
    set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
%     saveas(gcf,strcat(input_folder,'/',prefix,'_heatmap_byTissue.jpg'))
%     saveas(gcf,strcat(input_folder,'/',prefix,'_heatmap_byTissue_blur',options.PlotSmoothData,'.jpg')) % This is causing the input warning
end

%% Plot data ordered by slice (and vertebral level if opted for)
if bySlice==1 && (options.plots==2)
    figure('Name','By Vertebral Level','Renderer', 'painters', 'Position', [50 1000 630 700])
    subplot(4,1,[3,4])
    imagesc(heatmap)
    set(gca,'YTickLabel',[]); pbaspect([2 1 1])
    xlabel('{\bfTRs}'); set(gca,'FontSize',12)
%     ylabel('\leftarrow Inferior                     Superior \rightarrow')
    colormap gray
    caxis([c1 c2])
    % Add physio
    subplot(411); hold on
    plot(phys.(options.Traces(1)),'c','LineWidth',1.5); xlim([0 length(phys.(options.Traces(1)))])
    if options.stim ~= "-"
        ybound=ylim; 
        stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
        plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
    end
    ylabel(label.(options.Traces(1)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    set(gca,'XTickLabel',[],'FontSize',12); hold off
    subplot(412); hold on
    plot(phys.(options.Traces(2)),'g','LineWidth',1.5); xlim([0 length(phys.(options.Traces(1)))])
    if options.stim ~= "-"
        ybound=ylim; 
        stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
        plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
    end
    ylabel(label.(options.Traces(2)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    set(gca,'XTickLabel',[],'FontSize',12); hold off
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
%             line([nCol-10 nCol], [level-0.5 level-0.5], 'Color','y','LineWidth',1)
        end
        % Create a colorbar of the vertebral levels and add to figure
        vertebralLevels=ones(size(levChange,1)+1,3);
        vertebralLevels(1:end,1)=[levChange(:,1); levChange(end,1)+1];
        vertebralLevels(1:end,3)=[levChange(:,2); size(heatmap,1)];
        vertebralLevels(2:end,2)=vertebralLevels(1:end-1,3)+1;
        vertebralColorbar=zeros(size(heatmap,1),1);
        for v=1:size(vertebralLevels,1)
            vertebralColorbar(vertebralLevels(v,2):vertebralLevels(v,3))=vertebralLevels(v,1);
        end
        % Plot vertebral level next to heatmap
        subplot('Position',[0.11 0.125 0.019 0.348])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 maxLevel]) 
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    end
%     saveas(gcf,strcat(input_folder,'/',prefix,'_heatmap_bySlice.jpg'))
%     saveas(gcf,strcat(input_folder,'/',prefix,'_heatmap_bySlice_blur',string(options.PlotSmoothData),'.jpg'))
end

%% PHASE / FREQUENCY / POWER
% Fourier transform & calculations for phys data
for ph=1:2
    fs=1/TR; % 0.5 for TR=2
    phys_FT(:,ph)=fft(phys.(options.Traces(ph))); % Fourier transformed phys regressor
    n=length(phys_FT(:,ph));
    phys_power(:,ph)=abs(phys_FT(:,ph)).^2/n;                % non-zero centered power
    phys_powershift(:,ph)=abs(fftshift(phys_FT(:,ph))).^2/n; % zero-centered power
end
freq = (0:n-1)*(fs/n);      % non-zero centered frequency range
fshift=(-n/2:n/2-1)*(fs/n); % zero-centered frequency range
nyquist=fs/2;
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
if (options.plots==2)
    % Plot carpetplot frequency
    figure('Name','Power','Renderer', 'painters', 'Position', [50 1000 630 700])
    subplot(4,1,[3,4])
    imagesc(powerPlot(:,1:round(n/2))); set(gca,'xtick',[]); colormap pink

    % %%%%
    % size(powerPlot) 205 long
    % size(phys_power) 205 long
    % freq is 0 to .5 and is 205 long (#TRs)
    % %%%%
    caxis([0 0.25]); 
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
    %         for level=levChange(:,2)
    %             line([nCol/2-5 nCol/2], [level-0.5 level-0.5], 'Color','white')
    %         end
        end
        if useLevels==1
            % Plot vertebral level next to heatmap
            subplot('Position',[0.11 0.125 0.019 0.348])
            imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
            caxis([0 maxLevel]) 
        end
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
    title((options.Traces(1)))
    xlim([0 nyquist]) % limit plot by Nyquist frequency
    subplot(412)
    plot(freq,phys_power(:,2),'g','LineWidth',1.5)
    title((options.Traces(2)))
    xlim([0 nyquist])
end
freqmap=powerPlot(:,1:round(n/2)); % For function output
%% Run GLM
if (all(options.moco ~= "-")) && (all(options.GLMtask ~= "-"))
    if (all(options.moco ~= "-")) && (length(options.moco) == 2)
        % Make slicewise design matrices and demean each
        X=cell(size(heatmap,1),1);
        for s=1:size(heatmap,1)
            X{s}=[ones(size(heatmap,2),1) physconv.(options.GLMtask(1)) physconv.(options.GLMtask(2))...
                slicewise_moco_params_X(s,:)' slicewise_moco_params_Y(s,:)'];
            % Demean the design matrix
            for i=2:size(X{s},2)
                X{s}(:,i)=(X{s}(:,i)-mean(X{s}(:,i)))./range(X{s}(:,i));
            end
        end
        Bfit = zeros(nRow,size(X{1},2));
        tstats=zeros(nRow,size(X{1},2)-1);
        contrast=[0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
        % Run GLM with slicewise design matrices
        for row=1:nRow
            pin_mult=pinv(X{row})*pinv(X{row})';
            pinv_X=pinv(X{row});
            % Calculate beta coefficients
            Y=heatmap(row,:)';
            beta=pinv_X*Y;
            Bfit(row,:)=beta;
            % Calculate t-statistics
            r = Y-X{row}*beta;
            sumsq=r'*r; 
            DOF=size(X{row},1)-size(X{row},2); % Number of timepoints minus regressors
            rvar = sumsq/DOF;
            for i=1:size(contrast,1)
                c=contrast(i,:);
                tstats(row,i)=(c*beta)/sqrt(rvar*c*(pin_mult)*c');
            end
        end
    elseif length(options.moco) == 1
        X=[ones(size(heatmap,2),1) physconv.(options.GLMtask(1)) physconv.(options.GLMtask(2)) motion];
        % Demean the design matrix
        for i=2:size(X,2)
            X(:,i)=(X(:,i)-mean(X(:,i)))./range(X(:,i));
        end
        % Loop through each voxel Y individually when calculating beta=pinv(X)*Y
        % the length of beta will be the number of regressors. Then calculate t-stat
        Bfit = zeros(nRow,size(X,2));
        tstats=zeros(nRow,size(X,2)-1);
        pin_mult=pinv(X)*pinv(X)';
        pinv_X=pinv(X);
        % Define contrast vectors
        if size(motion,2)==6
            contrast=[0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0; 0 0 0 0 1 0 0 0 0;...
            0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1]; 
        elseif size(motion,2)==2
            contrast=[0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
        end
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
    end
else
    tstats=[];
end
%% Motion, phys, GLM plot
if all(options.moco ~= "-") && (all(options.GLMtask ~= "-")) && (options.plots==2)
    % Define colormaps
    greenMap = [zeros(256,1), linspace(0,1,256)', zeros(256,1)];
    cyanMap = [zeros(256,1), linspace(0,1,256)', linspace(0,1,256)'];
    magMap = [linspace(0,1,256)', zeros(256,1), linspace(0,1,256)'];
    rotationsMap = [linspace(0,1,256)', linspace(0,0.5686,256)', zeros(256,1)];
    translationsMap = [linspace(0,1,256)', zeros(256,1), zeros(256,1)];
    figure('Name','GLM, Convolved Regressors, and Motion Traces','Renderer', ...
        'painters','Position',[50 1 762 884])
%%% Positions:
    if size(motion,2)==6
        left=0.10; heatmap_bot=0.5131; heatmap_w=0.5147; heatmap_h=heatmap_w/2; phys_h=0.0904;
    elseif size(motion,2)==2
        left=0.10; heatmap_bot=0.3631; heatmap_w=0.7147; heatmap_h=heatmap_w/2; phys_h=0.0904;
    end
%%% Physio ( Position: [x0 y0 width height] )
    subplot('Position',[left heatmap_bot+heatmap_h+0.12 heatmap_w phys_h]); hold on
    plot(physconv.(options.GLMtask(1)),'c','LineWidth',1.5); xlim([0 length(physconv.(options.GLMtask(1)))])
    if options.stim ~= "-"
        ybound=ylim; 
        stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
        plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
    end
    xlim([0 length(physconv.(options.GLMtask(2)))]); set(gca,'xtick',[],'FontSize',12)
    ylabel(options.GLMtask(1),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold'); hold off
    subplot('Position',[left heatmap_bot+heatmap_h+0.02 heatmap_w phys_h]); hold on
    plot(physconv.(options.GLMtask(2)),'g','LineWidth',1.5); xlim([0 length(physconv.(options.GLMtask(2)))])
    if options.stim ~= "-"
        ybound=ylim; 
        stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
        plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
    end
    xlim([0 length(physconv.(options.GLMtask(2)))]); set(gca,'xtick',[],'FontSize',12)
    ylabel(options.GLMtask(2),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold'); hold off
%%% PLOT HEATMAP
    two_sd=2*std2(heatmap);
    subplot('Position', [left heatmap_bot heatmap_w heatmap_h])
    imagesc(heatmap)
    set(gca,'YTickLabel',[],'FontSize',11);
    caxis([c1 c2])
    xlabel('{\bfTRs}')
    colormap gray
    % Add small indicators of where the vertebral levels are
    if (bySlice==1) && (useLevels==1)
            levChange=(min(voxelDir(:,4)):max(voxelDir(:,4))-1); l=1;
            levChange=[levChange' zeros(1,range(voxelDir(:,4)))'];
            for i=2:length(voxelDir)
                if voxelDir(i-1,4)~=voxelDir(i,4)
                    levChange(l,2)=i;
                    l=l+1;
                end
            end
%             for level=levChange(:,2)
%                 line([nCol-10 nCol], [level-0.5 level-0.5], 'Color','white')
%             end
        % Plot vertebral level next to heatmap
        subplot('Position',[left-0.01 heatmap_bot 0.01 heatmap_h])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 maxLevel]) 
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    elseif bySlice==0
        % Plot tissue type next to heatmap [x0 y0 width height]
        subplot('Position',[left-0.01 heatmap_bot 0.01 heatmap_h])
        imagesc(tissueColorbar); colormap(gca,greengrayMap)
        caxis([0 max(tissueColorbar(:,1))]);
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    end
%%% Motion regressors
    if all(options.moco ~= "-") && (length(options.moco) == 2)
        subplot('Position',[left heatmap_bot-0.1631 heatmap_w phys_h]); hold on % X motion
        for i=1:size(mocoX,2)
            if all(options.slices >= 0)
                s=i-1; % Adjust slice # because of zero indexing
                % Include motion traces of slices being plotted
                if (s>topSlice) || (s<bottomSlice)
                    continue
                end
            end
            plot(mocoX(4:end,i),'LineWidth',0.25,'Color',[195/255 196/255 192/255])
        end
        plot(mean(mocoX(4:end,:),2),'LineWidth',1.5,'Color',[1 0.5686 0])
        if (abs(min(min(mocoX(4:end,:))))<1) && (abs(max(max(mocoX(4:end,:))))<1)
            ylim([-1 1])
        end
        if options.stim ~= "-"
            ybound=ylim; 
            stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
            plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
        end
        xlim([1 size(heatmap,2)]); set(gca,'FontSize',12,'XTickLabel',[],'xtick',[]);
        ylabel({'{\bfX Motion}','[mm]'},'rotation',90); hold off
        subplot('Position',[left heatmap_bot-0.2631 heatmap_w phys_h]); hold on % Y motion
        for i=1:size(mocoY,2)
            if all(options.slices >= 0)
                s=i-1; % Adjust slice # because of zero indexing
                % Include motion traces of slices being plotted
                if (s>topSlice) || (s<bottomSlice)
                    continue
                end
            end
            plot(mocoY(4:end,i),'LineWidth',0.25,'Color',[195/255 196/255 192/255])
        end
        plot(mean(mocoY(4:end,:),2),'LineWidth',1.5,'Color',[1 0 0])
        if (abs(min(min(mocoY(4:end,:))))<1) && (abs(max(max(mocoY(4:end,:))))<1)
            ylim([-1 1])
        end
        if options.stim ~= "-"
            ybound=ylim; 
            stim_height=range(ylim)*stim_percent; % Define height of binary stim vector to be 10% of y range
            plot(stim*stim_height+ybound(1),'LineWidth',2,'Color',stim_gray)
        end
        ylabel({'{\bfY Motion}','[mm]'},'rotation',90); xlim([1 size(heatmap,2)]); 
        set(gca,'FontSize',12,'XTickLabel',[],'xtick',[]); hold off
    elseif size(motion,2)==6
        motion_h=0.0705;
        subplot('Position',[left 0.02+(0.005+motion_h)*5 heatmap_w motion_h]) % Rx 
        plot(motion(:,1),'Color',[1 0.5686 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
        ylabel(options.mocoLabel(1),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
        subplot('Position',[left 0.02+(0.005+motion_h)*4 heatmap_w motion_h]) % Ry 
        plot(motion(:,2),'Color',[1 0.5686 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
        ylabel(options.mocoLabel(2),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
        subplot('Position',[left 0.02+(0.005+motion_h)*3 heatmap_w motion_h]) % Rz 
        plot(motion(:,3),'Color',[1 0.5686 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
        ylabel(options.mocoLabel(3),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
        subplot('Position',[left 0.02+(0.005+motion_h)*2 heatmap_w motion_h]) % Tx 
        plot(motion(:,4),'Color',[1 0 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
        ylabel(options.mocoLabel(4),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
        subplot('Position',[left 0.02+(0.005+motion_h)*1 heatmap_w motion_h]) % Ty 
        plot(motion(:,5),'Color',[1 0 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
        ylabel(options.mocoLabel(5),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
        subplot('Position',[left 0.02 heatmap_w motion_h]) % Tz 
        plot(motion(:,6),'Color',[1 0 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); 
        ylabel(options.mocoLabel(6),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    elseif size(motion,2)==2
        motion_h=0.0705; % phys_h is 0.0904
        subplot('Position',[left heatmap_bot-0.1631 heatmap_w phys_h]) % X 
        plot(motion(:,1),'Color',[1 0.5686 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); ylim([-1.2 1.2])
        ylabel(options.mocoLabel(1),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
        subplot('Position',[left heatmap_bot-0.2631 heatmap_w phys_h]) % Y 
        plot(motion(:,2),'Color',[1 0 0],'LineWidth',1.5); xlim([0 length(motion)]); set(gca,'xtick',[],'FontSize',12); ylim([-1.2 1.2])
        ylabel(options.mocoLabel(2),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    end
%%% GLM tstats ( Position: [x0 y0 width height] )
    tstat_w=0.0419;
    subplot('Position',[left+heatmap_w-0.04+tstat_w heatmap_bot 0.0419 heatmap_h]) % CO2
    imagesc(abs(tstats(:,1))); set(gca,'xtick',[],'ytick',[]); colormap(gca,cyanMap); caxis([0 5])
    subplot('Position',[left+heatmap_w-0.04+tstat_w*2 heatmap_bot 0.0419 heatmap_h]) % HR
    imagesc(abs(tstats(:,2))); set(gca,'xtick',[],'ytick',[]); colormap(gca,greenMap); caxis([0 5])
    if size(motion,2)==6
        subplot('Position',[left+heatmap_w-0.04+tstat_w*3 heatmap_bot 0.1257 heatmap_h]) % MOTION x6
        imagesc(abs(tstats(:,3:5))); set(gca,'xtick',[],'ytick',[]);  colormap(gca,rotationsMap); caxis([0 5])
        title('t-statistics','HorizontalAlignment','right','FontSize',12)
        subplot('Position',[left+heatmap_w-0.04+tstat_w*6 heatmap_bot 0.1257 heatmap_h]) % MOTION x6
        imagesc(abs(tstats(:,6:8))); set(gca,'xtick',[],'ytick',[],'FontSize',12);  colormap(gca,translationsMap); caxis([0 5])
    elseif size(motion,2)==2
        subplot('Position',[left+heatmap_w-0.04+tstat_w*3 heatmap_bot 0.0419 heatmap_h]) % X
        imagesc(abs(tstats(:,3))); set(gca,'xtick',[],'ytick',[]);  colormap(gca,rotationsMap); caxis([0 5])
        title('t-statistics','HorizontalAlignment','right','FontSize',12)
        subplot('Position',[left+heatmap_w-0.04+tstat_w*4 heatmap_bot 0.0419 heatmap_h]) % Y
        imagesc(abs(tstats(:,4))); set(gca,'xtick',[],'ytick',[],'FontSize',12);  colormap(gca,translationsMap); caxis([0 5])
    end
end
%% Reorganize data and plot by t-statistic magnitude - GLMtstats(1)
if all(options.moco ~= "-") && (all(options.GLMtask ~= "-")) && (options.plots==2)
    % GLMtstats(1) heatmap
    temp_sorter=[abs(tstats(:,1)) heatmap];
    heatmap_1_sort=sortrows(temp_sorter,1,'descend');
    heatmap_1_sort=heatmap_1_sort(:,2:end);
    % GLMtstats(1) tstats
    temp_sorter=[abs(tstats(:,1)) tstats];
    tstats_1_sort=sortrows(temp_sorter,1,'descend');
    tstats_1_sort=tstats_1_sort(:,2:end);
%%% Positions:
    heatmap_w=0.75; heatmap_h=heatmap_w/2; tstat_w=0.0519; heatmap_bot=0.12;
    % Physio (plotting nonconvolved even though convolved were used for the GLM...)
    figure('Name',strcat("Plot data by ",options.GLMtask(1)," t-statistic magnitude"),'Renderer', 'painters', 'Position', [50 1000 800 700])
    subplot('Position',[left heatmap_h+0.35 heatmap_w heatmap_h/2]); hold on
    plot(phys.(options.GLMtask(1)),'c','LineWidth',1.5); xlim([0 length(phys.(options.GLMtask(1)))]); set(gca,'xtick',[],'FontSize',12)
    if options.stim ~= "-"
        ybound=ylim;
        ylim(ybound); stim_vec=stim+ybound(1);
        plot(stim_vec,'LineWidth',2,'Color',[195/255 196/255 192/255])
    end
    ylabel(label.(options.GLMtask(1)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right'); hold off
    subplot('Position',[left heatmap_h+0.15 heatmap_w heatmap_h/2])
    plot(phys.(options.GLMtask(2)),'g','LineWidth',1.5); xlim([0 length(phys.(options.GLMtask(2)))]); set(gca,'xtick',[],'FontSize',12)
    ylabel(label.(options.GLMtask(2)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    % Heatmap
    subplot('Position',[left heatmap_bot heatmap_w heatmap_h])
    imagesc(heatmap_1_sort)
    set(gca,'YTickLabel',[],'FontSize',12);
    caxis([c1 c2])
    xlabel('{\bfTRs}')
    colormap gray
    % GLM tstats ( Position: [x0 y0 width height] )
    subplot('Position',[left+heatmap_w-0.05+tstat_w*1 heatmap_bot tstat_w heatmap_h]) % CO2
    imagesc(abs(tstats_1_sort(:,1))); set(gca,'xtick',[],'ytick',[]); colormap(gca,cyanMap); caxis([0 5]); title(options.GLMtask(1))
    subplot('Position',[left+heatmap_w-0.05+tstat_w*2 heatmap_bot tstat_w heatmap_h]) % HR
    imagesc(abs(tstats_1_sort(:,2))); set(gca,'xtick',[],'ytick',[]); colormap(gca,greenMap); caxis([0 5]); title(options.GLMtask(2))
end
%% Reorganize data and plot by t-statistic magnitude - GLMtstats(2)
if all(options.moco ~= "-") && (all(options.GLMtask ~= "-")) && (options.plots==2)
    % HR heatmap
    temp_sorter=[abs(tstats(:,2)) heatmap];
    heatmap_2_sort=sortrows(temp_sorter,1,'descend');
    heatmap_2_sort=heatmap_2_sort(:,2:end);
    % HR tstats
    temp_sorter=[abs(tstats(:,2)) tstats];
    tstats_2_sort=sortrows(temp_sorter,1,'descend');
    tstats_2_sort=tstats_2_sort(:,2:end);
    clear temp_sorter
    % Physio
    figure('Name',strcat("Plot data by ",options.GLMtask(2)," t-statistic magnitude"),'Renderer', 'painters', 'Position', [50 1000 800 700])
    subplot('Position',[left heatmap_h+0.35 heatmap_w heatmap_h/2]); hold on
    plot(phys.(options.GLMtask(1)),'c','LineWidth',1.5); xlim([0 length(phys.(options.GLMtask(1)))]); set(gca,'xtick',[],'FontSize',12)
    if options.stim ~= "-"
        ybound=ylim;
        ylim(ybound); stim_vec=stim+ybound(1);
        plot(stim_vec,'LineWidth',2,'Color',[195/255 196/255 192/255])
    end
    ylabel(label.(options.GLMtask(1)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right'); hold off
    subplot('Position',[left heatmap_h+0.15 heatmap_w heatmap_h/2])
    plot(phys.(options.GLMtask(2)),'g','LineWidth',1.5); xlim([0 length(phys.(options.GLMtask(2)))]); set(gca,'xtick',[],'FontSize',12)
    ylabel(label.(options.GLMtask(2)),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    % Heatmap
    subplot('Position',[left heatmap_bot heatmap_w heatmap_h])
    imagesc(heatmap_2_sort)
    set(gca,'YTickLabel',[],'FontSize',12);
    caxis([c1 c2])
    xlabel('{\bfTRs}')
    colormap gray
    % GLM tstats ( Position: [x0 y0 width height] )
    subplot('Position',[left+heatmap_w-0.05+tstat_w*1 heatmap_bot tstat_w heatmap_h]) % CO2
    imagesc(abs(tstats_2_sort(:,1))); set(gca,'xtick',[],'ytick',[]); colormap(gca,cyanMap); caxis([0 5]); title(options.GLMtask(1))
    subplot('Position',[left+heatmap_w-0.05+tstat_w*2 heatmap_bot tstat_w heatmap_h]) % HR
    imagesc(abs(tstats_2_sort(:,2))); set(gca,'xtick',[],'ytick',[]); colormap(gca,greenMap); caxis([0 5]); title(options.GLMtask(2))
end

%% Motion trace plot
if all(options.moco ~= "-") && (options.plots==2)
    figure('Name','Motion','Renderer', 'painters', 'Position', [50 1000 683 700])
    subplot(4,1,[1,2])
    imagesc(heatmap); xlabel('{\bfTRs}');
    set(gca,'YTickLabel',[],'FontSize',12); pbaspect([2 1 1])
    caxis([c1 c2]); colormap gray
    if length(options.moco) == 2
        subplot(413) % Plot X
        hold on
        for i=1:size(mocoX,2)
            if all(options.slices >= 0)
                s=i-1; % Adjust slice # because of zero indexing
                % Include motion traces of slices being plotted
                if (s>topSlice) || (s<bottomSlice)
                    continue
                end
            end
            plot(mocoX(4:end,i),'LineWidth',0.25,'Color',[195/255 196/255 192/255])
        end
        plot(mean(mocoX(4:end,:),2),'LineWidth',1.5,'Color',[1 0.5686 0])
        if (abs(min(min(mocoX(4:end,:))))<1) && (abs(max(max(mocoX(4:end,:))))<1)
            ylim([-1 1])
        end
        xlim([1 size(heatmap,2)]); set(gca,'FontSize',12,'XTickLabel',[]);
        ylabel({'{\bfX Motion}','[mm]'},'rotation',90); hold off
        subplot(414) % Plot Y
        hold on
        for i=1:size(mocoY,2)
            if all(options.slices >= 0)
                s=i-1; % Adjust slice # because of zero indexing
                % Include motion traces of slices being plotted
                if (s>topSlice) || (s<bottomSlice)
                    continue
                end
            end
            plot(mocoY(4:end,i),'LineWidth',0.25,'Color',[195/255 196/255 192/255])
        end
        plot(mean(mocoY(4:end,:),2),'LineWidth',1.5,'Color',[1 0 0])
        if (abs(min(min(mocoY(4:end,:))))<1) && (abs(max(max(mocoY(4:end,:))))<1)
            ylim([-1 1])
            fprintf('hi')
        end
        xlim([1 size(heatmap,2)]); set(gca,'FontSize',12,'XTickLabel',[]); 
        ylabel({'{\bfY Motion}','[mm]'},'rotation',90); hold off
    elseif size(motion,2)==6
        subplot(413) % Plot rotations
        hold on
        plot(motion(:,1),'Color',[1 0.5686 0],'LineWidth',1)
        plot(motion(:,2),'Color',[1 0.5686 0],'LineWidth',1)
        plot(motion(:,3),'Color',[1 0.5686 0],'LineWidth',1)
        ylim([-0.06 0.06]); xlim([1 length(motion)]); set(gca,'FontSize',12,'XTickLabel',[])
        if all(options.mocoLabel==["Rx" "Ry" "Rz" "Tx" "Ty" "Tz"])
            ylabel({'{\bfRotations}','[rads]'},'rotation',90)
        else
            ylabel(append(options.mocoLabel(1),",",options.mocoLabel(2),",",options.mocoLabel(3)),'FontWeight','bold')
        end
        hold off
        subplot(414) % Plot translations
        hold on
        plot(motion(:,4),'Color',[1 0 0],'LineWidth',1.5)
        plot(motion(:,5),'Color',[1 0 0],'LineWidth',1.5)
        plot(motion(:,6),'Color',[1 0 0],'LineWidth',1.5)
        ylim([-6 6]); xlim([1 length(motion)])
        if all(options.mocoLabel==["Rx" "Ry" "Rz" "Tx" "Ty" "Tz"])
            ylabel({'{\bfTranslations}','[mm]'},'rotation',90)
        else
            ylabel(append(options.mocoLabel(1),",",options.mocoLabel(2),",",options.mocoLabel(3)),'FontWeight','bold')
        end
        xlabel('{\bfTRs}'); set(gca,'FontSize',12)
        hold off
    elseif size(motion,2)==2
        % Plotting only 2 motion parameters 
        % summary tsv X and Y sct_fmri_moco outputs
        subplot(413) % Plot X average
        hold on
        plot(motion(:,1),'Color',[1 0.5686 0],'LineWidth',1)
        ylim([-1 1]); xlim([1 length(motion)]); set(gca,'XTickLabel',[]);
        ylabel({'{\bfX}','[mm]'},'rotation',90)
        subplot(414) % Plot Y average
        plot(motion(:,2),'Color',[1 0 0],'LineWidth',1.5)
        ylim([-1 1]); xlim([1 length(motion)])
        ylabel({'{\bfY}','[mm]'},'rotation',90)
        xlabel('{\bfTRs}'); set(gca,'FontSize',12)
        hold off
    end
    if (bySlice==1) && (useLevels==1)
        % Plot vertebral level next to heatmap
        subplot('Position',[0.111 0.5482 0.019 0.3768])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 maxLevel]) 
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    elseif bySlice==0
%         if useLevels==1
%             % Draw line b/t tissue types ---- why is this here I'm confused
%             subplot(4,1,[1,2])
%             line([0 nRow], [gmEnds+0.5 gmEnds+0.5],'Color','white','LineWidth',0.7) 
%         end
        % Plot tissue next to heatmap [x0 y0 width height]
        subplot('Position',[0.111 0.5482 0.019 0.3768])
        imagesc(tissueColorbar); colormap(gca,greengrayMap)
        caxis([0 max(tissueColorbar(:,1))]);
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[])
    end
%     saveas(gcf,strcat(input_folder,'/',prefix,'_heatmap_motion_blur',options.PlotSmoothData,'.jpg'))
end
%% DVARS plot
if (options.plots==2)
    dvars=load('dvars.txt');
    dvars(1)=NaN;
    figure('Name','DVARS', 'Renderer', 'painters','Position', [50 1000 630 500])
    subplot(311)
    plot(1:n,dvars,'b','LineWidth',1.5); xlim([0 length(dvars)]); ylim([0 max(dvars)])
    ylabel('DVARS','rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold')
    set(gca,'XTickLabel',[],'xtick',[])
    subplot(3,1,[2,3])
    imagesc(heatmap)
    set(gca,'YTickLabel',[]); pbaspect([2 1 1])
    xlabel('{\bfTRs}'); caxis([c1 c2]); colormap gray
    if (bySlice==1) && (useLevels==1)
        % Plot vertebral level next to heatmap
        subplot('Position',[0.11 0.124 0.019 0.488])
        imagesc(vertebralColorbar); colormap(gca,blueLightBlueMap)
        caxis([0 maxLevel]) 
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[],'FontSize',12)
    elseif bySlice==0
        % Plot tissue next to heatmap [x0 y0 width height]
        subplot('Position',[0.11 0.124 0.019 0.488])
        imagesc(tissueColorbar); colormap(gca,greengrayMap)
        caxis([0 max(tissueColorbar(:,1))]);
        set(gca,'XTickLabel',[],'xtick',[],'YTickLabel',[],'ytick',[],'FontSize',12)
    end
end
%% Write out files if requested saveas(gcf,strcat(input_folder,'/',prefix,'_heatmap_byTissue.jpg'))
if write_out==1
    if bySlice==1
        order='bySlice';
    elseif bySlice==0
        order='byTissue';
    end
    writematrix(heatmap,[input_folder '/' prefix order '_heatmap.txt'],'Delimiter','tab')
    writematrix(voxelDir,[input_folder '/' prefix order 'voxelDirectory.txt'],'Delimiter','tab')
    if tstats
        writematrix(voxelDir,[input_folder '/' prefix order '_tstats.txt'],'Delimiter','tab')
    end
end
fprintf(['\nFile(s) saved to: ' input_folder '\n'])
%% Remove added paths
rmpath(input_folder)
rmpath(trace_loc)

fprintf('\n...done!\n')
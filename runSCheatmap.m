%% Test example data (add input_folder then run this section)
% Add the full path to your output folder from x.heatmapPrep as a string
% E.g., '~/Downloads/spinalcordplot-main/exampleData/heatmap_output'
input_folder=

% Confirm phys_loc and mLoc paths are correct
phys_loc='~/Downloads/spinalcordplot-main/exampleData/phys';
mLoc='~/Downloads/spinalcordplot-main/exampleData/motionTraces.txt';

% Adjust this for tissue type / longitudinal level organization (1 or 0)
bySlice=1;

% Adjust this for whether to indicate vertebral levels on plot (1 or 0)
useLevels=1;

% Adjust this to indicate whether you would like files written to
% input_folder (1 or 0)
write_out=0;

% Don't edit between lines ------------------------------------------------
TR=2;
prefix='sub-example';
% -------------------------------------------------------------------------
%% Run this section to visualize just the heatmap
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,[],"basic",1);

%% Run this section to visualize basic heatmap figures
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,phys_loc);

%% Run this section to visualize other physiological traces
traces=["RVT" "O2"];
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,phys_loc,"Traces",traces);

%% Run this section to adjust heatmap contrast
% Edit variable 'c' to adjust contrast (default: 0.4)
c=0.1; 
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,phys_loc,"cBound",c);

%% Run this section to also visualize motion data
% Variable 'mLabel' used to change motion labels if input motion file is 
% not in the default organization. Example data is provided in the default
% organization.
mLabel=["Rx" "Ry" "Rz" "Tx" "Ty" "Tz"];
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,phys_loc,"moco",mLoc,"mocoLabel",mLabel);

%% Run this section to run / visualize a simple GLM
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,phys_loc,"moco",mLoc,"GLMtask",["CO2" "HR"]);

%% Run this section to plot smoothed data
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,phys_loc,"PlotSmoothData",1);


%% Test example data
% Add the full path to your output folder from x.heatmapPrep as a string
% E.g., ''~/Downloads/spinalcordplot-main/heatmap_output'
input_folder=

% Adjust this for tissue type / longitudinal level organization (1 or 0)
bySlice=1;

% Adjust this for whether to indicate vertebral levels on plot (1 or 0)
useLevels=1;

% Don't edit below here ---------------------------------------------------
TR=2;
phys_prefix='sub-example';
phys_loc='~/Downloads/spinalcordplot-main/exampleData/phys';
moco_loc='~/Downloads/spinalcordplot-main/exampleData/motionTraces.txt';
[heatmap,voxelDir,tstats]=SCheatmap(input_folder,bySlice,useLevels,TR,phys_prefix,phys_loc,moco_loc);
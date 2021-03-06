%% Test example data (add input_folder then run this section)
% Add the full path to your output folder from x.heatmapPrep as a string
% E.g., '~/Downloads/spinalcordplot-main/exampleData/heatmap_output'
input_folder=
% Add the full path to the location of traces as a string
% E.g., '~/Downloads/spinalcordplot-main/exampleData/traces'
trace_loc=


% Paths to motion traces
mocoX=strcat(trace_loc,"/sub-example_moco-X.txt");
mocoY=strcat(trace_loc,"/sub-example_moco-Y.txt");

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
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,[],"plots",1);

%% Run this section to visualize basic heatmap figures
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc);

%% Run this section to visualize other physiological traces
traces=["RVT" "O2"];
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc,"Traces",traces);

%% Run this section to adjust heatmap contrast
% Edit variable 'c' to adjust contrast (default: 0.4)
c=0.1; 
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc,"cBound",c);

%% Run this section to also visualize motion data
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc,"moco",[mocoX,mocoY]);

%% Run this section to run & visualize a simple GLM & stimulus vector
stim=strcat(trace_loc,"/sub-example_stim.txt");
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc,"moco",[mocoX,mocoY],"GLMtask",["CO2" "HR"],"stim",stim);

%% Run this section to plot smoothed data
[heatmap,freqmap,voxelDir]=SCheatmap(input_folder,write_out,bySlice,useLevels,TR,prefix,trace_loc,"PlotSmoothData",1);


## Spinal Cord fMRI heatmaps
Repository containing preliminary code and example data to create heatmaps for the spinal cord.

---



#### Guide to this repository:
There are two scripts used to create the heatmaps. The bash shell script `x.heatmapPrep` requires the Spinal Cord Toolbox, AFNI, and FSL to be installed on your machine. The second script, `SCheatmap.m` is a MATLAB script, and requires MATLAB to be installed.

#### How to use
1. To test out this code, first clone or download the repository onto your machine.
2. Download this Box folder titled [*exampleData*](https://northwestern.box.com/s/e6apnxx1lly1r0k9uj8stp32esqvquq8) (or https://northwestern.box.com/s/e6apnxx1lly1r0k9uj8stp32esqvquq8) and add to the repository folder on your local machine.
3. Open terminal and navigate to where you stored the repository folder on your machine.
4. **Run the shell script:** Type `bash x.heatmapPrep -i ~/Downloads/spinalcordplot-main/exampleData -o heatmap_output` into your terminal. Adjust if the input path does not match the location / name of the folder. For help information, just type `bash x.heatmapPrep -h` into your terminal.
5. **Run the MATLAB script:** Navigate within MATLAB to the downloaded / cloned repository folder and open the file `runSCheatmap.m`. This script should run the function `SCheatmap()`, and output various plots. For help information, type `help SCheatmap` into the MATLAB command window.


*Notes*:
- This code has only been tested on a computer running MacOS Catalina, MATLAB 2019b.
- Using an external monitor might alter some of the output figure aspect ratios.
- Using bash vs. zsh shells may cause terminal messages with color to appear just as text. This is a minor bug. Program will still run.

*README last updated: 12/21/2020*

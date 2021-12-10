## SpinalCordPlot: Spinal Cord fMRI heatmaps
Repository containing preliminary code and example data to create heatmaps for the spinal cord.

---

[![Twitter URL](https://img.shields.io/twitter/follow/KJHemm?style=social)](https://twitter.com/KJHemm)
[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/BrightLab-ANVIL/spinalcordplot/blob/main/LICENSE)
<br>

#### Guide to this repository:
There are two scripts used to create the heatmaps. The bash shell script `x.heatmapPrep` requires the Spinal Cord Toolbox, AFNI, and FSL to be installed on your machine. The second script, `SCheatmap.m` is a MATLAB script, and requires MATLAB to be installed.

#### How to use
1. To test out this code, first clone or download the repository onto your machine.
2. Download this Box folder titled [*exampleData*](https://northwestern.box.com/s/hhao1e7n8z4h7p90uks1d4z5is4ehcl2) (or https://northwestern.box.com/s/hhao1e7n8z4h7p90uks1d4z5is4ehcl2) and add to the repository folder on your local machine.
3. Open terminal and navigate to where you stored the repository folder on your machine.
4. **Run the shell script:** Type `zsh x.heatmapPrep -i ~/Downloads/spinalcordplot-main/exampleData/label -f ~/Downloads/spinalcordplot-main/exampleData/func.nii.gz -o heatmap_output -s` into your terminal. Adjust the input path for the proper location on your computer. For this example, not including CSF is recommended. For help information, just type `zsh x.heatmapPrep -h` into your terminal.
5. There will be one manual step while running the shell script. An fsleyes window with a spinal cord mask will pop-up. In edit mode, erase the top and bottom slices in the mask and save / overwrite the file, then close the window.
6. **Run the MATLAB script:** Navigate within MATLAB to the downloaded / cloned repository folder and open the file `runSCheatmap.m`. Each section within this script should run the function `SCheatmap()`, and output various plots. Run sections individually. For help information, type `help SCheatmap` into the MATLAB command window.


*Notes*:
- Latest testing on Mac running MacOS Catalina, MATLAB 2019b, SCT 5.2.0, AFNI 21.1.01, FSL 6.0.
- Using an external monitor might alter some of the output figure aspect ratios.
- Using bash vs. zsh shells may cause terminal messages with color to appear just as text. This is a minor bug. Program will still run.

*Citation*:
K. J. Hemmerling and M. G. Bright, "A visualization tool for assessment of spinal cord functional magnetic resonance imaging data quality," 2021 43rd Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC), 2021, pp. 3391-3394, doi: 10.1109/EMBC46164.2021.9630903.

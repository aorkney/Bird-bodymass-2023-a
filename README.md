Figure 1: Figure 1 is a phylogram illustrating the phylogenetic structure of body mass and an ancestral state reconstruction over the 149 bird species used in this study. 
There are also a collection of subplots which illustrate increasing variance in skeletal proportions in more massive bird taxa. Allometric trends have been parsed-out prior to this analysis. 
These subplots demonstrate that some modular regions, such as the leg, exhibit a proximo-distal gradient of increasing variance- consistent with the sequence of embryogenesis. 
However, some modular regions- such as the wing- do not show this pattern. 
Run the script 'Fig_var_stru_08_08_2023.R' to produce this figure.

Figure 2: Figure 2 is a collection of line plots which will help the user visualise how the strength of evolutionary integration between bones within different regions/modules of the avian skeleton change as avian body mass increases. This is achieved by generating many taxonomic subsets of birds with different mean body mass values, and repeating analyses to determin evolutionary integration.

a) Run script 'Residual_Csize_10_08_2022.R' to generate a dataset of avian bone sizes (centroid size computed from landmark constellations), corrected for patterns of allometric scaling (size-dependent scaling; this must be removed because the role of body mass is being investigated and size-dependent scaling could obfuscate such patterns).

b) Run script 'Integration_v_body_mass_11_08_2022.R' to run the analysis.

c) Run script 'Plot_integration_v_body_mass_14_11_2022.R' to plot the analysis output (Figure 2).

# LC connectivity analysis

Analysis of Lobula Columnar (LC) neuron connectivity using [Drosophila Hemibrain connectivity data](https://doi.org/10.7554/eLife.57443), exported from [neuPrint](http://neuprint.janelia.org).


## LC_T-bars
Matlab script to visualize unique T-bars for each LC neuron using xyz coordinates from the hemibrain dataset. Script also requires Matlab Statistics and Machine Learning Toolbox to perform kmeans clustering.


## LC_VLPN_connectivity
Analysis of ventrolateral protocerebrum neurons (VLPN) that directly integrate from LC outputs in the optic glomeruli. Colormap used is included and originally created by [Thyng et al 2016](http://dx.doi.org/10.5670/oceanog.2016.66).


## Similarity_Distance
Analysis of number of downstream neurons (VLPN) integrating from optic glomeruli versus the distance or similarity between glomeruli. Similarity is computed using peak calcium responses provided in /Input/lc_round1_compiled_2019.08.06.mat


## T2_analysis
Jupyter notebook for exporting putative T2 neuron connectivity to LC neurons.
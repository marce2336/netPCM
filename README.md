# netPCM (network-based phenotypic characterization of mutants)

## Description

Collection of scripts used to reproduce the simulations described in the publication "Integration of lipidomics data into a metabolic model accurately characterizes the effects of mutations in _Arabidopsis thaliana’s_ lipid metabolism". Also included are the source data and scripts necessary to construct the figures shown in the publication.

Last update: 2024-04-19

This repository is administered by Sandra M. Correa (Cordoba@mpimp-golm.mpg.de, sandra.correa.cordoba.1@uni-potsdam.de), Bioinformatics Group, University of Potsdam - Germany


## Installation

### Requirements

- A functional MATLAB installation (version R2019b it's recommended due to its compatibility with IBM CPLEX solver)

- The COBRA toolbox for MATLAB <https://opencobra.github.io/cobratoolbox/stable/installation.html>

- IBM CPLEX Studio IDE 12.10.0 <https://www.ibm.com/products/ilog-cplex-optimization-studio>, <http://ibm.biz/CPLEXonAI>

- R <https://cran.r-project.org/> and RStudio <https://posit.co/downloads/>

### Dependencies (Recommended Software)

- The COBRA toolbox: the instructions to install the COBRA toolbox version for MATLAB can be found in <https://opencobra.github.io/cobratoolbox/stable/installation.html>

- IBM CPLEX solver: the solver and corresponding license can be obtained from <https://www.ibm.com/products/ilog-cplex-optimization-studio>.

- R and RStudio: to create the figures shown in the publication, the following libraries must be installed in advance

        install.packages("readxl")
        install.packages("pheatmap")
        install.packages("ggplot2")
        install.packages("viridis")

- Flux samples: the sets of flux samples generated after implementing the Artificially Centered Hit-and-Run (ACHR) algorithm can be downloaded as .zip file from the figshare repository: https://doi.org/10.6084/m9.figshare.25612167.v1. Once the download is complete, unzip the file and move all its contents to the folder located in the path: '..\SamplingResults\Random_samples'. To reproduce the simulations from scratch, ignore this step and run the 'processSamplingResults' function as described below. If the last option is selected, keep in mind that each simulation round can take up to three hours depending on the characteristics of the equipment used.


## Usage

1. Initialize MATLAB with administrator privileges

2. Install the software listed in section 'Requirements'

3. Download the contents of this repository to a known location. It is optional to add the path with sub-folders.

4. Initialize the CobraToolbox:

          initCobraToolbox()

5. Change the cobra solver to 'cplexlp':

          changeCobraSolver('cplexlp')

6. Type in the Command Window any of the commands suggested below



### Perform flux sampling for Arabidopsis T-DNA lines

#### Description:
The sampling procedure and the computation of optimal growth values for Arabidopsis T-DNA lines can be done by running the function 'processSamplingResults'. The script and dependencies are located in the path '..\SimulateGrowthMTs'.

#### Usage:
          [statistics, processedData, growthMTs] = processSamplingResults(flagDataset, flagPruning)

#### Input:
- flagDataset	 -  ('all') Default. Load information for the complete set of mutants.
                  ('mpi') Load information only for mutants measured in-house.
                  ('lusk') Load information only for mutants published by Lusk et al. (doi:https://doi.org/10.1093/pcp/pcac088).

- flagPruning	 -  (1) Default. Eliminate from lipid profiles the species that are not included in the Plant Lipid Module (doi:https://doi.org/10.1038/s41467-023-40644-9).
                  (0) Keep all species listed in the lipidome profiles

#### Example:
          [statistics, processedData, growthMTs] = processSamplingResults('all', 1)

#### Output:
- growthMTs - Optimal growth values computed for Arabidopsis T-DNA lines

- processedData - struct with the fields:
  - processedProfiles  - fold change computed for lipid relative abundances (m = RA-T-DNA/RA-WT)
  - referenceFluxes	 - incoming reactions and fluxes for each of the lipid species measured
  - modelIrrev  - irreversible version of the model built for Arabidopsis wild-type Col-0 plants
  - flagRelaxation  - flag that indicates if simulations are performed relaxing bounds with standard deviation (1 - bounds are relaxed)
  - modelsMT_w - models built for Arabidopsis T-DNA lines constrained with flux sums and lipid abundances

- statistics - struct with the following fields:
  - meanMTs_wR  - mean flux values for reactions underlying the mutated locus
  - oneSampleTtest_wR  - results of one-sample t-Test performed for mean flux values of T-DNA lines
  - twoSampleTtest_wR  - results of two-sample t-Test for mean flux values of wild type Col-0 and T-DNA lines
  - tsTtestUnique_wR - results two-sample t-Test for reactions catalyzed by unique enzyme product
  - tsTtestIsoenzymes_wR - results two-sample t-Test for reactions catalyzed by isoenzymes

#### Description files generated:
The files generated during the simulations are stored in several sub-folders located in the path '..\SamplingResults'.
- FilesSamplingProcedure - set of files that contain the 10.000 samples obtained for the T-DNA lines and wild-type Col-0 plants.
- filesWarmupPoints - warmup point sets created for hit-and-run sampling by combining orthogonal and random points.
- FVA_results - flux vectors with results from flux variability analysis performed for T-DNA lines and wild-type Col-0 plants.
- Models_sampling - models built for T-DNA lines and wild-type Col-0 plants.
- Random_samples - samples for T-DNA lines and wild-type Col-0 plants loaded for analysis.
- ReducedModel - version of the models built for T-DNA lines and wild-type Col-0 plants eliminating reactions not carrying flux.
- SampledFluxes - flux samples for the reactions underlying the mutated loci for T-DNA lines and wild-type Col-0 plants.
- Ttest -  results of mean of fluxes ('Results_classification_mutants.xlsx'), and t-Test analysis ('Results_statistic_tests.xlsx') performed to compare mean of fluxes among T-DNA lines and wild-type Col-0 plants.


### Perform sampling for Arabidopsis T-DNA lines adding constraints to lipid pools

#### Description:
The script 'constrainUnmeasuredSpecies' is used to identify mutants whose simulation results do not match the phenotypic descriptions available in the literature. This is followed by the inclusion in the model of an additional set of constraints on the level of unmeasured lipid pools. Subsequently, the algorithm performs the sampling procedure again. Before running this script, the flux sampling must be performed for all T-DNA lines. To do this, the 'processSamplingResults' function must be executed as explained above.

#### Usage:
          [statistics, mutantsClasses_wR, listOutliers, growthMTs, processedData] = constrainUnmeasuredSpecies(flagDataset, flagSave, flagPruning)

#### Input:
- flagDataset  -  ('lusk') Perform random sampling for Lusk dataset.
                  ('mpi') Perform random sampling for in-house dataset.

- flagSave  -  (0) Default. Don´t save output from processing sampling results.
               (1) Save output from processing sampling results.

- flagPruning	 -  (1) Default. Eliminate from lipid profiles the species that are not included in the Plant Lipid Module (doi:https://doi.org/10.1038/s41467-023-40644-9).
                  (0) Keep all species listed in the lipidome profiles.

#### Examples:
          [statistics, mutantsClasses_wR, listOutliers, growthMTs, processedData] = constrainUnmeasuredSpecies('lusk', 1, 1)
          [statistics, mutantsClasses_wR, listOutliers, growthMTs, processedData] = constrainUnmeasuredSpecies('mpi', 1, 1)

#### Output:
- growthMTs	 -  optimal growth values computed for Arabidopsis T-DNA lines.

- processedData  -  struct with the fields:
  - processedProfiles - fold change computed for lipid relative abundances (m = RA-T-DNA/RA-WT).
  - referenceFluxes	  - incoming reactions and fluxes for each of the lipid species measured.
  - modelIrrev - irreversible version of the model built for Arabidopsis wild-type Col-0 plants.
  - flagRelaxation - flag that indicates if simulations are performed relaxing bounds with standard deviation (1 - bounds are relaxed).
  - modelsMT_w - models built for Arabidopsis T-DNA lines constrained with flux sums and lipid abundances.

- statistics  -  struct with the fields:
  - meanMTs_wR - mean flux values for reactions underlying the mutated locus.
  - oneSampleTtest_wR - results of one-sample t-Test performed for mean flux values of T-DNA lines.
  - twoSampleTtest_wR - results of two-sample t-Test for mean flux values of wild-type Col-0 plants and T-DNA lines.
  - tsTtestUnique_wR - results of two-sample t-Test for reactions catalyzed by unique enzyme product.
  - tsTtestIsoenzymes_wR - results of two-sample t-Test for reactions catalyzed by isoenzymes.

- mutantsClasses_w  -  summary of results of mean fluxes computed for Arabidopsis wild-type Col-0 plants and T-DNA lines.

- listOutliers  -  list of T-DNA lines classified as outliers due to their disagreement among predicted fluxes and published phenotypic information.

#### Description files generated:
The results generated during the simulations are stored as .xlsx files in the folder 'Ttest' located in the path '..\SamplingResults\Ttest':
- Results_statistic_tests_Outl_(lusk/mpi) - results of the one- and two-sample t-Tests performed for a selected group of Arabidopsis T-DNA lines whose models were constrained at the level of lipid pools.
- Results_classification_mutants_Outl_(lusk/mpi) - classification of the Arabidopsis T-DNA lines according to differences found in the mean of fluxes compared to Col-0 plants.


### Perform differential flux analysis for selected Arabidopsis T-DNA lines

#### Description:
The script 'differentialFluxAnalysis' is used to find out if the T-DNA lines for which absence of phenotypic changes is predicted, may exhibit differences in other pathways indirectly or not related to the metabolic steps affected by the mutation.

#### Usage:
          [models, samples, resultsTtests, summaryFCs] = differentialFluxAnalysis

#### Output:
- models  -  models built for selected Arabidopsis T-DNA lines constrained with flux sums and lipid abundances.

- samples  -  flux samples obtained via ACHR for selected Arabidopsis T-DNA lines.

- resultsTtests	 -  results of t-Test performed to identify reactions that exhibit differential fluxes.

- summaryFCs  -  summary of reactions of terpenoid pathway predicted to be up-regulated in _GPAT1_ (_AT1G06520_) T-DNA line.

#### Description files generated:
The results generated during the simulations are stored as .xlsx files in the folder 'Ttest' located in the path '..\SamplingResults\Ttest':
- barplots_differential_fluxes - summary of results of differential flux analysis for a group of four selected Arabidopsis T-DNA lines. The reactions that exhibit differential fluxes are classified into different categories. The data is formatted for use in creating Figure 4 and Supplemental Figure 2.
- Validated_TDNAs_full_rescue - results of the validation of reactions that exhibit differential fluxes with transcripts data. The data is formatted for use in creating Figure 5.


### Construction of Figures 2A-B and Supplemental Figures 1A-B:
The scripts used below are located in the path '..\Figures_source_data\Figure_2_&_Supplemental_Figure_1'

1. Execute the MATLAB script 'analyzeLipidomes.m' by entering the command:

       [statsLipidsOutliers,statsLipidsConcordants,statsLipidsUnknown] = analyzeLipidomes

3. Execute the MATLAB script 'formatData4heatmaps.m' entering the command:

       [prunedFC_mpi,prunedFC_lusk,prunedScaffold_mpi,prunedScaffold_lusk] = formatData4heatmaps

5. Run the R script 'HeatMapGeneration.R'


### Construction of Figure 4 and Supplemental Figure 2:
1. Execute the MATLAB script 'differentialFluxAnalysis.m' located in the path '..\SimulateGrowthMTs', by entering the command:

        [models, samples, resultsTtests, summaryFCs, summaryGeneExpression] = differentialFluxAnalysis
   The script creates the file 'barplots_differential_fluxes.xlsx', that is saved in the path '..\SamplingResults\Ttest'.

2. Run the R script 'BarPlots_differential_fluxes.R' located in the path '..\Figures_source_data\Figure_4_&_Supplemental_Figure_2'


### Construction of Figure 5:
1. If the Matlab script 'differentialFluxAnalysis.m' has not been executed, follow the instructions described above.
   The script generates as output the file 'Validated_TDNAs_full_rescue.xlsx' located in the path: '..\SamplingResults\Ttest'.

3. Run the R script 'Horizontal_barPlot.R' located in the path '..\Figures_source_data\Figure_5'


### Construction of Supplemental Figure 3:
The principal component analysis was performed using the webserver MetaboAnalyst (Chong et al., 2019; https://www.metaboanalyst.ca/). The data stored in the files 'lusk_MTs.csv' and 'mpi_MTs.csv' were uploaded to the web server. Data were auto scaled by mean-centering and dividing by the standard deviation of each variable. The corresponding figures and analysis reports were exported from the web server and can be found in the path: '..\Figures_source_data\Supplemental_figure_3'.


## Contributors

Correa S. (Cordoba@mpimp-golm.mpg.de, sandra.correa.cordoba.1@uni-potsdam.de), Bioinformatics Group, University of Potsdam - Germany.

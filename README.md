# PBS
Authors: Xiting Yan, Anqi Liang, Hongyu Zhao and Geoffrey Chupp

These are the instructions to run the R version of the PBS program (PBS-R.ZIP). This version is intended for more computational experienced biologists, bioinformaticians or computational biologists. Another version for users with less experiences is still under development.

The PBS-R program described here reflects the version of the methodology described and used in the Yan et al 2016 manuscript. For details about the method and the content of the output please see the manuscript.


## Overview
Distance based unsupervised clustering of the gene expression data is commonly used to identify heterogeneity in samples. Due to the high noise level in the gene expression data and correlation between genes, traditional distances such as Euclidean distance may not be appropriate at representing the true biological differences between samples. Pre-defined biological pathways have been associated with disease phenotypes in a way that different genes in the same pathway can be perturbed in different subjects that have similar clinical outcomes. Therefore assessing the biological differences between samples based on pathways instead of single genes may enhance the robustness and accuracy of the clustering results. 

PBS is a software package that assesses the biological differences between samples using gene expression data by assuming that ontologically defined biological pathways in biologically similar samples have similar behavior. This package includes three parts: simulation, pathways and pbs_distance. The simulation codes generate the simulated data as described Yan et al. (2015). The pathways part provides different types of pre-defined pathways prepared in the required format for the pbs_distane codes. The pbs_distancec part calculate the pathway based distance score between subjects using the gene expression data and pre-defined pathways. 

## Instructions
You need to install R release 3.0 or later.

1. Copy the PBS-R.ZIP file to your computer. 

2. Unzip the file PBS-R.ZIP using the option to create subdirectories.
  This should create the following files and subdirectories:

  * PBS program and functions in R (all the PBS code is conatined there):

    ```
    PBS-R/PBS.1.0.R
    PBS-R/simulation.R
    ```

  * Directory with example datasets, gct file:
    ```
    PBS-R/Datasets/        
                    asthma.gct
                    simul1_independent.gct
                    simul1_depedent.gct
    ```
  * Directory with gene set databases, gmt files:
    ```
    PBS-R/GeneSetDatabases/
                            C1.gmt
                            C2.gmt
                            C3.gmt
                            C4.gmt
                            Lung_Boston_poor_outcome.gmt
                            Lung_Michigan_poor_outcome.g
    ```
  * Directories with results of running the examples described in the paper:
    ```
    PBS-R/asthma_C2/
                      asthma
    ```

  * One page R scripts to run the examples described in the paper:
    ```
    GSEA/GSEA-P-R/
                  Run.Gender_C1.R
                  Run.Gender_C2.R
                  Run.Leukemia_C1.R
                  Run.Lung_Boston_C2.R
                  Run.Lung_Stanford_C2.R
                  Run.Lung_Michigan_C2.R
                  Run.Lung_Boston_outcome.R
                  Run.Lung_Michigan_outcome.R
                  Run.P53_C2.R
    ```

To run, for example, the simulated dataset 1 with independent setting with the C2 gene set database go to the file PBS-R/Run.Simul1_C2.R and change the file pathnames to reflect the location of the GSEA directory in your machine. For example if you expanded the ZIP file under your directory "C:/my_directory" you need to change the line: 
```
GSEA.program.location <- "d:/CGP2005/GSEA/GSEA-P-R/GSEA.1.0.R"  
```
To:
```
GSEA.program.location <- "c:my_directory/GSEA/GSEA-P-R/GSEA.1.0.R"
```
 And the same change to each pathname in that file: you need to replace "d:/CGP2005" with "C"/my_directory".

 You may also want to change the line:

doc.string            = "Leukemia_C1",

To:

doc.string            = "my_run_of_Leukemia_C1",

or any other prefix label you want to give your results. This way you won't overwrite the original results that come in those directories and can use them for comparison with the results of you own run. 

After the pathnames have been changed to reflect the location of the directories in your machine to run GSEA program just open the R GUI and paste the content of the Run.<example>.R files on it.  Fro example to run the Leukemia vs. C1 example use the contents of the file "Run.Leukemia_C1.R" The program is self-contained and should run and produce the results under the directory "C:my_directory/GSEA/GSEA-P-R/Leukemia_C1". These files are set up with the parameters used in the examples of the paper (e.g. to produce detailed results for the significant and top 20 gene sets). You may want to start using these parameters and change them only when needed and when you get mnore experience with the program. For details of what are the effect of changing some of the parameters see the Supporting Information document.

If you want to run a completely new dataset the easiest way is:

i) Create a new directory: e.g. GSEA/GSEA-P-R/my_dataset, where you can store the inputs and outputs of running GSEA on those files. 
ii) Convert manually your files to *.gct (expression dataset) and *.cls (phenotype labels)
iii) Use Run.Leukemia_C1.R as a template to make a new script to run your data.
iv) Change the relevant pathnames to point to your input files in directory my_dataset. Change the doc.string to an approprote prefix name for your files.
v) Cut and paste the contents of this new script file in the R GUI to run it. The results will be stored in my_directory.

The GSEA-P-R program reads input files in *.gct, *.cls and *.gmt formats. As you can see from the examples's files these are simple tab separated ASCII files. If your datasets are not in this format you can use a text editor to convert them. If you start with a tab separated ASCII file tipically the conversion would consist in  modifying the header lines on top of the file.

If you have questions or problems running or using the program please  send them to gsea@broad.mit.edu. Also lets us know if you find GSEA a useful tool in your work.





##References
Yan X, Chu JH, Gomez J, Koenigs M, Holm C, He X, Perez MF, Zhao H, Mane S, Martinez FD, Ober C, Nicolae DL, Barnes KC, London SJ, Gilliland F, Weiss ST, Raby BA, Cohn L, Chupp GL. 2015. [Noninvasive analysis of the sputum transcriptome discriminates clinical phenotypes of asthma](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4451618/). Am J Respir Crit Care Med. 191(10):1116-25. 

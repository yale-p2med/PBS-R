my.dir<-"/Users/yanxiting/Documents/Research/Xiting/Pathway_Clustering/PBS-R"

setwd(my.dir)
source("PBS.1.0.R")
data.filepath<-"./Datasets/independent_simulated_data1_datamatrix.gct"
msigdb.filepath<-"./GeneSetDatabases/c2.cp.v3.0.symbols_mapped_to_HuGene_1_0_st.chip.gm.gmt"
pathway.name.prefix<-"KEGG_"
output.dir<-"./results_simul_independent_c2"
result.name.prefix<-"independent_simul1"
pbs.cal(data.filepath,msigdb.filepath,pathway.name.prefix,output.dir,result.name.prefix)

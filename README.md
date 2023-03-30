# mouse TEspeX  

This the transposon analysis developped here:  

https://github.com/fansalon/TEspeX  


1- The mouse TE counting is started with the mouse_te_v0.01.pbs command after putting the different fastq files in a folder  
2- The RNA-Seq_mouse pipeline from the Genomics core is usually ran previous to the TE counting, allowing the merging of the two datasets   
3- The STAR results from the RNA-Seq_mouse pipeline is placed inside the folder  
4- The differential expression of genes and TEs is made by DESeq2 using the deseq2_wald_two_groups.R  



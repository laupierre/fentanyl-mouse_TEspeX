# mouse TEspeX  

This the TEspeX transposon analysis developped here:  

https://github.com/fansalon/TEspeX  


1- The mouse TE counting is started with the mouse_te_v0.01.pbs command after putting the different fastq files in a folder  
2- The RNA-Seq_mouse pipeline is usually ran previous to the TE counting, allowing the merging of the two datasets   
3- The STAR results from the RNA-Seq_mouse pipeline is placed inside the same folder as in 1)  
4- In this example, we used the fentanyl mouse RNA=Seq experiment. The differential expression of genes and TEs is made by DESeq2 using the deseq2_wald_two_groups.R.  

The results of the 3 comparisons is obtained:
hippocampus_deseq2_FentanylvsControl_differential_expression.xlsx (F vs C)
hippocampus_deseq2_WithdrawalvsControl_differential_expression.xlsx (W vs C)
hippocampus_deseq2_WithdrawalvsFentanyl_differential_expression.xlsx (W vs F)



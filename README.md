# 16s_absolute_quantification
Code to produce absolute abundance of microbial communities using 16s amplicon sequencing. 
1. Create environment for subsequent steps, cutadapt to remove primer
2. DADA2 to denoise, create ASV table, and using silva database to assign taxonomy
3. Using the ASV table and taxomony table generated in 2 to identify spike-in, rescaling the rest of samples

Source documentations:
https://benjjneb.github.io/dada2/tutorial.html
https://astrobiomike.github.io/amplicon/dada2_workflow_ex

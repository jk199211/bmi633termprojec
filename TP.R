####################################################

### Simple implementation of PRS calculation    #### 

### 1. first using old static p - value threshold ##

### 2. then using improved bonferroni method #######

### Members: James Kambere and Purushottam Panta ###

####################################################

### installing required libraries

    install.packages('R.utils')
    install.packages("BEDMatrix")
    install.packages("LDlinkR")
    install.packages("tidyverse")
    install.packages("data.table")
    
### first set we will load summary statistics data

    setwd("C:/Users/AGILE SYSTEMS/OneDrive - LUANAR/UKY/BMI633/termproject")
    library(data.table)
    
### Read in file
    
    ## attached library R.utils for reading .gz files
    
    library(R.utils)
    dat <- fread("Height.gwas.txt.gz")
    nrow(dat)

### Filter out SNPs with low MAF and INFO
  
    result <- dat[INFO > 0.8 & MAF > 0.01]
    
    ## view results
    
    result
    nrow(result)

### Output the gz file
     
    fwrite(result, "Height.gz", sep="\t")

### QC
    
    ## attaching tidyverse for sub setting data frames
    
    library(tidyverse)
    hieght_dt <-fread("Height.gz")

    ### Removing duplicates SNPs
  
    no_duplicates <- hieght_dt[!duplicated(SNP)]
    nrow(no_duplicates)
    fwrite(no_duplicates,"Hieght.QC.gz",sep="\t")

### updating beta in the file

    dat <- fread("Hieght.QC.gz")
    fwrite(dat[,BETA:=log(OR)], "Height.QC.Transformed", sep="\t")
    
### getting phenotype data from .fam, .bed and .bim
    
    ## creating a list of individuals from given .bed file
    
    ## attaching BEDMAtrix package
    
    library(BEDMatrix)
    
    ## this path can be changed to suit the location of phenotype files
    
    path <- system.file("extdata", "EUR.bed",package = "BEDMatrix")
    
    # Create a BEDMatrix object the example .bed file
    
    m1 <- BEDMatrix(path)
    mysnps = data.frame(m1[1:2, 1:529493])
    
    ## getting a list of snps
    
    get_name = colnames(mysnps)
    
    ## converting snps to data frame
    
    data_frame_ofnames = data.frame(get_name)
    with_exploded_names = data.frame(do.call("rbind", strsplit(as.character(data_frame_ofnames$get_name), "_", fixed = TRUE)))
    
    ## transposing the bed
    
    new_data = as.matrix(mysnps)
    totranspose <- as.matrix(new_data)
    tansposed <- data.frame(t(as.matrix(totranspose)))
    
    ## view results, the effect allele is concatinated to the SNPs
    
    tansposed
    
    ## sub setting
    
    ## here we are selecting person 1 and 2 only 
    
    person1 = data.frame(tansposed[,1])
    person2 = data.frame(tansposed[,2])
    
    datafile <- cbind(with_exploded_names,person1,person2)
    
    ## view the list
    
    datafile
    
    ## Remove names
    
    names(datafile) <- NULL
    
    ## add real names
    
    names(datafile) <- c("SNP","A1","PER001","PER002")
    
    ## view file
    
    datafile
    
    ## removing N/As (another QC)
    
    ## first convert to data table
    
    new_file = data.table(datafile)
    
    ## view it
    
    new_file
    
    ## remove n/a
    
    new_file_withoutna_1 = new_file[!is.na(new_file$PER001),]
    new_file_withoutna_2 = new_file[!is.na(new_file$PER001),]
    
    ## complete phenotype file
    
    new_file_withoutna_2
    
### Linkage disequilibrium for r - squared
    
    qc_transformed  = fread("Height.QC.Transformed")
    new_file = subset(qc_transformed, select = c(SNP))
    new_data = as.data.frame(new_file)
    
    ## getting SPNs into a vector 
    
    sps_list = NULL
    
    ## test with 5000 snps
    
    ## might take more time it uses an online tool. sometime tocken generation may be required
    
    for(i in 1:5000){
      
      spns = new_data[i,1]
      print(spns)
      
      sps_list = rbind(sps_list, spns)
      
    }
    
    final = c(sps_list)
    
    ## view list
    
    final
    
    ## installing required packages 
    
    ## attaching a clumping tool LDlinkR
    
    library(LDlinkR)
    
    ## performing clumping
    
    ## here we are using correlation coefficient of 0.2
    
    remaining_alleles = SNPclip(final, "YRI", "0.2", "0.01", token = "37cf85ff085f")
    
    remaining_alleles_2 = data.frame(remaining_alleles)
    remaining_final = remaining_alleles_2[which(remaining_alleles_2$Details == "Variant kept."),]
    remaining_final
    nrow(remaining_final)
    
    ## drop names 
  
    final_trimmed_spns = data.table(subset(remaining_final, select = c(RS_Number)))
    final_trimmed_spns
    nrow(final_trimmed_spns)
    
    ## add column name
    
    names(final_trimmed_spns) = c("SNP")
    final_trimmed_spns
    
    ## merging with the other files
    
    require(dplyr)
    
    all_snps = as.data.frame(qc_transformed)
    prunned_spns = data.table(final_trimmed_spns)
    
    snps.compare_to = (all_snps$SNP %in% prunned_spns$SNP )
    final_crumped_list = all_snps[snps.compare_to,]
    final_crumped_list
    
    
### p - value threshold

    to_order <- fread("Height.QC.Transformed")
    ordered <- to_order[order(P)]
    n = nrow(ordered)
    
    ## old way of threshold for comparison sake
    
    pval = 0.05/n
    pval
    
    pruned_data_set = data.table(final_crumped_list)
    subset_pruned3 = pruned_data_set[P <= pval]
    nnew = nrow(subset_pruned3)
    subset_prunned_final = data.table(subset_pruned3)
    
    subset_prunned_final
    
### old style of calculating PRS 
   
    ## concatenating with the phenotype file
    
    prs_snp_file = merge(subset_prunned_final,new_file_withoutna_2,by="SNP")
    
    ## calculating the weight of each snps
    
    ## sample data, rows are SNPs, columns are individuals
    
    id_names <- paste0("PER00",1:2)
    
    ## the names of the new columns that have been multiplied by beta
    
    new_id_names <- paste0("new_",id_names)
    
    prs_snp_file[, (new_id_names) := .SD * BETA, .SDcols = id_names]
    prs_snp_file
    
    ## getting the sum of scores
    
    prs_old_way = prs_snp_file[, colSums(.SD), .SDcols = new_id_names]
    prs_old_way
    
### new way with improved bonferronis to get more SNPs as possible
    
    data_snp = data.table(final_crumped_list)
    ordered_pvals <- data_snp [order(data_snp $P)]
    crumped_snps_list = data.frame(ordered_pvals[,c(8)])
    m = nrow(crumped_snps_list)
    m
    
    num.rows = nrow(ordered_ps)
    
    df = NULL
    
    for(i in 1:m){
      
      pvalue = crumped_snps_list[i,1]
      imp.bon = i*0.05/num.rows
      
      if(pvalue <= imp.bon){
        
        print(pvalue)
        df = rbind(df, data.frame(pvalue))
        
      }
      
    }
    
    ## sub setting crumped SPNs
    
    require(dplyr)
    
    fullmodel = as.data.frame(final_crumped_list)
    crumpedmodel = as.data.frame(df)
    
    snps.compare = (fullmodel$P %in% crumpedmodel$pvalue)
    significant_snps = fullmodel[snps.compare,]
    
    final_snp = data.table(significant_snps)
    
    ## check to confirm number
    
    n = nrow(final_snp)
    n
    
### calculation of PRS
    
    ## concatenating with the phenotype file
    
    prs_snp_file2 = merge(final_snp,new_file_withoutna_2,by="SNP")
    
    prs_snp_file2
    
    ## calculating the weight of each snps
    
    ## sample data, rows are SNPs, columns are individuals
    
    id_names <- paste0("PER00",1:2)
    
    ## the names of the new columns that have been multiplied by beta
    
    new_id_names <- paste0("new_",id_names)
    
    prs_snp_file2[, (new_id_names) := .SD * BETA, .SDcols = id_names]
    prs_snp_file2
    
    ## getting the sum of scores
    
    prs_new_way = prs_snp_file2[, colSums(.SD), .SDcols = new_id_names]
    prs_new_way
    
### graphical comparison
    
    library(ggplot2)
    
    datf <- data.frame(method=c("B-OLD", "B-NEW"),
                     SNPs=c(19, 46))
    ggplot(data=datf, aes(x=method, y=SNPs)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_minimal()
    
### end of programe ####
    
    
  
    
    


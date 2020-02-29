#A Multidimensional Systems Biology Analysis of Cellular Senescence in Ageing and Disease
#http://genomics.senescence.info/cells/

#Script to scramble GTEx accessions to find how significant CellAge DEGs with age differentially
#expressed in different tissues are

library(dplyr)

#Cutoff for significant log2C
log2_threshold=log2(1.5)
pval_filter=0.05

#Number of simulations
z=0
# zmax=10000 #Set max number (10,000)
zmax=10000

#Set working directory to wherever all GTEx data is
# setwd('/home/nis/roberto/Scripts/Scrambled_expression_simulations/all_GTEx_expression_results')
# add cellage build 1 here
# cellage=read.csv('/home/nis/roberto/Scripts/Scrambled_expression_simulations/cellage_pedro_original.csv', header=TRUE, stringsAsFactors = FALSE)%>%dplyr::select(gene_name, senescence_effect)%>%dplyr::filter(senescence_effect!='Unclear')
# add signatures of CS from https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13041 here
# CS_overexpressed=read.csv('/home/nis/roberto/Scripts/Scrambled_expression_simulations/signatures_of_CS/Kasit overexpressed CS signatures.csv', header=TRUE)%>%dplyr::select(external_gene_name)
# CS_underexpressed=read.csv('/home/nis/roberto/Scripts/Scrambled_expression_simulations/signatures_of_CS/Kasit underexpressed CS signatures.csv', header=TRUE)%>%dplyr::select(external_gene_name)
cs=rbind(CS_overexpressed, CS_underexpressed)
cs=unique(cs)
colnames(cs)='gene_name'
#Make dataframe of file names
files=as.data.frame(list.files())#

#Make dataframe of tissues
tissues=c()#
for (i in 1:nrow(files)){#
  x=unlist(strsplit(as.character(files$`list.files()`[i]), '[_]'))[2]#
  tissues=rbind(tissues, x)#
}#

tissues=as.data.frame(tissues)#

#Read in files, calculate LogFC50, divide tables into over and underexpressed
for (i in 1:nrow(tissues)){
  #Add underscores to tissues with spaces
  x=unlist(strsplit(as.character(tissues$V1[i]), ' '))#
  if (length(x)>1){#
    x=paste(x[1], x[2], sep='_')#
  }else{#
    x=as.character(tissues$V1[i])#
  }#
  
  #Read in files, add in Tissue Name
  y=(read.csv(paste('GTEx', tissues$V1[i], 'full_result_gene_name.csv', sep='_'), stringsAsFactors = FALSE))%>%dplyr::mutate(tissue=tissues$V1[i])#
  
  #Remove duplicate rows from when we convert Ensembl to gene symbol
  y = y[!duplicated(y$external_gene_name),]#
  
  #Move gene symbol to front of table, remove ensembl and entrez ID
  y=y%>%dplyr::select(external_gene_name, everything(), -ensembl_id, -entrezgene)#
  
  #Assign each table to a different variable
  assign(paste(x, 'DEGs', sep='_'), y) #
}

scrambled_overlaps_overexpressed=c()
scrambled_overlaps_underexpressed=c()

#CS DEGs
scrambled_overlaps_CS_overexpressed=c()
scrambled_overlaps_CS_underexpressed=c()

while(z+1<=zmax){
  #set.seed
  set.seed(z+2)
  scrambled_tissue_overexpressed2=c()
  scrambled_tissue_underexpressed2=c()
  saveme=c()
  for (i in 1:nrow(tissues)){
    x=unlist(strsplit(as.character(tissues$V1[i]), ' ')) #Change spaces to underscores in tissue names
    if (length(x)>1){ 
      x=paste(x[1], x[2], sep='_')
    }else{
      x=as.character(tissues$V1[i])
    }
    
    #scramble
    scrambled_tissue=c()
    
    scrambled_tissue=cbind(sample((get(paste(x, 'DEGs', sep='_')))$external_gene_name),#Scramble gene names
                           get(paste(x, 'DEGs', sep='_'))[,2:ncol(get(paste(x, 'DEGs', sep='_')))])
    saveme=rbind(saveme, scrambled_tissue)
    #overexpressed
    scrambled_tissue_overexpressed2=rbind(scrambled_tissue_overexpressed2, scrambled_tissue%>%dplyr::filter(logFC_50>0)) #bind tissue[i] to dataframe containing all scrambled tissue data
    
    #underexpressed
    scrambled_tissue_underexpressed2=rbind(scrambled_tissue_underexpressed2, scrambled_tissue%>%dplyr::filter(logFC_50<0)) #bind tissue[i] to dataframe containing all scrambled tissue data
  }
  #save scrambled results to server
  colnames(saveme)[1]='scrambled_gene'
  #Write individual simulations here
  # write.table(x=saveme, file=paste('/home/nis/roberto/Scripts/Scrambled_expression_simulations/scrambled_data/simulation_', z+1, '.txt', sep=''),row.names=FALSE,quote=FALSE,sep='\t')
  
  #Filter for p-val and log2FC
  scrambled_tissue_overexpressed2=scrambled_tissue_overexpressed2%>%dplyr::filter(adj.P.Val<pval_filter, logFC_50>log2_threshold)
  
  #Absolute value for underexpressed
  scrambled_tissue_underexpressed2=scrambled_tissue_underexpressed2%>%dplyr::filter(adj.P.Val<pval_filter, abs(logFC_50)>log2_threshold)
  
  #Change the name of the first column to gene_symbol
  colnames(scrambled_tissue_overexpressed2)[1]='gene_symbol' 
  colnames(scrambled_tissue_underexpressed2)[1]='gene_symbol' 
  
  #Bind CellAge overlaps to overexpressed table
  x=as.data.frame(table(scrambled_tissue_overexpressed2$gene_symbol))%>%dplyr::arrange(-Freq)%>%
    dplyr::filter(Var1%in%cellage$gene_name)#Filter out genes that are not significantly expressed
  x=as.data.frame(table(x$Freq))%>%dplyr::mutate(test=z+1)
  scrambled_overlaps_overexpressed=rbind(scrambled_overlaps_overexpressed, x)
  
  #Bind CS overlaps to scrambled_overlaps_CS_overexpressed
  x=as.data.frame(table(scrambled_tissue_overexpressed2$gene_symbol))%>%dplyr::arrange(-Freq)%>%
    dplyr::filter(Var1%in%cs$gene_name)#Filter out genes that are not significantly expressed
  x=as.data.frame(table(x$Freq))%>%dplyr::mutate(test=z+1)
  scrambled_overlaps_CS_overexpressed=rbind(scrambled_overlaps_CS_overexpressed, x)
  
  #Bind CellAge overlaps to underexpressed table
  x=as.data.frame(table(scrambled_tissue_underexpressed2$gene_symbol))%>%dplyr::arrange(-Freq)%>%
    dplyr::filter(Var1%in%cellage$gene_name)#Filter out genes that are not significantly expressed
  x=as.data.frame(table(x$Freq))%>%dplyr::mutate(test=z+1)
  scrambled_overlaps_underexpressed=rbind(scrambled_overlaps_underexpressed, x)
  
  #Bind CS overlaps to scrambled_overlaps_CS_underexpressed
  x=as.data.frame(table(scrambled_tissue_underexpressed2$gene_symbol))%>%dplyr::arrange(-Freq)%>%
    dplyr::filter(Var1%in%cs$gene_name)#Filter out genes that are not significantly expressed
  x=as.data.frame(table(x$Freq))%>%dplyr::mutate(test=z+1)
  scrambled_overlaps_CS_underexpressed=rbind(scrambled_overlaps_CS_underexpressed, x)
  
  #loop
  z=z+1
  print(z)
}

colnames(scrambled_overlaps_overexpressed)=c('tissues overexpressing', 'number CellAge', 'test')
colnames(scrambled_overlaps_underexpressed)=c('tissues underexpressing', 'number CellAge', 'test')
colnames(scrambled_overlaps_CS_overexpressed)=c('tissues overexpressing', 'number CS', 'test')
colnames(scrambled_overlaps_CS_underexpressed)=c('tissues underexpressing', 'number CS', 'test')

#Get unique number of genes, find how many cellage and CS genes are actually in any of the expression data
PC_genes=unique(saveme%>%dplyr::select(scrambled_gene))

#get cellage number
cellage_number=nrow(PC_genes%>%dplyr::filter(scrambled_gene%in%cellage$gene_name))

#get CS number
CS_number=nrow(PC_genes%>%dplyr::filter(scrambled_gene%in%cs$gene_name))

# % simulations overlapping
##overexpressed CellAge
scrambled_overlaps_overexpressed$`tissues overexpressing`=as.numeric(as.character(scrambled_overlaps_overexpressed$`tissues overexpressing`))

#For every tissue overlap number, find % of CellAge genes differentially expressed in multiple tissues
overexpressed_percent=c()
##for minimum number of tissues overlapping to maximum number of tissues overlapping
for (i in min(scrambled_overlaps_overexpressed$`tissues overexpressing`):max(scrambled_overlaps_overexpressed$`tissues overexpressing`)){
  ##overlap_sum is the sum of cellage genes in each tissue overlap number 
  overlap_sum=sum((scrambled_overlaps_overexpressed%>%dplyr::filter(scrambled_overlaps_overexpressed$`tissues overexpressing`==i))$'number CellAge')
  overlap_perc=overlap_sum/(cellage_number*zmax)*100
  overexpressed_percent=rbind(overexpressed_percent, c(i, overlap_perc))
}

colnames(overexpressed_percent)=c('Number_tissues', 'perc_cellage')

# ##overexpressed CS
scrambled_overlaps_CS_overexpressed$`tissues overexpressing`=as.numeric(as.character(scrambled_overlaps_CS_overexpressed$`tissues overexpressing`))

#For every tissue overlap number, find % of cs genes differentially expressed in multiple tissues
overexpressed_CS_percent=c()
##for minimum number of tissues overlapping to maximum number of tissues overlapping
for (i in min(scrambled_overlaps_CS_overexpressed$`tissues overexpressing`):max(scrambled_overlaps_CS_overexpressed$`tissues overexpressing`)){
  ##overlap_sum is the sum of cs genes in each tissue overlap number 
  overlap_sum=sum((scrambled_overlaps_CS_overexpressed%>%dplyr::filter(scrambled_overlaps_CS_overexpressed$`tissues overexpressing`==i))$'number CS')
  overlap_perc=overlap_sum/(CS_number*zmax)*100
  overexpressed_CS_percent=rbind(overexpressed_CS_percent, c(i, overlap_perc))
}

colnames(overexpressed_CS_percent)=c('Number_tissues', 'perc_CS')

##underexpressed CellAge
scrambled_overlaps_underexpressed$`tissues underexpressing`=as.numeric(as.character(scrambled_overlaps_underexpressed$`tissues underexpressing`))

#For every tissue overlap number, find % of CellAge genes differentially expressed in multiple tissues
underexpressed_percent=c()
##for minimum number of tissues overlapping to maximum number of tissues overlapping
for (i in min(scrambled_overlaps_underexpressed$`tissues underexpressing`):max(scrambled_overlaps_underexpressed$`tissues underexpressing`)){
  ##overlap_sum is the sum of cellage genes in each tissue overlap number 
  overlap_sum=sum((scrambled_overlaps_underexpressed%>%dplyr::filter(scrambled_overlaps_underexpressed$`tissues underexpressing`==i))$'number CellAge')
  overlap_perc=overlap_sum/(cellage_number*zmax)*100
  underexpressed_percent=rbind(underexpressed_percent, c(i, overlap_perc))
}

colnames(underexpressed_percent)=c('Number_tissues', 'perc_cellage')

# ##underexpressed CS
scrambled_overlaps_CS_underexpressed$`tissues underexpressing`=as.numeric(as.character(scrambled_overlaps_CS_underexpressed$`tissues underexpressing`))

#For every tissue overlap number, find % of cs genes differentially expressed in multiple tissues
underexpressed_CS_percent=c()
##for minimum number of tissues overlapping to maximum number of tissues overlapping
for (i in min(scrambled_overlaps_CS_underexpressed$`tissues underexpressing`):max(scrambled_overlaps_CS_underexpressed$`tissues underexpressing`)){
  ##overlap_sum is the sum of cs genes in each tissue overlap number 
  overlap_sum=sum((scrambled_overlaps_CS_underexpressed%>%dplyr::filter(scrambled_overlaps_CS_underexpressed$`tissues underexpressing`==i))$'number CS')
  overlap_perc=overlap_sum/(CS_number*zmax)*100
  underexpressed_CS_percent=rbind(underexpressed_CS_percent, c(i, overlap_perc))
}

colnames(underexpressed_CS_percent)=c('Number_tissues', 'perc_CS')

#convert percent into p-value and adjusted p-value
##CellAge overexpressed simulation
overexpressed_percent=as.data.frame(overexpressed_percent)
overexpressed_percent=overexpressed_percent%>%dplyr::mutate(pval=overexpressed_percent$'perc_cellage'/100)
overexpressed_percent$adj_pval=p.adjust(overexpressed_percent$pval, 'BH')
overexpressed_percent$DEG='Overexpressed'

#CS overexpressed simulation
overexpressed_CS_percent=as.data.frame(overexpressed_CS_percent)
overexpressed_CS_percent=overexpressed_CS_percent%>%dplyr::mutate(pval=overexpressed_CS_percent$'perc_CS'/100)
overexpressed_CS_percent$adj_pval=p.adjust(overexpressed_CS_percent$pval, 'BH')
overexpressed_CS_percent$DEG='Overexpressed'

#convert percent into p-value and adjusted p-value
##CellAge underexpressed simulation
underexpressed_percent=as.data.frame(underexpressed_percent)
underexpressed_percent=underexpressed_percent%>%dplyr::mutate(pval=underexpressed_percent$'perc_cellage'/100)
underexpressed_percent$adj_pval=p.adjust(underexpressed_percent$pval, 'BH')
underexpressed_percent$DEG='Underexpressed'

#CS underexpressed simulation
underexpressed_CS_percent=as.data.frame(underexpressed_CS_percent)
underexpressed_CS_percent=underexpressed_CS_percent%>%dplyr::mutate(pval=underexpressed_CS_percent$'perc_CS'/100)
underexpressed_CS_percent$adj_pval=p.adjust(underexpressed_CS_percent$pval, 'BH')
underexpressed_CS_percent$DEG='Underexpressed'

#Print results
print('CellAge')
overexpressed_percent
underexpressed_percent

print('CS signatures')
overexpressed_CS_percent
underexpressed_CS_percent

print('Total number of simulations')
print(max(scrambled_overlaps_overexpressed$test))

#Save results CellAge
# write.table(overexpressed_percent,
#           '/home/nis/roberto/Scripts/Scrambled_expression_simulations/overexpressed_simulations_cellage.txt', row.names = FALSE,sep='\t')
# write.table(underexpressed_percent,
#           '/home/nis/roberto/Scripts/Scrambled_expression_simulations/underexpressed_simulations_cellage.txt', row.names=FALSE,sep='\t')

#Save results signatures of CS
# write.table(overexpressed_CS_percent,
#           '/home/nis/roberto/Scripts/Scrambled_expression_simulations/overexpressed_simulations_CS.txt', row.names = FALSE,sep='\t')
# write.table(underexpressed_CS_percent,
#           '/home/nis/roberto/Scripts/Scrambled_expression_simulations/underexpressed_simulations_CS.txt', row.names=FALSE,sep='\t')



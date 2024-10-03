mtd<-fread('metadata/proteomics_samples_batch_metadata.csv')

ratios<-fread('ref-data/proteomics_raw_data/abundances_ratio.csv')

abunds<-fread('ref-data/proteomics_raw_data/abundances.csv')

samples<-colnames(ratios)[2:ncol(ratios)]


#Ratio by batch
batch_infos<-fread('metadata/batch_infos.csv')

ggplot(batch_infos,aes(x=as.factor(batch),y=avg_ratio))+geom_violin()+geom_boxplot() 
ggplot(batch_infos,aes(x=as.factor(batch),y=log2(avg_abund)))+geom_violin()+geom_boxplot() 
ggplot(batch_infos,aes(x=as.factor(batch),y=med_ratio))+geom_boxplot(outlier.shape = NA)

ggplot(dcast(batch_infos,Accession~batch,value.var = 'Abundances (Grouped): Pool'))+
  geom_point(aes(x=b1,y=b2))

ratiosf<-ratios[pct.det>0.25] #5420
abundsf<-abunds[pct.det>0.25]

rna_mat<-data.frame(fread('/Alignment_final/featureCounts/featureCounts_clean.txt'),row.names = 'V1')
rna_cpms<-edgeR::cpm(rna_mat)
rna_mtd<-fread('/Alignment_final/featureCounts/RNAseqPhen_aug2021.csv')

prot_infos<-fread('metadata/proteins_infos.csv')

prot_infosf<-prot_infos[Accession%in%abundsf$Accession]
sum(rownames(rna_mat)%in%prot_infosf$`Ensembl Gene ID`) #5k
rna_cpmsf<-rna_cpms[rownames(rna_cpms)%in%prot_infosf$`Ensembl Gene ID`,]

prot_infosff<-prot_infosf[`Ensembl Gene ID`%in%rownames(rna_matf)]
abundsfa<-merge(abundsf,prot_infosff)


#Gene level matrix

library(biomaRt)
attrs<-GetBiomartAttrs(GetMartGenes())
attrs[str_detect(name,'uniprot')]
attrs[str_detect(name,'refse')]
attrs[str_detect(name,'entrez')]

trans_iso<-data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version','uniprot_isoform'),
                                             filters = c('uniprot_isoform'), 
                                             values = prot_infos$Accession ,
                                             mart = GetMartGenes()))
trans_uni<-data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version','uniprotswissprot'),
                                     filters = c('uniprotswissprot'), 
                                     values = prot_infos$Accession ,
                                     mart = GetMartGenes()))
trans_trb<-data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version','uniprotsptrembl'),
                                     filters = c('uniprotsptrembl'), 
                                     values = prot_infos$Accession ,
                                     mart = GetMartGenes()))

#Annotate Accession
prot_infos[,source:=ifelse(Accession%in%trans_uni$uniprotswissprot,'uniprotswissprot',ifelse(Accession%in%trans_iso$uniprot_isoform,'uniprot_isoform',ifelse(Accession%in%trans_trb$uniprotsptrembl,'uniprotsptrembl',NA)))]
prot_infos[is.na(source)]
trans_entrez<-data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version','entrezgene_id'),
                                     filters = c('entrezgene_id'), 
                                     values = unique(prot_infos[is.na(source)&`Ensembl Gene ID`=='']$`Entrez Gene ID`) ,
                                     mart = GetMartGenes()))


prot_infos[is.na(source)&str_detect(`Ensembl Gene ID`,'ENSG')] #225/347

trans_entrezf<-unique(trans_entrez,by='entrezgene_id')

#Add ensembl id from entrez
prot_infos[source=='entrezgene_id',ensembl_gene_id_ver:=trans_entrezf[prot_infos[source=='entrezgene_id']$`Entrez Gene ID`,on='entrezgene_id']$ensembl_gene_id_version]

trans_unif<-unique(trans_uni,by='uniprotswissprot')
prot_infos[source=='uniprotswissprot',ensembl_gene_id_ver:=trans_unif[prot_infos[source=='uniprotswissprot']$Accession,on='uniprotswissprot']$ensembl_gene_id_version]


trans_isof<-unique(trans_iso,by='uniprot_isoform')
prot_infos[source=='uniprot_isoform',ensembl_gene_id_ver:=trans_isof[prot_infos[source=='uniprot_isoform']$Accession,on='uniprot_isoform']$ensembl_gene_id_version]

trans_trbf<-unique(trans_trb,by='uniprotsptrembl')
prot_infos[source=='uniprotsptrembl',ensembl_gene_id_ver:=trans_trbf[prot_infos[source=='uniprotsptrembl']$Accession,on='uniprotsptrembl']$ensembl_gene_id_version]

prot_infos[!is.na(ensembl_gene_id_ver)]#8749/8871!


#Save gene annotations####
hgnc_symbols<-merge(unique(TransEnsemblVerstoSymbol(prot_infos$ensembl_gene_id_ver)[hgnc_symbol!='']),
                    unique(TransEnsembltoSymbol(prot_infos$ensembl_gene_id)[hgnc_symbol!='']),all=T)

prot_genes_infos<-merge(prot_infos,hgnc_symbols,all.x=T)
fwrite(prot_genes_infos,'metadata/proteins_genes_infos.csv')


#Aggregate matrix by ensgenome
abundsf_dt<-melt(abundsf,variable.name = 'sample_id',id.vars = c('Accession','pct.det'),value.name = 'abundance')
abundsf_dta<-merge(abundsf_dt,prot_genes_infosf)
abundsf_dta[,aggr.abund:=sum(abundance,na.rm = T),by=.(ensembl_gene_id_ver,sample_id)]
abunds_gl<-unique(abundsf_dta,by=c('ensembl_gene_id_ver','sample_id'))
ab_gene_mat<-dcast(abunds_gl,ensembl_gene_id_ver~sample_id,value.var = 'aggr.abund')
fwrite(ab_gene_mat,'ref-data/proteomics_raw_data/abundances_gene_level_aggregated_filtered.csv')

batch_infosfa<-merge(batch_infos,unique(abundsf_dta[,.(Accession,ensembl_gene_id_ver)]))
batch_infosfa[,aggr.pool.abund:=sum(`Abundances (Grouped): Pool`,na.rm = T),by=.(ensembl_gene_id_ver,batch)]

#Calculate new ratio
abunds_gl_a<-merge(abunds_gl,mtd[,.(sample_id,batch)],by='sample_id')
abunds_gl_ap<-merge(abunds_gl_a[,-'Accession'],unique(batch_infosfa,by=c('ensembl_gene_id_ver','batch')),by=c('ensembl_gene_id_ver','batch'))

ggplot(abunds_gl_ap)+geom_density(aes(x=aggr.ratio,col=batch))


ratio_gene_mat<-dcast(abunds_gl_ap,ensembl_gene_id_ver~sample_id,value.var = 'aggr.ratio')

fwrite(ratio_gene_mat,'ref-data/proteomics_raw_data/abundances_ratio_gene_level_aggregated_filtered.csv')


log2_gene_mat<-dcast(abunds_gl_ap,ensembl_gene_id_ver~sample_id,value.var = 'log2.ratio')

fwrite(log2_gene_mat,'ref-data/proteomics_raw_data/log2_abundances_ratio_gene_level_aggregated_filtered.csv')

fwrite(abunds_gl_ap,'ref-data/proteomics_raw_data/gene_level_proteomics_data_tidy.csv')
abunds_gl_ap<-fread('ref-data/proteomics_raw_data/gene_level_proteomics_data_tidy.csv')

# PCA
log2_gene_mat<-as.matrix(data.frame(log2_gene_mat,row.names = 'ensembl_gene_id_ver'))
log2_gene_mat<-log2_gene_mat[which(apply(log2_gene_mat, 1, var)!=0),]

pca<-prcomp(t(log2_gene_mat[rowSums(is.na(log2_gene_mat))==0,]),scale.=TRUE)

mtd[,sample_id_corr:=ifelse(str_detect(sample_id,'^[0-9]'),ps('X',sample_id),sample_id)]
pca_dt<-merge(data.table(pca$x,keep.rownames = 'sample_id_corr'),mtd)
pctpcs<-pctPC(pca)

pc_x<-'PC1'
pc_y<-'PC2'

p1<-ggplot(pca_dt)+geom_point(aes_string(x=pc_x,y=pc_y,col='cell_type'))+
  labs(x=paste0(pc_x,' (',round(pctpcs[pc_x]*100),'%)'),
       y=paste0(pc_y,' (',round(pctpcs[pc_y]*100),'%)'))+
  theme_bw()

p2<-ggplot(pca_dt)+geom_point(aes_string(x=pc_x,y=pc_y,col='genotype'))+
  labs(x=paste0(pc_x,' (',round(pctpcs[pc_x]*100),'%)'),
       y=paste0(pc_y,' (',round(pctpcs[pc_y]*100),'%)'))+
  theme_bw()


p3<-ggplot(pca_dt)+geom_point(aes_string(x=pc_x,y=pc_y,col='ipsc_type'))+
  labs(x=paste0(pc_x,' (',round(pctpcs[pc_x]*100),'%)'),
       y=paste0(pc_y,' (',round(pctpcs[pc_y]*100),'%)'))+
  theme_bw()



#DEA Analysis
library(prolfqua)
abunds_gl_ap<-fread('ref-data/proteomics_raw_data/gene_level_proteomics_data_tidy.csv')
sample_mtd<-fread('metadata/proteomics_samples_batch_metadata.csv') # add sex age etc

abunds_gl_ap<-merge(abunds_gl_ap,sample_mtd)

abunds_gl_apf<-abunds_gl_ap[ipsc_type=='population'&genotype%in%c('APOE33','APOE44')]

#Remove proteins that are not annotated
abunds_gl_apf<-abunds_gl_apf[!is.na(ensembl_gene_id_ver)][ensembl_gene_id_ver!='']

#Remove proteins identified only by a single peptide.
abunds_gl_apf<-abunds_gl_apf[`# Peptides`>1]


#Only astrocytes
abunds_gl_apf_astro<-abunds_gl_apf[cell_type=='Astrocyte']

atable <- AnalysisTableAnnotation$new()
atable$fileName = "sample_id"
atable$hierarchy[["protein_id"]] <- c("ensembl_gene_id_ver")
atable$factors[["genotype"]] = "genotype"
atable$factors[["batch"]] = "batch"
atable$factors[["sex"]] = "gender"

atable$set_response("aggr.abund")

config <- AnalysisConfiguration$new(atable)
adata <- setup_analysis(abunds_gl_apf_astro, config)


lfqdata <- prolfqua::LFQData$new(adata, config)
#Change 0 to NA
lfqdata$remove_small_intensities()

lfqplotter <- lfqdata$get_Plotter()
density_nn <- lfqplotter$intensity_distribution_density()
lfqplotter$NA_heatmap()

#Normalize protein intensities 
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq #transformed data

pl <- transformed$get_Plotter()
density_norm <- pl$intensity_distribution_density()
gridExtra::grid.arrange(density_nn, density_norm)

norm_trans_mat<-as.matrix(data.frame(data.table(transformed$to_wide()$data)[,-'isotopeLabel'],row.names = 'protein_id'))
#Rename with original sample name
mtd2<-data.table(transformed$factors())
colnames(norm_trans_mat)<-mtd2[colnames(norm_trans_mat),on='sample_name_corr']$sample_id

fwrite(data.table(norm_trans_mat,keep.rownames = 'protein_id'),fp(out,'robscaled_log2_abund_astro.csv.gz'))


#Linear model
transformed$config$table$get_response() # "transformedIntensity"

formula_Condition <-  strategy_lm("transformedIntensity ~ genotype+batch+sex")

modelName  <- "apoe_batch_sex"

Contrasts <- c("44vs33" = "genotypeAPOE44 - genotypeAPOE33")
mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_Id = transformed$config$table$hierarchy_keys() )

#Compute contrasts
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()
contr <- prolfqua::ContrastsModerated$new(contr)
contrdf <- contr$get_contrasts()
contrdt<-data.table(contrdf)

v2<-ggplot(contrdf)+geom_point(aes(x=diff,y=-log10(p.value),col=modelName))

prot_genes_infosf<-fread('metadata/proteins_genes_infos_filitered_unique_accession.csv')

contrdt<-merge(contrdt,unique(prot_genes_infosf[,protein_id:=ensembl_gene_id_ver][,.(protein_id,hgnc_symbol,Accession)],by='protein_id'))

#volcano plot
contrdt[,top10:=rank(p.value)<10]
ggplot(contrdt,aes(x=diff,y=-log10(p.value)))+
  geom_point(aes(col=-log10(p.value)))+
  scale_color_gradient2(high = '#56B1F7',mid ='#80B1F2' ,midpoint = 2,low = '#132B43')+
  geom_text_repel(aes(label=ifelse(top10,hgnc_symbol,'')),
                  max.overlaps = 2000,
                  size=3)+theme_minimal()


#top50 heatmap
library(pheatmap)

mat_astro<-as.matrix(data.frame(fread( 'outputs/01-iPSCpops_differential_proteomic_analysis/robscaled_log2_abund_astro.csv.gz'),row.names = 'protein_id'))

mat_mcc<-as.matrix(data.frame(fread( 'outputs/01-iPSCpops_differential_proteomic_analysis/robscaled_log2_abund_mcc.csv.gz'),row.names = 'protein_id'))


mtd<-fread('metadata/proteomics_samples_batch_metadata.csv')
mtd[,sample_id_corr:=ifelse(str_detect(sample_id,'^[0-9]'),ps('X',sample_id),sample_id)]
fwrite(mtd,'metadata/proteomics_samples_batch_metadata.csv')

res_de[,change:=ifelse(diff>0,'up in APOE44','down in APOE44')]

res_de[,top50:=rank(p.value)<=25,by=c('cell_type','change')]

deps_annot50<-unique(res_de[,deps:=paste(unique(cell_type[p.value<0.05]),collapse='&'),by=hgnc_symbol][(top50)][,.(hgnc_symbol,deps)])


#Only astrocytes
mtdf<-mtd[colnames(mat_astro),on='sample_id_corr']
expr_mat50<-mat_astro[res_de[(top50)&cell_type=='Astrocyte']$protein_id,]
rownames(expr_mat50)<-res_de[(top50)&cell_type=='Astrocyte']$hgnc_symbol

expr_mat50na<-expr_mat50
expr_mat50na[is.na(expr_mat50na)]<-min(expr_mat50,na.rm = T)*0.5

col_breaks<-(-30:30)/20

pheatmap(expr_mat50na,
               scale = 'row',
               annotation_col = data.frame(mtdf,row.names = 'sample_id_corr')[colnames(expr_mat50na),c('genotype'),drop=F],
              # annotation_row = data.frame(deps_annot50,row.names = 'hgnc_symbol'),
               color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                     "RdBu")))(length(col_breaks)-1),fontsize=6,
               breaks =col_breaks )



#fgsea
library(fgsea)


pathways_dir <- "/projectnb/tcwlab/MSigDB/"
gmt_mtd<-fread(file.path(pathways_dir,'gmt_metadata.csv')) 
pathways_info<-fread(file.path(pathways_dir,'all_CPandGOs_genesets_metadata.csv.gz'))

res_de<-fread('outputs/01-iPSCpops_differential_proteomic_analysis/res_batch_and_sex_corr_wald_test_moderated_44vs33_prolfqua_trans_on_raw_abund_gene_level_all_cell_type.csv.gz')


res_gsea_all<-Reduce(rbind,lapply(unique(res_def$cell_type), function(ct){
  
  message(paste('testing enrichement in',ct))
  
  res_de<-res_def[cell_type==ct]
  #calculate/extract the gene stats
  gene_stats<-setNames(res_de$statistic,res_de$hgnc_symbol)
  
  #run fgsea for every pathway source
  res_gsea<-Reduce(rbind,lapply(pathways_to_test, function(p){
    
    message(paste('testing enrichement for',gmt_mtd[name==p]$desc))
    
    pathways<- gmtPathways(file.path(pathways_dir,gmt_mtd[name==p]$gmt))
    
    res<-fgsea(pathways,
               stats=gene_stats,minSize=10,maxSize=2000,scoreType='std',nPermSimple = 10000)
    
    return(res[,source:=p])
    
  }))
  return(res_gsea[,cell_type:=ct])
  
}))

res_gsea_all<-merge(res_gsea_all,pathways_info,by=c('pathway'))[order(cell_type,source,subcat,pval)]


ggplot(res_gsea_all[padj<0.05])+geom_bar(aes(x=NES>0,fill=subcat))+facet_wrap('cell_type')+
  labs(y='number of enriched terms (padj<0.05)')

#Top pathways
res_gsea_top30<-res_gsea_all[padj<0.05,.SD[head(order(padj),40)],by=.(subcat,cell_type)]

emmaplot(head(res_gsea_top30[cell_type=='MCC'][order(pval)],100),label.size=2,cols=c('blue','white','red'))
#only up
emmaplot(head(res_gsea_all[padj<0.05&cell_type=='MCC'][order(pval)][NES>0],100),label.size=2,cols=c('blue','white','red'))

emmaplot(head(res_gsea_top30[cell_type=='Astrocyte'][order(pval)],100),label.size=2,cols=c('blue','white','red'))
#rep to eostradiol..because sex bias ?
table(mtd[ipsc_type=='population'&cell_type=='Astrocyte']$genotype,mtd[ipsc_type=='population'&cell_type=='Astrocyte']$gender)

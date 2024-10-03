library(prolfqua)

#Test different normalisation strategy

abunds_gl_ap<-fread('ref-data/proteomics_raw_data/gene_level_proteomics_data_tidy.csv')
sample_mtd<-fread('metadata/proteomics_samples_batch_metadata.csv') # add sex age etc
sample_mtd[,batch:=paste0('b',batch)]
abunds_gl_ap<-merge(abunds_gl_ap,sample_mtd)
 
#Remove all proteins that are identified only by a single peptide.
abunds_gl_apf<-abunds_gl_apf[`# Peptides`>1]


#44 vs 33 T test
abunds_gl_apf[,pct.det.ct:=1-(sum(is.na(log2.ratio))+sum(aggr.ratio==0,na.rm = T))/.N,by=.(cell_type,ensembl_gene_id_ver)]

abunds_gl_apf[,to.test:=pct.det.ct>=0.5&sum(aggr.ratio[genotype=='APOE44']>0,na.rm = T)>2&
                sum(aggr.ratio[genotype=='APOE33']>0,na.rm = T)>2,
              by=.(ensembl_gene_id_ver,cell_type)]

min_ra<-min(abunds_gl_apf$aggr.ratio[abunds_gl_apf$aggr.ratio!=0],na.rm = T)
abunds_gl_apf[log2.ratio==-Inf,log2.ratio:=log2(aggr.ratio+min_ra/2)]

ggplot(abunds_gl_apf)+geom_density(aes(x=log2.ratio,col=batch))

abunds_gl_apf[(to.test)&!is.na(log2.ratio),pval:=t.test(log2.ratio[genotype=='APOE44'],
                                                        log2.ratio[genotype=='APOE33'],alternative ='two.sided',var.equal =T)$p.value,
              by=.(ensembl_gene_id_ver,cell_type)]

abunds_gl_apf[(to.test)&!is.na(log2.ratio),t.stat:=t.test(log2.ratio[genotype=='APOE44'],
                                                          log2.ratio[genotype=='APOE33'],alternative ='two.sided',var.equal =T)$stat,
              by=.(ensembl_gene_id_ver,cell_type)]


abunds_gl_apf[(to.test),padj:=p.adjust(pval),by=.(cell_type)]

abunds_gl_apf[,avglog2_33:=mean(na.omit(log2.ratio[genotype=='APOE33'])),by=.(ensembl_gene_id_ver,cell_type)]
abunds_gl_apf[,avglog2_44:=mean(na.omit(log2.ratio[genotype=='APOE44'])),by=.(ensembl_gene_id_ver,cell_type)]
abunds_gl_apf[,log2FC:=avglog2_44-avglog2_33,by=.(ensembl_gene_id_ver,cell_type)]


res_de_log2<-unique(abunds_gl_apf[(to.test)][,.(cell_type,ensembl_gene_id_ver,t.stat,pval,padj,avglog2_33,avglog2_44,log2FC)])

#add gene name
res_de_log2<-merge(res_de_log2,unique(prot_genes_infosf[,.(ensembl_gene_id_ver,hgnc_symbol,Accession)],by='ensembl_gene_id_ver'))

fwrite(res_de_log2,fp(out,'res_t.test_44vs33_log2_expr_gene_level.csv.gz'))

res_de_log2<-fread(fp(out,'res_t.test_44vs33_log2_expr_gene_level.csv.gz'))


#Astrocytes
abunds_gl_apf_astro<-abunds_gl_apf[cell_type=='Astrocyte']
abunds_gl_apf_astro[,norm.abund:=aggr.ratio*median(aggr.pool.abund),by='ensembl_gene_id_ver']

plot(density(na.omit(abunds_gl_apf_astro$norm.abund)))


atable <- AnalysisTableAnnotation$new()
atable$fileName = "sample_id"
atable$hierarchy[["protein_id"]] <- c("ensembl_gene_id_ver")
atable$factors[["genotype"]] = "genotype"
atable$factors[["batch"]] = "batch"
atable$factors[["sex"]] = "gender"


atable$set_response("norm.abund")

#create the AnalysisConfiguration which needs the AnalysisTableAnnotation  and the AnalysisParameters.
config <- AnalysisConfiguration$new(atable)
config$table$fileName
adata <- setup_analysis(abunds_gl_apf_astro, config)

#Create the LFQData class instance and remove zeros from data (MaxQuant encodes missing values with zero).
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities()
lfqdata$factors()

#Visualization of not normalized data
lfqplotter <- lfqdata$get_Plotter()
density_nn <- lfqplotter$intensity_distribution_density()
lfqplotter$NA_heatmap()
lfqdata$get_Summariser()$plot_missingness_per_group()
lfqplotter$missigness_histogram()
#CV, Mean etc stat
stats <- lfqdata$get_Stats()
stats$violin()
stats$density_median()

#Normalize protein intensities 
lt <- lfqdata$get_Transformer()
transformed <- lt$log2()$robscale()$lfq
#?lt$log2()
transformed$config$table$is_response_transformed
pl <- transformed$get_Plotter()
density_norm <- pl$intensity_distribution_density()
gridExtra::grid.arrange(density_nn, density_norm)
#pl$pairs_smooth()

#fitting linear models 
transformed$config$table$get_response() # "transformedIntensity"

formula_Condition <-  strategy_lm("transformedIntensity ~ genotype+batch")

# specify model definition
modelName  <- "apoe44"
unique(transformed$data$genotype)
# [1] "APOE44" "APOE33"

Contrasts <- c("44vs33" = "genotypeAPOE44 - genotypeAPOE33")
mod <- prolfqua::build_model(
  transformed$data,
  formula_Condition,
  subject_Id = transformed$config$table$hierarchy_keys() )
mod$anova_histogram("FDR.Pr..F.")
#Compute contrasts
#non moderated
contr <- prolfqua::Contrasts$new(mod, Contrasts)
v1 <- contr$get_Plotter()$volcano()
v1$p.value
#moderated
contr <- prolfqua::ContrastsModerated$new(contr)
contrdf <- contr$get_contrasts()
head(contrdf)
contrdt<-data.table(contrdf)
unique(contrdf$modelName)
v2<-ggplot(contrdf)+geom_point(aes(x=diff,y=-log10(p.value),col=modelName))
v1$p.value/v2
#save
#add gene name
prot_genes_infosf<-fread('metadata/proteins_genes_infos_filitered_unique_accession.csv')

contrdt<-merge(contrdt,unique(prot_genes_infosf[,protein_id:=ensembl_gene_id_ver][,.(protein_id,hgnc_symbol,Accession)],by='protein_id'))
fwrite(contrdt,fp(out,'res_batch_corr_wald_test_moderated_44vs33_prolfqua_trans_on_norm_abund_gene_level_astro.csv.gz'))

#same norm mat
norm_trans_mat<-as.matrix(data.frame(data.table(transformed$to_wide()$data)[,-'isotopeLabel'],row.names = 'protein_id'))
#rename with original sample name
mtd2<-data.table(transformed$factors())
mtd2[,sample_name_corr:=str_replace_all(sampleName,'~','.')]
colnames(norm_trans_mat)<-mtd2[colnames(norm_trans_mat),on='sample_name_corr']$sample_id

fwrite(data.table(norm_trans_mat,keep.rownames = 'protein_id'),fp(out,'robscaled_log2gisnormalized_abund_astro.csv.gz'))


#volcano plot
contrdt[,top10:=rank(p.value)<10]
ggplot(contrdt,aes(x=diff,y=-log10(p.value)))+
  geom_point(aes(col=-log10(p.value)))+
  scale_color_gradient2(high = '#56B1F7',mid ='#80B1F2' ,midpoint = 2,low = '#132B43')+
  geom_text_repel(aes(label=ifelse(top10,hgnc_symbol,'')),
                  max.overlaps = 2000,
                  size=3)+theme_minimal()

#
#save both
res_deps<-rbind(contrdt[,cell_type:='Astrocyte'][,-'top10'],
                fread('outputs/01A-test_normalisation_and_modeling/res_batch_corr_wald_test_moderated_44vs33_prolfqua_trans_on_norm_abund_gene_level_mcc.csv.gz')[,cell_type:='MCC'])
fwrite(res_deps,fp(out,'res_batch_corr_wald_test_moderated_44vs33_prolfqua_trans_on_norm_abund_gene_level_all_celltype.csv.gz'))

fwrite(res_deps[cell_type=='MCC'],fp(out,'res_batch_corr_wald_test_moderated_44vs33_prolfqua_trans_on_norm_abund_gene_level_mcc.csv'))
fwrite(res_deps[cell_type=='Astrocyte'],fp(out,'res_batch_corr_wald_test_moderated_44vs33_prolfqua_trans_on_norm_abund_gene_level_astrocyte.csv'))

#save tidy format with GIS norm abundance
abunds_gl_apf<-rbind(abunds_gl_apf_mcc,abunds_gl_apf_astro)
fwrite(abunds_gl_apf,'outputs/01A-test_normalisation_and_modeling/gis_normalized_filtered_data_tidy.csv.gz')
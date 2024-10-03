library(WGCNA)
library(edgeR)
library(DESeq2)

#Run WGCNA on RNA then on Proteo and compare module conservation using eigengenes
#WGCNA::modulePreservation() assesses network module preservation across cohorts. We also used this function to assess the effect of missing values on the consensus network. Zsummary composite preservation scores were obtained using the consensus network as the template versus each other cohort or missing value threshold tested, with 500 permutations. Random seed was set to 1 for reproducibility, and the quickCor option was set to 0. We also assessed network module preservation using synthetic eigenproteins. In brief, protein module members in the consensus network template with a kME.intramodule among the top 20th percentile were assembled into a synthetic module in each target cohort, and synthetic modules with at least four members were used to calculate synthetic weighted eigengenes representing the variance of all members in the target network across case samples via the WGCNA::moduleEigengenes() function. Statistics and correlation scatter plots involving target cohort traits were then calculated and visualized'.
#method stats in details here : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3024255/


#APOE4 related RNA module annotation
apoe4_rna<-c('yellow','red','lightcyan','darkgreen','orange','green','darkslateblue','darkgoldenrod','magenta')
res_annot<-fread('outputs/02B-module_conservation/RNA_modules/wgcna_RNA_modules_neurons_terms_enriched_padj0.1_overlap5.csv.gz')
res_annot[query.%in%apoe4_rna][,head(.SD[order(padj)]$term,20),by='query']


sel_pathws<-list(yellow=c('GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY','REACTOME_LONG_TERM_POTENTIATION',
                          'GOBP_SYNAPTIC_VESICLE_PRIMING'),#up
                 red=c('GOBP_RNA_SPLICING','GOBP_RNA_PROCESSING'),#up, 
                 lightcyan=c('GOMF_EXTRACELLULAR_MATRIX_BINDING','NABA_CORE_MATRISOME',
                             'GOBP_COLLAGEN_FIBRIL_ORGANIZATION'), #up
                 darkgreen=c('GOCC_SYNAPTIC_VESICLE_MEMBRANE',
                             'GOBP_SYNAPTIC_SIGNALING','GOBP_NEUROTRANSMITTER_TRANSPORT'),#up, remaining are down
                 orange=c('GOBP_NEGATIVE_REGULATION_OF_BIOSYNTHETIC_PROCESS','GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY'),
                 green=c('KEGG_LYSOSOME','GOCC_SECRETORY_GRANULE','GOBP_LIPID_METABOLIC_PROCESS'),
                 darkslateblue=c('KEGG_RIBOSOME','GOBP_CYTOPLASMIC_TRANSLATION','REACTOME_CELLULAR_RESPONSE_TO_STARVATION'),
                 darkgoldenrod=c('REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR','GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS'),
                 magenta=c('GOBP_ATP_METABOLIC_PROCESS','REACTOME_GLUCOSE_METABOLISM'))

sel_pathws<-rbindlist(lapply(apoe4_rna, function(m)data.table(term=sel_pathws[[m]],query=m)))


ggplot(merge(res_annot,sel_pathws)[,query:=factor(query,levels =apoe4_rna )])+
  geom_col(aes(x=-log10(padj),y=term,fill=query))+facet_wrap('query',ncol = 1,scales = 'free')+
  scale_fill_manual(values = apoe4_rna)+theme_minimal()+
  theme(axis.text.y = element_text(size = 7))



####Proteomic Module
res_enr<-fread('outputs/02B-module_conservation/wgcna_prot_modules_MCC_pathways_enrichment.csv.gz')
res_enr5<-res_enr[padj<0.25,head(.SD[order(pval)],5),by='query']

enr5_mat<-dcast(res_enr5,query~pathway,value.var = 'pval')
enr5_mat<-(-log10(enr5_mat))

col_breaks<-c(0,0.2,0.4,0.8,1,1.5,2,3,4,6,8,10,12,15,18,20,25,30)

pdf(fp(out,'fig1-heatmap_regulons_GOKEGG.pdf'),width = 5,height = 8.5)
print(pheatmap(t(top5go_mat),
               color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                     "RdYlBu")))(length(col_breaks)-1),
               fontsize=6,
               breaks =col_breaks ))
dev.off()


####Modules conservation

CreateJobForRfile('scripts/02Biii-run_module_conservation_rna_protein.R')
RunQsub('scripts/02Biii-run_module_conservation_rna_protein.R')

CreateJobForRfile('scripts/02Biii-run_module_conservation_rnaNeuroProp_protein.R')
RunQsub('scripts/02Biii-run_module_conservation_rnaNeuroProp_protein.R',job_name = 'conservationRNANeuroProp')


mp<-readRDS('outputs/02B-module_conservation/RNA_modules_neuro_prop/rna_module_conservation_in_prot.rds')

resmp = merge(data.table(mp$preservation$observed[[1]][[2]][, 1:2],keep.rownames = 'module'), data.table(mp$preservation$Z[[1]][[2]][, 1:2],keep.rownames = 'module'))
resmp<-merge(resmp,
             data.table(module=c('yellow','red','darkgreen', 'orange','green','darkslateblue','deeppink'),
                        activity=c('synaptic_plasticity','rna_splicing','ER_stress','ribosome/translation','lysosome/lipid_catabolism/cilium','catabolism/vesicles','synapse')),
             all.x=T,by='module')

resmp<-resmp[!module%in%c('random','grey')]
resmp[,module:=factor(module,levels = module[order(-medianRank.pres)])]


ggplot(resmp,aes(y=Zsummary.pres,x=module))+
  geom_hline(yintercept = 2,col='blue',linetype='dashed')+
  geom_hline(yintercept = 10,col='darkgreen',linetype='dashed')+
  geom_point(aes(col=module,size=moduleSize))+
  scale_color_manual(values = ifelse(levels(resmp$module)=='lightcyan','lightcyan3',levels(resmp$module)))+
  scale_size(limits = c(5,1000))+scale_x_discrete(guide = guide_axis(angle=60))+
  geom_text_repel(aes(label=ifelse(Zsummary.pres>3,activity, '')),size=3)+
  theme_bw()
run_prolfqua <- function(df, variable, model_name, contrast) {
  
  #Define Column names for configuration
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "Sample"
  atable$workIntensity = "Intensity"
  atable$hierarchy[["protein_Id"]] = "protein_Id"
  atable$factors[["Group"]] = variable
  
  config <- prolfqua::AnalysisConfiguration$new(atable)
  
  # Build LFQData object
  analysis_data <- prolfqua::setup_analysis(df, config)
  lfqdata <- prolfqua::LFQData$new(analysis_data, config)

  #Remove NAs
  #lfqdata$remove_small_intensities()
  
  #Normalize protein intensities 
  lt <- lfqdata$get_Transformer()
  transformed <- lt$log2()$robscale()$lfq
  transformed$config$table$is_response_transformed
  
  write.csv(as.data.frame(transformed$data), paste0(main_dir, "/output/DEA/", model_name, "_prolfqua_normalized.csv"), row.names = F, quote = F)
  
  formula_Condition <-  strategy_lm("transformedIntensity ~ Group")

  # specify model definition
  modelName  <- model_name
  
  # Filter for proteins that are present in at least two groups
  valid_proteins <- transformed$data %>%
    dplyr::group_by(protein_Id, Group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::count(protein_Id) %>%
    dplyr::filter(n >= 2) %>%
    dplyr::pull(protein_Id)
  
  filtered_data <- transformed$data %>%
    dplyr::filter(protein_Id %in% valid_proteins)
  
  mod <- prolfqua::build_model(
    filtered_data,
    formula_Condition,
    subject_Id = transformed$config$table$hierarchy_keys() )
  
  ##### Compute Contrasts ##### 
  #Method1
  
  contr <- prolfqua::Contrasts$new(mod, contrast)
  contrdf <- contr$get_contrasts()
  write.csv(contrdf, paste0(main_dir, "/output/DEA/", model_name, "_prolfqua_contrast_dea.csv"), row.names = F, quote = F)
  
  #Method2
  contr2 <- prolfqua::ContrastsModerated$new(contr)
  contrdf2 <- contr$get_contrasts()
  write.csv(contrdf2, paste0(main_dir, "/output/DEA/", model_name, "_prolfqua_contrast_moderated_dea.csv"), row.names = F, quote = F)
}


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
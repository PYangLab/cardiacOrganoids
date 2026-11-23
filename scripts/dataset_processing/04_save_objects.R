# split seurat objects into time points and chambers / individual objects and save:

# SINGLE CELL DATASETS:

giacomelli <- seu_list_sc_filtered_annotated$giacomelli

zhang <- seu_list_sc_filtered_annotated$zhang

schmidt <- seu_list_sc_filtered_annotated$schmidt
schmidt.atria <- schmidt[,schmidt$orig.ident == "Atria_d9p5_TNI_222826.e5"]
schmidt.avc <- schmidt[,schmidt$orig.ident %in% c("AVC_d9p5_TNI_222832.e7","AVC_d9p5_TNT2_222827.e6")]
schmidt.lv <- schmidt[,schmidt$orig.ident %in% c("LV_d7p5_TNI_222823.e5","LV_d7p5_TNT2_222828.e6")]
schmidt.rv <- schmidt[,schmidt$orig.ident %in% c("RV_d9p5_TNI_222824.e5","RV_d9p5_TNT2_222829.e6")]
schmidt.oft <- schmidt[,schmidt$orig.ident %in% c("OFT_d9p5_TNI_222825.e5","OFT_d9p5_TNT2_222830.e6")]

ho <- seu_list_sc_filtered_annotated$ho
hoD8 <- ho[,ho$day == "d8"]
hoD14 <- ho[,ho$day == "d14"]
hoD21 <- ho[,ho$day == "d21"]

abilez <- seu_list_sc_filtered_annotated$abilez

feng <- seu_list_sc_filtered_annotated$feng
fengD15.atria <- feng[,feng$treatment == "RA+" & feng$stage == "D15"]
fengD30.atria <- feng[,feng$treatment == "RA+" & feng$stage == "D30"]
fengD15.ventricle <- feng[,feng$treatment == "RA-" & feng$stage == "D15"]
fengD30.ventricle <- feng[,feng$treatment == "RA-" & feng$stage == "D30"]

zou <- seu_list_sc_filtered_annotated$zou
zouD25.pd <- zou[,zou$condition == "SCPD_D25.5"]
zouD25.pd$replicate <- NULL # remove empty metadata slot to avoid downstream errors
zouD42.pd <- zou[,zou$condition == "SCPD_D42.5"]
zouD42.pd$batch <- zouD42.pd$replicate
zouD42.npd <- zou[,zou$condition == "SCWOPD_D42.5"]
zouD42.npd$batch <- zouD42.npd$replicate

ohern <- seu_list_sc_filtered_annotated$ohern
ohernD15.mac <- ohern[,ohern$condition == "hHO_MPs"]
ohernD15.non.mac <- ohern[,ohern$condition == "Control"]
ohernD26.mac <- ohern[,ohern$condition == "hHO_MPs_D26"]
ohernD26.non.mac <- ohern[,ohern$condition == "hHO_D26"]

sc_dataset_list <- list(giacomelli = giacomelli,
                     zhang = zhang,
                     hoD8 = hoD8,
                     hoD14 = hoD14,
                     hoD21 = hoD21,
                     abilez = abilez,
                     fengD15.atria = fengD15.atria,
                     fengD15.ventricle = fengD15.ventricle,
                     fengD30.atria = fengD30.atria,
                     fengD30.ventricle = fengD30.ventricle,
                     schmidt.atria = schmidt.atria,
                     schmidt.avc = schmidt.avc,
                     schmidt.lv = schmidt.lv,
                     schmidt.oft = schmidt.oft,
                     schmidt.rv = schmidt.rv,
                     zouD25.pd = zouD25.pd,
                     zouD42.pd = zouD42.pd,
                     zouD42.npd = zouD42.npd,
                     ohernD15.mac = ohernD15.mac,
                     ohernD15.non.mac = ohernD15.non.mac,
                     ohernD26.mac = ohernD26.mac,
                     ohernD26.non.mac = ohernD26.non.mac
)

sc_dataset_list <- Map(function(obj, dataset_name) {
  if (grepl("hca", dataset_name, ignore.case = TRUE)) {
    return(obj)
  }
  
  celltype <- as.character(obj$celltype)
  
  celltype_counts <- table(celltype)
  rare_celltypes <- names(celltype_counts)[celltype_counts < 20]
  
  celltype[
    celltype == "Lymphoid" |
      celltype %in% rare_celltypes
  ] <- "unassigned"
  
  obj$celltype <- celltype
  return(obj)
}, sc_dataset_list, names(sc_dataset_list))

# saveRDS(sc_dataset_list, "~/Cardiac_organoids/processed_data/sc_organoids.rds")




# SINGLE NUCLEI DATASETS:

pocock <- seu_list_sn_filtered_annotated$pocock

voges <- seu_list_sn_filtered_annotated$voges

landau <- seu_list_sn_filtered_annotated$landau
landau.mac <- landau[,landau$condition == "mac"]
landau.non.mac <- landau[,landau$condition == "non"]

sn_dataset_list <- list(voges = voges,
                     pocock = pocock,
                     landau.mac = landau.mac,
                     landau.non.mac = landau.non.mac)

sn_dataset_list <- Map(function(obj, dataset_name) {
  if (grepl("hca", dataset_name, ignore.case = TRUE)) {
    return(obj)
  }
  
  celltype <- as.character(obj$celltype)
  
  celltype_counts <- table(celltype)
  rare_celltypes <- names(celltype_counts)[celltype_counts < 20]
  
  celltype[
    celltype == "Lymphoid" |
      celltype %in% rare_celltypes
  ] <- "unassigned"
  
  obj$celltype <- celltype
  return(obj)
}, sn_dataset_list, names(sn_dataset_list))

# saveRDS(sn_dataset_list, "~/Cardiac_organoids/processed_data/sn_organoids.rds")
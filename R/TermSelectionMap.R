# 1. Term Selection Map: 
#     inputs: list of terms
#     output: matrix of terms annotated for "RNA","Protein","Immune",  "DNA/rep",  "Cancer",  "Metab", "ECM", "Stress", "EXSig",  "Apoptosis", "Angio"
# TermSelectionMap()

TermSelectionMap=function(topGSTerms, termi=1){
  ProTerms=c("TRANSLATION", "RIBOSOME", "43S", "PROTEOL", "PEPTIDE", "EIF", "GOLGI", "ELONGATION", "UB", "UBIQUITIN", "LYSOSOME")
  RNAnames=c("_RNA", "TRANSCRIPT", "POL", "PROMOTER", "SPLIC", "EXON", "INTRON", "NONSENSE_MEDIATED", "3_UTR", "5_UTR", "GENE_EXPRESSION", "PPARA", "MRNA_", "_MRNA")
  
  ImmTerms=c("IFN", "IL[0-9]","IL_[0-9]", "TNF", "CYTOKINE", "NKCELL", "ANTIGEN", "CHEMOKINE", "TCR", "BCR", "CYTOTOXIC",  "IMMUNE", "TH1", "STAT[0-9]", "INFLAMMATORY", "COMPLEMENT", "COAGULATION", "ALLOGRAFT","INTERFERON", "IGG","IGA_", "BCELL", "TCELL", "LYMPHOCYTE", "CTL", "TH2", "GRANZYME", "NEUTROPHIL", "T_CELL", "LEUKOCYTE", "B_CELL", "TREG", "INFLUENZA", "FC_", "IGF", "CLASSIC", "COMP_", "THROMBIN", "DEFENSIN", "ASTHMA", "LUPUS", "PDGF",  "TOLL", "ILS", "CCR", "MHC", "THELP", "TCAPOPTOSIS", "GMCSF", "GM_CSF", "CXCR", "CXCL", "MAST_CELL", "GRANULOCYTE", "TCAPOPTOSIS", "MONOCYTE", "INTEGRIN", "IMMUNO", "PD1_SI", "VIRION", "AUTOIMMUNE", "ALLOGRAFT","LYM_PATH",
             "CTL_PATHWAY","EPHA4_PATHWAY", "VIRAL", "PLATELET", "NFAT_", "HIV_", "GRAFT", "TCPTP", "TID_",
             "INTERLEUKIN", "RUNX3", "HEMATOPOIETIC") #TID?
  DNAnames=c("DNA","G0", "G1", "G2", "MITO", "MEIO", "REPLIC", "REPAIR", "CELL_CYCLE", "TELOMER", "CENTROSOME", "HDAC", "AURORA", "CDC", "MCM_PATHWAY", "RANMS","PLK1",
             "ANAPHASE", "APC", "M_PHASE", "CYCLIN_A", "PROPHASE", "CHROMOSOMES", "S_PHASE")
  CancerTerms=c( "RAS_", "ECADHERIN", "ARMS", "EGFR", "TRKA", "CTNNB1", "WNT", "NFKB", "APC_", "MYC_", "P63", "PIP2","ATM_","NGF","ERK_",
                 "ERBB", "NOTCH", "JAK", "_ATR_", "BRCA", "CATENIN", "PI3K", "AKT", "MAPK", "FGFR", "TGF", "_ATM_", "YAP","RAB_", "CSK_PATHWAY",
                 "HEDGEHOG", "FOX", "HER2" ,"_AR_", "MTOR", "P38", "P73", "TRAF6", "PI_3K_", "HNF3B", "PHOSPHOLIPASE", "CADHERIN", "KIT_PATHWAY", "P53_", "HORMONE_RECEPTOR", "SHC", "PTEN", "CTCF",
                 "_P53", "AXIN", "CANCER", "ESR_MEDIATED_SIGNALING", "ESTROGEN")
  EpiGen=c("METHYLATION", "CTCF")
  
  #SignallingTerms=c("SIGNALING", "PATHWAY", "ERK_")
  ECMTerms=c("ECM", "EXTRACELLULAR", "MMP", "COLLAGEN", "MATRISOME", "KERATAN", "BASEMENT", "FRS2_", "ADHERENS", "PROTEOGLYCAN", "CELL2CELL", "ACTIN", "JUNCTION", "ADHESION", "INTEGRIN", 
             "LAMININ", "KERATIN", "AMYLOID", "CORNIFIED_ENVELOPE")
  StressTerms=c("STRESS", "HIF", "HYPOXIA", "UPR", "UNFOLDED_PROTEIN", "ARF6", "ARF4")
  TFTerms=c("TF", "AR_TF")
  MetabTerms=c("METABOL", "GLUCOSE", "GLYCO", "PENTOSE", "TCA_", "PYRUVATE", "LIPID", "MITOCHON", "OXIDAT", "CARBOHYDRATE", "CREB", "CATABOL", "CHOLESTER", "DIGESTION", "GLUCURONIDATION", "TRIACYLGLYCEROL", "BIOSYNTHESIS", "GLUCONEO", "DIABETES", "RECYCLING", "BETA_CELL", "PTP1B", "INSULIN", "ELECTRON")
  #UBTerms=c("UBIQ", "PROTEASOME", "AUTOPHAGY", "INSULIN")
  APOPTOSIS=c("APOPTOSIS", "BCL", "D4GDI", "CELL_DEATH", "SENESCENCE")
  ANGIO=c("VEGF")
  EXTRASIG=c("NEUROTRANSMITTER", "TRANSPORTER", "SECRETED", "GABA_", "CALCIUM", "K_CHANNEL","TRANSMEMBRANE",
             "INTEGRIN", "GPCR", "MUSCLE", "NMDA_RCEPTOR_GLUTAMATE", "DARP32", "OPIOD", "ENDOCYTOSIS",  "G_ALPHA", 
             "GAP_JUNCTION", "ION_CHANNEL", "RECEPTOR_TYROSINE_KINASES", "CELLULAR_RESPONSES", "AXON_GUIDANCE",
             "CELL_SURFACE_INTERACTIONS", "PCP")
  DIFFERENTIATION=c("SEMA3A")
  
  ImmTerms2=topGSTerms[unlist(sapply(ImmTerms, function(x) grep(x, topGSTerms)))]
  ProTerms2=topGSTerms[unlist(sapply(ProTerms, function(x) grep(x, topGSTerms)))]
  DNATerms2=topGSTerms[unlist(sapply(DNAnames, function(x) grep(x, topGSTerms)))]
  RNATerms2=topGSTerms[unlist(sapply(RNAnames, function(x) grep(x, topGSTerms)))]
  Methyl=topGSTerms[unlist(sapply(EpiGen, function(x) grep(x, topGSTerms)))]
  CancerTerms2=topGSTerms[unlist(sapply(CancerTerms, function(x) grep(x, topGSTerms)))]
  #SignallingTerms2=topGSTerms[unlist(sapply(SignallingTerms, function(x) grep(x, topGSTerms)))]
  ECMTerms2=topGSTerms[unlist(sapply(ECMTerms, function(x) grep(x, topGSTerms)))]
  StressTerms2=topGSTerms[unlist(sapply(StressTerms, function(x) grep(x, topGSTerms)))]
  MetabTerms2=topGSTerms[unlist(sapply(MetabTerms, function(x) grep(x, topGSTerms)))]
  #UBTerms2=topGSTerms[unlist(sapply(UBTerms, function(x) grep(x, topGSTerms)))]
  APOP=topGSTerms[unlist(sapply(APOPTOSIS, function(x) grep(x, topGSTerms)))]
  ANG=topGSTerms[unlist(sapply(ANGIO, function(x) grep(x, topGSTerms)))]
  EXS=topGSTerms[unlist(sapply(EXTRASIG, function(x) grep(x, topGSTerms)))]
  
  matrix2=matrix(0, ncol=12, nrow=length(topGSTerms))
  rownames(matrix2)=topGSTerms
  colnames(matrix2)=c("RNA","Protein","Immune",  "DNA/rep",  "Cancer","Metab", "ECM", "Stress", "EXSig",  "Apoptosis", "Angio", "Epi")
  matrix2[match(ImmTerms2, rownames(matrix2)), 3]=termi
  matrix2[ match(DNATerms2, rownames(matrix2)), 4]=termi
  matrix2[ match(RNATerms2, rownames(matrix2)),1]=termi
  matrix2[ match(ProTerms2, rownames(matrix2)),2]=termi
  matrix2[ match(CancerTerms2, rownames(matrix2)),5]=termi
  matrix2[ match(ECMTerms2, rownames(matrix2)),7]=termi
  matrix2[ match(StressTerms2, rownames(matrix2)),8]=termi
  matrix2[ match(MetabTerms2, rownames(matrix2)),6]=termi
  matrix2[ match(EXS, rownames(matrix2)),9]=termi
  matrix2[ match(APOP, rownames(matrix2)),10]=termi
  matrix2[ match(ANG, rownames(matrix2)),11]=termi
  matrix2[ match(Methyl, rownames(matrix2)),11]=termi
  
  matrix2
}

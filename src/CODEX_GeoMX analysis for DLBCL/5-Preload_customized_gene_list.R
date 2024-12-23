
# define term and gene
Macro_Cytoplasmic_Translation  <- c("RPL4", "RPL30", "RPL32", "RPL31", "RPL34", "RPLP1", "RPLP0", 
                                    "RPL10A", "RPL8", "RPL9", "RPL6", "RPL7", "RPL7A", "RPS14", 
                                    "RPS17", "RPS16", "RPL18A", "RPS19", "RPS18", "RACK1", 
                                    "RPLP2", "RPL38", "RPL37", "RPS11", "RPS10", "RPS9", "RPL21", 
                                    "RPS8", "RPL22", "RPS6", "RPS3A", "RPSA", "RPL37A", "RPL24", 
                                    "RPL27", "RPL26", "RPL29", "RPL28", "RPL12", "RPL11", "RPS15A", 
                                    "RPL14", "RPS3", "RPL13", "RPL15", "RPS2", "RPL18", "RPS27A", 
                                    "RPL17", "RPL19", "RPL41", "RPL35A", "RPS26", "RPS25", "RPS28", 
                                    "RPS29", "RPL27A", "RPS20", "FAU", "RPS21", "RPS24")
Macro_Peptide_Biosynthetic_Process <- c("RPL4", "RPL30", "RPL32", "RPL31", "RPL34", "RPLP1", 
                                        "RPLP0", "MRPL37", "RPL10A", "RPL8", "RPL9", "MRPL35", 
                                        "RPL6", "RPL7", "RPS14", "GGTLC1", "RPL7A", "RPS17", "RPS16", 
                                        "RPL18A", "RPS19", "RPS18", "RPLP2", "RPL38", "RPL37", 
                                        "RPS11", "RPS10", "RPS9", "RPL21", "RPS8", "RPL22", 
                                        "RPS6", "RPSA", "RPS3A", "RPL37A", "RPL24", "RPL27", "RPL26", 
                                        "RPL29", "RPL28", "RPL12", "RPL11", "MRPL10", "RPS15A", "RPL14", 
                                        "RPS3", "RPL13", "RPL15", "RPS2", "RPL18", "RPS27A", "RPL17", 
                                        "RPL19", "RPL41", "RPL35A", "RPS26", "RPS25", "RPS28", 
                                        "TNIP1", "RPS29", "RPL27A", "RPS20", "FAU", "RPS21", "RPS24")

Macro_merged_translation<- unique(Macro_Cytoplasmic_Translation,Macro_Peptide_Biosynthetic_Process )

Macro_Macrophage_differentiation <- c("CSF1", "ID2", "RIPK1", "C1QC")
Macro_Gene_expression <- c("RPL4", "RPL30", "RPL32", "RPL31", "RPL34", "RPLP1", 
                           "RPLP0", "MRPL37", "RPL8", "RPL10A", "RPL9", "MRPL35", 
                           "RPL6", "RPL7", "RPS14", "RPL7A", "RPS17", "RPS16", 
                           "RPL18A", "RPS19", "RPS18", "MAGOH", "RPLP2", "RPL38", 
                           "RPL37", "RPS11", "RPS10", "RPS9", "RPL21", "RPS8", 
                           "RPL22", "RPS6", "RPSA", "RPS3A", "HNRNPUL1", 
                           "DNAJB11", "RPL37A", "RPL24", "SRSF3", "RPL27", 
                           "RPL26", "RPL29", "RPL28", "RPL12", "RPL11", "MRPL10", 
                           "EXOSC10", "RPS15A", "RPL14", "RPS3", "RPL13", "RPL15", "RPS2", 
                           "RPL18", "RPS27A", "RPL17", "RPL19", "RPL41", "GALNT2", 
                           "ALYREF", "RPL35A", "RPS26", "SETX", "RPS25", "RPS28", 
                           "TNIP1", "RPS29", "RPL27A", "RPS20", "FAU", "RPS21", "RPS24")

Macro_TOR_signaling <- c("PRR5", "RPS6KB1", "RPS6", "LAMTOR4", "CARD11")
Macro_Regulation_T_Cell_Differentiation <- c("IL4I1", "DUSP10", "HLA.DRA", "SOX12", "HLA.G", "LILRB4", "HLA.DRB1")
Macro_MHC_class_II <- c("HLA.DRA", "IFI30", "HLA.DOA", "CTSD", "HLA.DQB2", "HLA.DRB1","CD74")
Macro_NF_kappaB <- c("ALK", "NPM1", "TFRC", "CARD9", "IL18", "CIB1", "TNFRSF11A", "RFPL4A", "TRIM6", "DDRGK1", "RPS3", "FLOT2", "MAP3K13", "PPIA", "CARD11", "IL18R1")
Macro_Response_To_Hypoxia <- c("RPTOR", "SUV39H2", "PTGIS", "KCND2", "TERT", "MLST8", "RORA", "ADO", "HIF1A", "BMP7")
Macro_Negative_Regulation_Of_Inflammatory_Response <- c("GHSR", "IL22", "APCS", "PTGIS", "MACIR", "MVK", "RORA", "NLRC3", "CDH5", "LRFN5", "MDK", "ELF4", "MKRN2", "TRIM65")
Macro_Signal_Transduction_In_Absence_Of_Ligand <- c("IL1A", "CSF2", "TERT", "IFI6", "BCL2")
# tumor to macrophage
Macro_Negative_Regulation_Of_Tumor_Necrosis_Factor_Production <- c("GHSR", "ARG2", "ELF4", "EHHADH", "NLRC3", "BPI", "DICER1", "LILRA4", "HAVCR2")
Macro_Positive_Regulation_Of_T_Cell_Activation <- c("SMARCE1", "CD274", "SMARCB1", "TFRC", "HLA.A", "AIF1", "CCDC88B", "HLA.DRA", 
                                                    "FLOT2", "HLA.DOA", "BRD7", "B2M", "HLA.DRB1", "HLA.DQB2")
Macro_MHC_Class_II_Protein_Complex_Assembly <- c("HLA.DRA", "HLA.DOA", "HLA.DQB2", "HLA.DRB1")
Macro_customized_M1M2 <- c("C1QA", "CD163", "MRC1", "STAT6", "IRF4", "ARG1", 
                           "STAT1", "IL12A", "IL23A", "NFKB1", "NOS2")
# tumor signatures
Tumor_Regulation_Of_Macrophage_Proliferation <- c("IL34", "MAPK1", "PTK2")
Tumor_Regulation_Of_Macrophage_Chemotaxis <- c("UXT", "MDK", "IL34", "MTUS1", "MAPK1", "PTK2")
Tumor_Antigen_Processing_And_Presentation_Of_Exogenous_Peptide_Antigen_Via_MHC_Class_II <- c("HLA.DMB", "FCER1G", "HLA.DRA", "IFI30", "HLA.DOA", "CTSD", 
                                                                                             "HLA.DQB2", "HLA.DRB1", "HLA.DPA1", "DNM2")
Tumor_Peptide_Antigen_Assembly_With_MHC_Protein_Complex <- c("HLA.DMB", "HLA.DRA", "CALR", "HLA.DOA", "HLA.DQB2", 
                                                             "HLA.DRB1", "TAPBP", "HLA.DPA1")
Tumor_Negative_Regulation_Of_Response_To_Type_II_Interferon <- c("NR1H2", "OTOP1", "PARP14")
Tumor_IFNG_cascade <- c( "AKT1", "AKT2", "AKT3", "BCL2L1", "CBL", "CBLB", "CBLC", "CCND1", "CCND2", "CCND3", "CISH", "CLCF1",
                         "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9",
                         "CNTF", "CNTFR", "CREBBP", "CRLF2", "CSF2", "CSF2RA", "CSF2RB", "CSF3", "CSF3R", "CSH1", "CTF1", 
                         "EP300", "EPO", "EPOR", "GH1", "GH2", "GHR", "GRB2", "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", 
                         "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNAR1", "IFNAR2", "IFNB1",
                         "IFNE", "IFNG", "IFNGR1", "IFNGR2", "IFNK", "IFNL1", "IFNL2", "IFNL3", "IFNLR1", "IFNW1", "IL10",
                         "IL10RA", "IL10RB", "IL11", "IL11RA", "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IL13", "IL13RA1", 
                         "IL13RA2", "IL15", "IL15RA", "IL19", "IL2", "IL20", "IL20RA", "IL20RB", "IL21", "IL21R", "IL22", 
                         "IL22RA1", "IL22RA2", "IL23A", "IL23R", "IL24", "IL26", "IL2RA", "IL2RB", "IL2RG", "IL3", "IL3RA", 
                         "IL4", "IL4R", "IL5", "IL5RA", "IL6", "IL6R", "IL6ST", "IL7", "IL7R", "IL9", "IL9R", "IRF9", "JAK1", 
                         "JAK2", "JAK3", "LEP", "LEPR", "LIF", "LIFR", "MPL", "MYC", "OSM", "OSMR", "PIAS1", "PIAS2", "PIAS3", 
                         "PIAS4", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R5", "PIM1", "PRL",
                         "PRLR", "PTPN11", "PTPN6", "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS7", "SOS1", "SOS2", 
                         "SPRED1", "SPRED2", "SPRY1", "SPRY2", "SPRY3", "SPRY4", "STAM", "STAM2", "STAT1", "STAT2", "STAT3", 
                         "STAT4", "STAT5A", "STAT5B", "STAT6", "TPO", "TSLP", "TYK2")
Tumor_IRF_family <- c("IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9")
Tumor_Positive_Regulation_Of_Leukocyte_Mediated_Cytotoxicity <- c("KLRC2", "NOS2", "IL23R", "KLRC4", "HLA.B", "HLA.C", "CD1E", 
                                                                  "HLA.A", "HLA.G", "NCR3", "HLA.DRA", "ULBP1", "HLA.DRB1")
Tumor_customized_M1M2_formation <- c("IL4", "IL13", "IL10", "TGFB1", "TNF", "IFNG")
Tumor_Innate_Immune_System_R_HSA_168249 <- c("PGLYRP3", "NCF1", "NCF2", "WIPF1", "CLEC10A", "PROS1", "WIPF3", 
                  "TXN2", "LGALS3", "PSMD4", "GRAP2", "FTH1", "TOM1", "COTL1", 
                  "ATP6V1E1", "ACAA1", "GYG1", "B2M", "SKP1", "MEF2A", "ARSA", 
                  "MAP2K1", "CHRNB4", "IST1", "S100A1", "FBXW11", "DEFB116", 
                  "NANOS3", "PRKCD", "DEFB114", "HLA.B", "HLA.C", "CNPY3", 
                  "HLA.A", "PLPP5", "DOK3", "PLPP4", "ELMO1", "PSME2", "CD300E", 
                  "PRTN3", "UBE2V1", "TP53", "PPIA", "UBA52", "ARL8A", "CARD11", 
                  "CFB", "S100A7", "RAB5B", "PSMD14", "SHC1", "DHX9", "CUL1", 
                  "PLD4", "TARP", "CNN2", "SDCBP", "PRDX5", "PDZD2", "PCBP2", 
                  "CD300LB", "PI3", "NLRP1", "ELANE", "ATP6V1F", "SIGLEC16", 
                  "HSPA8", "XRCC6", "TRAPPC1", "KLRC2", "XRCC5", "MUC17", 
                  "ARPC5", "VNN1", "ARPC2", "ARPC3", "BCL2L1", "NFKBIB", 
                  "ARHGAP9", "PIGR", "IFNA4", "HSP90AB1", "DDX3X", "ARPC1B", 
                  "CTSZ", "HMGB1", "TCIRG1", "SLC2A5", "AP2A2", "ORM2", "GM2A", 
                  "AIFM2", "STBD1", "ITGAX", "RAC2", "CASP1", "CASP2", "CTSD", 
                  "HRAS", "CTSB", "ACTR3", "ACTR2", "ATP6V0B", "EGFL7", "DBNL", 
                  "FCER1G", "ANXA2", "GAA", "NFAM1", "NME2", "PPP2R5D", "KCNAB2", 
                  "MUC5AC", "APRT", "DNM2", "VAMP8", "HCK", "PKM", "PSMA2", 
                  "NPC2", "LCK", "PGAM4", "MYH9", "PAFAH1B2", "CD44", "TRIM56", 
                  "RAB7A", "LILRA6", "CD63", "CSTB", "USP14", "MOSPD2", "DEFB123", 
                  "AGPAT2", "USP18", "TSPAN32", "RELA", "CDC42", "IFI16", "PSMB3", 
                  "CXCR2", "POLR2E", "RPS27A", "LAIR1", "POLR2L", "MUC6", 
                  "DYNC1H1", "GSN", "TMEM30A", "PTGES2", "PLA2G2A", "CARD9", 
                  "PDAP1", "LILRB2", "DEFB132", "MUC5B", "TUBB4B", "LRG1", 
                  "PSMC6", "CFHR1", "PSMC1", "CAPZA1", "SAA2", "PKN1", "ACTR10", 
                  "LAMTOR3", "POLR3K")

CD4T_T_Cell_Proliferation <-c("CCDC88B", "IL4I1", "HLA.DMB", "TFRC", "GPNMB", "CD3E", "HLA.G", "LILRB4", "HLA.DRB1", "HLA.DPA1", "SYK", "EBI3", "RPS3", "HMGB1", "AIF1")
CD4T_Antigen_Receptor_signaling <- c("FYB1", "SYK", "HLA.A", "CD3E", "CD3D", "IGHG3", "CD79B", "LAT2", "IGHG4", "CD79A", "IGHG1", "IGKC", "CSK", "HLA.DRB1", "HLA.DPA1")
CD4T_Regulatory_T_Cell_Differentiation <- c("IL4I1", "HLA.DRA", "HLA.G", "LILRB4", "HLA.DRB1", "FANCA")
CD4T_Programmed_Cell_Death <-  c("UNC5A", "STAT3", "APIP", "HMGB1", "H1.2", "LMNB1", "H1.5", "SDCBP", "PSMC6", "PSMB3", "AKT2", "PSMD2", "PSMC1", "PSME2", "MAPT", "TP53", "ELANE", "BCL2L1", "PLEC")
CD4T_Cellular_Senescence <- c("CBX8", "PHC1", "CDKN1B", "H2AC6", "CBX2", "TINF2", "H2BC4", "STAT3", "H2AC19", "HMGA1", "ETS1", "PHC3", "H2BC17", "NFKB1", "TXN2", "H1.2", "CABIN1", "LMNB1", "H1.5", "H4C15", "H3C13", "UBE2S", "TP53", "H3.3B")
CD4T_Cell_Population_Proliferation <- c("BTG2", "VIPR1", "CDKN1B", "TFRC", "PDCD5", "IFI35", "FGF1", 
                                        "ING1", "ETS1", "CCAR1", "CRKL", "RPS4X", "AKT2", "FTH1", 
                                        "PIM2", "AMH", "CER1", "RPS9", "H2AC6", "PRMT1", "RPS6", 
                                        "P3H1", "BTC", "BIRC5", "PFDN1", "SGK1", "TP53", "CTBP2", 
                                        "MAZ", "COX17", "PLG", "TOB2", "AIF1", "EGFR", "GNAI2", "SDCBP", "RPS15A", 
                                        "GPNMB", "DDRGK1", "STX4", "ROMO1", "RPL17", "AMELX", "SLAMF1", "TRPM4", "MCTS1", "ZBTB17", "CGREF1", "XRCC6", 
                                        "NPM1", "NOS3", "STAT3", "OSGIN1", "HMGA1", "SOD2", "SLURP1", "EPOR", "DNAJA3", "AQP11", "CALR", "INHA", "FGFR4", "LGMN", "PF4")
CD4T_T_Cell_Mediated_Immunity  <- c("IL4I1", "CEACAM1", "HLA.G", "LILRB4", "HLA.B", "HLA.C", "HLA.DRA", "HLA.A", "HLA.G", "ULBP1", "HLA.DRB1", "HLA.E")
CD4T_Cytolysis <- c("CFHR1","ROMO1")
CD4T_T_Cell_Proliferation <- c("CCDC88B", "HLA.DMB", "TFRC", "SYK", "EBI3", "RPS3", "HMGB1", "CD3E", "AIF1", "HLA.DPA1")
CD4T_Regulation_Of_Viral_Process <- c("BTBD17", "APOBEC3H", "MX1", "IFIT5", "UBP1", "TRIM27", "TRIM31", "MID2")
CD4T_Positive_Regulation_Of_TNF <- c("IL33", "SYK", "STAT3", "CCL3", "TWIST1", "HMGB1", "AMH", "ORM2", "TMEM106A", "PF4")
CD4T_customized_Tcell_exhaustion <- c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "BTLA", "TIGIT", 
                                      "CD160", "CD244", "ENTPD1", "VSIR")
CD4T_Innate_Immune_System <- c("HSP90AB1", "DDX3X", "NCF1", "UBE2D2", "CTSZ", "WIPF3", "HMGB1", 
                               "SRP14", "SLC2A5", "ORM2", "TXN2", "AIFM2", "PSMD2", "STBD1", 
                               "FTH1", "CTSD", "GYG1", "B2M", "GOLGA7", "CTSB", "SKP1", 
                               "ACTR3", "ATP6V0B", "CHRNB4", "EGFL7", "ATP6V0E1", "SYK", 
                               "ANXA2", "FBXW11", "HLA.B", "NME2", "HLA.C", "PPP2R5D", 
                               "CNPY3", "HLA.A", "HLA.E", "PLPP5", "CEACAM1", "PKM", "PLPP4", 
                               "NPC2", "ELMO1", "PSME2", "PGAM4", "CHI3L1", "UBE2V1", "TP53", 
                               "PPIA", "PAFAH1B2", "ARL8A", "CD44", "CSTB", "SFTPA2", "USP18", 
                               "PLD3", "CDC42", "SDCBP", "PDZD2", "PSMB3", "PCBP2", "POLR2E", 
                               "ATP6V1H", "SLC15A4", "ELANE", "LAIR1", "POLR2L", "HSPA8", 
                               "XRCC6", "TMEM30A", "CARD9", "LYZ", "MUC5B", "TUBB4B", 
                               "NFKB1", "PSMC6", "ARPC2", "CFHR1", "PSMC1", "CAPZA1", 
                               "SAA2", "GRB2", "CD68", "ACTR10", "LGMN", "BCL2L1")
## lingua's gene list
CD4T_PancancerT_NaiveTcell <- c("IL7R", "CCR7", "SELL", "FOXP1", "KLF2", "KLF3", "LEF1", "TCF7", "ACTN1", "BTG1", "BTG2", "TOB1")
CD4T_PancancerT_Activation_Effector_function <- c(
  "FAS", "CD44", "CD69", "CD38", "NKG7", 
  "KLRB1", "KLRD1", "KLRG1", "CX3CR1", "CD300A", 
  "FGFBP2", "ID2", "ID3", "PRDM1", "RUNX3", 
  "TBX21", "ZEB2", "BATF", "NR4A1", "NR4A2", 
  "HOPX", "FOS", "FOSB", "FOSL2", "JUN", 
  "JUNB", "JUND", "STAT1", "STAT3", "EOMES", "AHR"
)
CD4T_PancancerT_Exhaustion <- c(
  "PDCD1", "LAYN", "HAVCR2", "LAG3", "CTLA4", 
  "TIGIT", "TOX", "VSIR", "BTLA", "ENTPD1"
)
CD4T_PancancerT_TCR_Signaling <- c(
  "CALM1", "CALM2", "CALM3", "CD4", "CAST", 
  "CD247", "CD3D", "CD3E", "CD3G", "CSK", 
  "DOK2", "FYN", "LCK", "NFATC1", "NFATC2", 
  "PLEK", "PTPN11", "PTPN13", "PTPN2", "PTPN22", 
  "PTPN4", "PTPN6", "PTPN7", "PTPRC", "PTPRCAP", 
  "S100A10", "S100A11", "S100A4", "S100A6", "ZAP70", 
  "DUSP1", "DUSP2", "DUSP4", "DUSP16", "LAT", 
  "FOS", "FOSB", "FOSL2", "JUN", "JUNB", 
  "JUND", "NR4A1", "NR4A2", "BATF", "IRF1", 
  "SH2D1A", "SH2D2A", "MAP2K3", "MAP3K4", "MAP3K8", 
  "MAP4K1", "NFKB2", "NFKBIA", "NFKBIZ", "REL", "RELB"
)
CD4T_PancancerT_Cytotoxicity <- c(
  "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", 
  "GNLY", "PRF1", "IFNG", "TNF", "SERPINB9", 
  "CTSA", "CTSB", "CTSC", "CTSD", "CTSH", 
  "CTSW", "CST7", "CAPN2", "PLEK"
)
CD4T_PancancerT_Cytokine_CytokineReceptor <- c(
  "CSF1", "CSF2", "ADAM19", "ADAM8", "ADAM12", 
  "CD70", "IL12RB2", "IL17A", "IL17F", "IL10RA", 
  "IL1R1", "IL1R2", "IL21R", "IL21", "IL26", 
  "IL2RA", "IL2RB", "IL2RG", "IL32", "IL6R", 
  "TGFB1", "TGFBR2", "TGFBR3"
)
CD4T_PancancerT_Chemokine_ChemokineReceptor <- c(
  "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", 
  "CXCR3", "CXCR4", "CXCR5", "CXCR6", 
  "CCL3", "CCL4", "CCL5", "CCL20", 
  "CXCL13", "CXCL8", "XCL1"
)
CD4T_PancancerT_Adhesion <- c(
  "ITGA6", "ITGA4", "ITGAE", "ITGAL", "ITGAM", 
  "ITGB1", "ITGB2", "ITGB7", "SELL", "SELPLG", 
  "S1PR1", "ICAM2"
)
CD4T_PancancerT_IFN_Response <- c(
  "STAT1", "STAT3", "MX1", "IRF1", "ISG15", 
  "ISG20", "IFITM1", "IFITM2", "IFITM3", "OAS1", 
  "OAS2", "OASL", "SOCS1", "SOCS3", "TRIM22", 
  "APOL6", "IFNAR2", "IFNGR1", "GBP1", "GBP2", 
  "GBP4", "GBP5", "BST2", "IFI16", "IFI35", 
  "IFI44L", "IFI6", "PARP8", "PARP9"
)
CD4T_PancancerT_Treg_Signature <- c(
  "FOXP3", "IKZF2", "IKZF4", "IL2RA", "ENTPD1", 
  "CCR4", "ICOS", "IL10RA", "TGFB1", "TIGIT", 
  "CTLA4", "LAG3", "HAVCR2", "PDCD1"
)
CD4T_PancancerT_Costimulatory_Molecules <- c(
  "TNFRSF25", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFRSF18", 
  "CD27", "CD28", "CD44", "CD48", "ICOS", 
  "CD2", "SLAMF1", "CD40LG", "CD84"
)
CD4T_PancancerT_Pro_Apoptosis <- c(
  "BAX", "BAG3", "CASP1", "CASP4", "CYCS", "BCL2L11"
)
CD4T_PancancerT_Anti_Apoptosis <- c(
  "BIN1", "BIN2", "BIRC3", "BCL2", "BCL2L1", "MCL1"
)
CD4T_PancancerT_StressResponse <- c("CDK4", "CDK6", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "FOS", "JUN", "MAPK1", "MAPK3", 
                                    "RPS27A", "UBA52", "UBB", "UBC", "MIR4738", "BMI1", "CCNA2", "CDC27", "CDK2", "CDKN1A", 
                                    "CDKN1B", "CEBPB", "MAPK14", "E2F1", "E2F2", "E2F3", "PHC1", "PHC2", "EZH2", "IFNB1", 
                                    "IGFBP7", "IL1A", "IL6", "CXCL8", "MDM2", "MDM4", "MAP3K5", "MOV10", "NFKB1", "MAPK7", 
                                    "MAPK8", "MAPK11", "MAPK9", "MAPK10", "MAP2K3", "MAP2K6", "MAP2K7", "RBBP4", "RBBP7", 
                                    "RELA", "RING1", "RNF2", "RPS6KA1", "RPS6KA2", "RPS6KA3", "MAP2K4", "STAT3", "TFDP1", 
                                    "TFDP2", "TP53", "TXN", "UBE2D1", "UBE2E1", "MAPKAPK3", "CBX4", "MAPKAPK5", "CDC23", 
                                    "EED", "CDC16", "CCNA1", "MAPKAPK2", "MAP4K4", "ANAPC10", "EHMT2", "UBE2C", "SCMH1", 
                                    "TNIK", "TNRC6B", "KDM6B", "CBX6", "SUZ12", "ANAPC15", "AGO1", "TNRC6A", "ANAPC2", 
                                    "ANAPC4", "MINK1", "FZR1", "ANAPC5", "ANAPC7", "ANAPC11", "CBX8", "TNRC6C", "ANAPC1", 
                                    "EHMT1", "PHC3", "CBX2", "ANAPC16", "AGO3", "AGO4", "CDC26", "MIR3605", "MIR6755", 
                                    "ATM", "CCNE1", "ERF", "ETS1", "ETS2", "HMGA1", "ID1", "LMNB1", "MRE11", "NBN", "RB1", 
                                    "SP1", "TERF1", "TERF2", "HIRA", "HMGA2", "CCNE2", "RAD50", "KAT5", "CABIN1", "ASF1A", 
                                    "POT1", "TINF2", "UBN1", "TERF2IP", "EP400", "ACD", "NUDT2", "AQP8", "ARNT", "ATOX1", 
                                    "ATP7A", "ATR", "BAG1", "CA9", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAT", "CREBBP", 
                                    "CYBA", "CYBB", "EEF1A1", "EP300", "EPAS1", "EPO", "FKBP4", "MTOR", "GPX1", "GPX2", 
                                    "GPX3", "GPX5", "GPX7", "GSK3B", "GSR", "GSTP1", "HIF1A", "HSBP1", "HSF1", "HSPA1A", 
                                    "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA7", "HSPA8", "HSPA9", 
                                    "HSP90AA1", "HSP90AB1", "DNAJB1", "NCF2", "NCF4", "NUP88", "NUP98", "P4HB", "PRDX1", 
                                    "PIK3C3", "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG1", "PSMA1", "PSMA2", "PSMA3", 
                                    "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", 
                                    "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10", "PSMC1", "PSMC2", "PSMC3", "PSMC4", 
                                    "PSMC5", "PSMC6", "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD5", "PSMD7", "PSMD8", 
                                    "PSMD9", "PSMD10", "PSMD11", "PSMD12", "PSMD13", "PSME1", "PSME2", "RANBP2", "RHEB", 
                                    "RPA1", "RPA2", "RPA3", "SOD1", "SOD2", "SOD3", "ST13", "HSPA13", "ELOC", "ELOB", 
                                    "PRDX2", "TPR", "TSC1", "TSC2", "DNAJC7", "TXNRD1", "UBE2D2", "UBE2D3", "UVRAG", 
                                    "VCP", "VEGFA", "VHL", "YWHAE", "SEM1", "NUP214", "AAAS", "ULK1", "CUL2", "RAE1", 
                                    "LAMTOR3", "DYNLL1", "BECN1", "MTMR3", "LIMD1", "ATG12", "ATG5", "PSMF1", "BAG5", 
                                    "BAG4", "BAG3", "BAG2", "PRDX6", "NUP155", "NUP93", "ATG13", "NUP58", "RB1CC1", 
                                    "PSMD6", "POM121", "NUP153", "CCS", "RBX1", "HDAC6", "DNAJB6", "PSME3", "PSMD14", 
                                    "RRAGB", "CITED2", "ATG7", "LAMTOR5", "TXNRD2", "RRAGA", "PTGES3", "NUP50", "HSPH1", 
                                    "PRDX3", "WDR45", "GABARAP", "GABARAPL2", "HSPA4L", "ATG14", "NUP205", "ATG4B", 
                                    "PSME4", "NUP210", "NUP160", "SIRT1", "NUP188", "NUP62", "GABARAPL1", "GABARAPL3", 
                                    "PRDX5", "TXN2", "CHMP2B", "HIGD1A", "WIPI2", "DNAJC2", "CHMP2A", "LAMTOR2", 
                                    "CHMP4A", "ERO1A", "PIK3R4", "NOX4", "HSPA14", "PRKAG2", "HIKESHI", "CHMP3", 
                                    "NUP54", "PRKAG3", "CYCS", "EGLN1", "LAMTOR1", "ATG16L1", "WIPI1", "AMBRA1", 
                                    "HIF1AN", "NDC1", "NUP133", "WDR45B", "NUP107", "RPTOR", "CCAR2", "RRAGD", 
                                    "ATG101", "RRAGC", "MLST8", "HIF3A", "MTMR14", "ATG3", "NUP37", "ATG9A", "NOX5", 
                                    "CHMP6", "NUP85", "MAP1LC3B", "SEH1L", "ATG10", "AKT1S1", "MAP1LC3A", "ATG4C", 
                                    "AJUBA", "ATG4D", "RPS19BP1", "CHMP7", "CHMP4C", "EGLN2", "EGLN3", "ATG4A", 
                                    "HSPA12B", "PSMB11", "WTIP", "CHMP4B", "NUP35", "DYNLL2", "PSMA8", "GPX6", 
                                    "HSPA12A", "ATG9B", "NUP43", "LAMTOR4", "MAP1LC3C", "GPX8", "NCF1", "POM121C", 
                                    "SOD2-OT1", "MIR1281", "MIR7703")
# Function to remove MHC molecules (genes starting with "HLA.")
remove_mhc <- function(gene_vector) {
  gene_vector[!grepl("^HLA\\.", gene_vector)]
}

# Apply the function to each vector
CD4T_T_Cell_Proliferation <- remove_mhc(CD4T_T_Cell_Proliferation)
CD4T_Antigen_Receptor_signaling <- remove_mhc(CD4T_Antigen_Receptor_signaling)
CD4T_Regulatory_T_Cell_Differentiation <- remove_mhc(CD4T_Regulatory_T_Cell_Differentiation)
CD4T_Programmed_Cell_Death <- remove_mhc(CD4T_Programmed_Cell_Death)
CD4T_Cellular_Senescence <- remove_mhc(CD4T_Cellular_Senescence)
CD4T_Cell_Population_Proliferation <- remove_mhc(CD4T_Cell_Population_Proliferation)
CD4T_T_Cell_Mediated_Immunity <- remove_mhc(CD4T_T_Cell_Mediated_Immunity)
CD4T_Cytolysis <- remove_mhc(CD4T_Cytolysis)
CD4T_Regulation_Of_Viral_Process <- remove_mhc(CD4T_Regulation_Of_Viral_Process)
CD4T_Positive_Regulation_Of_TNF <- remove_mhc(CD4T_Positive_Regulation_Of_TNF)
functionalenrichment <- list(
  Macro_merged_translation = Macro_merged_translation,
  Macro_Macrophage_differentiation = Macro_Macrophage_differentiation,
  Macro_Gene_expression = Macro_Gene_expression,
  Macro_TOR_signaling = Macro_TOR_signaling,
  Macro_Regulation_T_Cell_Differentiation = Macro_Regulation_T_Cell_Differentiation,
  Macro_MHC_class_II = Macro_MHC_class_II,
  Macro_NF_kappaB = Macro_NF_kappaB,
  Macro_Response_To_Hypoxia = Macro_Response_To_Hypoxia,
  Macro_Negative_Regulation_Of_Inflammatory_Response = Macro_Negative_Regulation_Of_Inflammatory_Response,
  Macro_Signal_Transduction_In_Absence_Of_Ligand = Macro_Signal_Transduction_In_Absence_Of_Ligand,
  # macro in tumor
  Macro_Negative_Regulation_Of_Tumor_Necrosis_Factor_Production = Macro_Negative_Regulation_Of_Tumor_Necrosis_Factor_Production,
  Macro_Positive_Regulation_Of_T_Cell_Activation = Macro_Positive_Regulation_Of_T_Cell_Activation,
  Macro_MHC_Class_II_Protein_Complex_Assembly =Macro_MHC_Class_II_Protein_Complex_Assembly,
  # Tumor
  Tumor_Regulation_Of_Macrophage_Proliferation = Tumor_Regulation_Of_Macrophage_Proliferation,
  Tumor_Regulation_Of_Macrophage_Chemotaxis = Tumor_Regulation_Of_Macrophage_Chemotaxis, 
  Tumor_Antigen_Processing_And_Presentation_Of_Exogenous_Peptide_Antigen_Via_MHC_Class_II = Tumor_Antigen_Processing_And_Presentation_Of_Exogenous_Peptide_Antigen_Via_MHC_Class_II,
  Tumor_Peptide_Antigen_Assembly_With_MHC_Protein_Complex = Tumor_Peptide_Antigen_Assembly_With_MHC_Protein_Complex,
  Tumor_Negative_Regulation_Of_Response_To_Type_II_Interferon = Tumor_Negative_Regulation_Of_Response_To_Type_II_Interferon,
  Tumor_IFNG_cascade = Tumor_IFNG_cascade,
  Tumor_IRF_family = Tumor_IRF_family,
  Tumor_Positive_Regulation_Of_Leukocyte_Mediated_Cytotoxicity = Tumor_Positive_Regulation_Of_Leukocyte_Mediated_Cytotoxicity,
  Tumor_Innate_Immune_System_R_HSA_168249 =Tumor_Innate_Immune_System_R_HSA_168249,
  #
  CD4T_T_Cell_Proliferation = CD4T_T_Cell_Proliferation,
  CD4T_Antigen_Receptor_signaling = CD4T_Antigen_Receptor_signaling,
  CD4T_Regulatory_T_Cell_Differentiation = CD4T_Regulatory_T_Cell_Differentiation,
  CD4T_Programmed_Cell_Death = CD4T_Programmed_Cell_Death,
  CD4T_Cellular_Senescence = CD4T_Cellular_Senescence,
  CD4T_Cell_Population_Proliferation = CD4T_Cell_Population_Proliferation,
  CD4T_T_Cell_Mediated_Immunity = CD4T_T_Cell_Mediated_Immunity,
  CD4T_Cytolysis = CD4T_Cytolysis,
  CD4T_Regulation_Of_Viral_Process = CD4T_Regulation_Of_Viral_Process,
  CD4T_Positive_Regulation_Of_TNF = CD4T_Positive_Regulation_Of_TNF,
  # CD4T_PancancerT_NaiveTcell = CD4T_PancancerT_NaiveTcell,
  # CD4T_PancancerT_Activation_Effector_function = CD4T_PancancerT_Activation_Effector_function,
  # CD4T_PancancerT_Exhaustion = CD4T_PancancerT_Exhaustion,
  # CD4T_PancancerT_TCR_Signaling = CD4T_PancancerT_TCR_Signaling,
  # CD4T_PancancerT_Cytotoxicity = CD4T_PancancerT_Cytotoxicity,
  # CD4T_PancancerT_Cytokine_CytokineReceptor = CD4T_PancancerT_Cytokine_CytokineReceptor,
  # CD4T_PancancerT_Chemokine_ChemokineReceptor = CD4T_PancancerT_Chemokine_ChemokineReceptor,
  CD4T_PancancerT_Adhesion = CD4T_PancancerT_Adhesion,
  # CD4T_PancancerT_IFN_Response = CD4T_PancancerT_IFN_Response,
  # CD4T_PancancerT_Treg_Signature = CD4T_PancancerT_Treg_Signature,
  # CD4T_PancancerT_Costimulatory_Molecules = CD4T_PancancerT_Costimulatory_Molecules,
  # CD4T_PancancerT_Pro_Apoptosis = CD4T_PancancerT_Pro_Apoptosis,
  # CD4T_PancancerT_Anti_Apoptosis = CD4T_PancancerT_Anti_Apoptosis,
  # CD4T_PancancerT_StressResponse = CD4T_PancancerT_StressResponse
  CD4T_customized_Tcell_exhaustion = CD4T_customized_Tcell_exhaustion,
  CD4T_Innate_Immune_System = CD4T_Innate_Immune_System
)

# table(annotation_merged$coreName)
# annotation_merged %>% 
#   filter(coreName == "DFCI_14.1") %>% 
#   ggplot(aes(x = X_cent, y = Y_cent, color = LMP1)) + 
#   geom_point() +
#   scale_color_gradient(low = "white", high = "red")
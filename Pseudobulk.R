# Load Packages
library(readxl)
library(gtsummary)
library(labelled)
library(Hmisc)
library(tidyselect)
library(tidyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(ComplexHeatmap)
library(lubridate)
library(data.table)
library(circlize)
library(dendextend)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(rstatix)
library(corrplot)
library(writexl)
library(plyr)
library(digest)
library(BiocManager)
library(ggVennDiagram)
library(flextable)
library(officer)
library(nationalparkcolors)

# Upload Data
single_pt <- read_excel("/Users/jroot/Library/CloudStorage/OneDrive-InsideMDAnderson/LabMembers/JessicaRoot/Isoplexis Analysis/Data/New_Metadata/metadata_single_patient.xlsx", sheet = 1)
df1 <- read_excel("/Users/jroot/Library/CloudStorage/OneDrive-InsideMDAnderson/LabMembers/JessicaRoot/Isoplexis Analysis/Data/New_Metadata/combined_metadata.xlsx", sheet = 1)

# Subset the Data
# Create a vector of Donor IDs from the metadata for the patient cohort
cohort <- c("PB1B","PB2B","PB3B","PB4B","PB5B","PB6B","PB7B","PB8B","PB9B","PB10B",
            "PB11B","PB12B","PB13B","PB14B","PB15B","PB16B","PB17B","PB18B","PB19B","PB20B","BM21B")

# Filter metadata 
new_metadata <- filter(df1, DonorID %in% cohort)

# Subset filtered metadata for only one cell subset - each patient has two cell subsets per sample
new_metadata <- new_metadata %>% 
  filter(CellSubset == "CD4+") 
new_metadata <- new_metadata %>% 
  mutate(Dx = recode(Dx, Mega.="Megakaryocytic"))
new_metadata <- new_metadata %>% 
  mutate(Dx = recode(Dx, AMOL="AML"))
new_metadata <- new_metadata %>% 
  mutate(Dx = recode(Dx, AMML="AML"))

# Venn Diagram
# create vectors for each cohort
baseline_bm <- c("PB1B","PB2B","PB3B","PB4B","PB5B","PB6B","PB7B","PB8B","PB9B","PB10B",
                 "PB16B","PB17B","PB18B","PB19B","PB20B","BM21B")
baseline_pb <- c("PB1B","PB2B","PB3B","PB4B","PB5B","PB6B","PB7B","PB8B","PB9B","PB10B",
                 "PB11B","PB12B","PB13B","PB14B","PB15B","PB16B","PB17B","PB18B","PB19B","PB20B")
paired_bm <- c("PB16B","PB17B","PB18B","PB19B","PB20B","BM21B")
paired_pb <- c("PB11B","PB12B","PB13B","PB14B","PB15B","PB16B","PB17B","PB18B","PB19B","PB20B")

x <- list(baseline_bm, baseline_pb, paired_bm, paired_pb)   # List the vectors

venn_cohort <- ggVennDiagram(x, category.names = c("Baseline BM",     # Plot venn diagram
                                                   "Baseline PB",
                                                   "Paired BM",
                                                   "Paired PB"),
                             label_size = 4,
                             label_alpha = 0,
                             label_color = "Black") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_fill_distiller(palette = "Blues") +
  scale_color_manual(values = c("Darkgrey","Darkgrey","Darkgrey","Darkgrey"))
plot(venn_cohort)

# Table Summary
single_pt <- single_pt %>% 
  mutate(Dx = recode(Dx, "Mega."="Megakaryocytic AML"))
single_pt <- single_pt %>% 
  mutate(Dx = recode(Dx, "AMOL"="AML"))
single_pt <- single_pt %>% 
  mutate(Dx = recode(Dx, "AMML"="AML"))
single_pt <- single_pt %>% 
  mutate(Resp = recode(Resp, "CR"="CR/CRi/HI"))

single_pt <- single_pt %>% select("Age","Sex","Race","Dx","sec","Karyotype_new","PRIOR HMA","PRIOR SCT","WBC","Blast","Tot Co","Resp",#"CR/CRi/PR vs NR",
                                  "Time to Resp-m")

# Adjust class of variables if necessary - this will determine how your variables are displayed in the table
single_pt$`PRIOR HMA` <- as.character(single_pt$`PRIOR HMA`)


var_label(single_pt) <- list(Sex="Gender", Dx="Acute Leukemia Subtype", `sec`="Secondary AML", `Karyotype_new`="Cytogenetic Group",
                             `PRIOR HMA`="Prior HMA", `PRIOR SCT` = "Prior SCT", Blast = "PB Blast %", `Tot Co` = "Total Aza/Nivo Cycles",
                             Resp="Response",`Time to Resp-m` = "Time to Response (months)")

single_pt %>%
  tbl_summary(statistic = all_continuous() ~ c("{min}, {median}, {max}")) %>%
  bold_labels() %>%
  modify_caption("Table 1. Clinical Characteristics") %>%
  as_flex_table() %>%
  set_table_properties(layout = "autofit") %>% 
  autofit(add_w = 0, add_h = 0) %>%
  theme_vanilla()

# Oncoprint
single_pt <- read_excel("/Users/jroot/Library/CloudStorage/OneDrive-InsideMDAnderson/LabMembers/JessicaRoot/Isoplexis Analysis/Data/New_Metadata/metadata_single_patient.xlsx", sheet = 1)

subset_metadata <- filter(single_pt, DonorID %in% cohort) 
names(subset_metadata)[names(subset_metadata) == 'Start'] <- "Start_date"
names(subset_metadata)[names(subset_metadata) == 'Off'] <- "Response_date"
names(subset_metadata)[names(subset_metadata) == 'Dx Date'] <- "Dx_date"
names(subset_metadata)[names(subset_metadata) == 'Died Day?'] <- "died_date"
names(subset_metadata)[names(subset_metadata) == 'StatusR'] <- "Response"
names(subset_metadata)[names(subset_metadata) == 'Dx'] <- "Diagnosis"
names(subset_metadata)[names(subset_metadata) == 'PRIOR SCT'] <- "Prior SCT"
names(subset_metadata)[names(subset_metadata) == 'PRIOR HMA'] <- "Prior HMA"

subset_metadata <- subset_metadata %>%                                         
  dplyr::select(DonorID, everything()) 
subset_metadata <- data.frame(subset_metadata)

names(subset_metadata)[names(subset_metadata) == 'Baseline_BM'] <- "Baseline BM"
names(subset_metadata)[names(subset_metadata) == 'Baseline_PB'] <- "Baseline PB"
names(subset_metadata)[names(subset_metadata) == 'Paired_BM'] <- "Post BM"
names(subset_metadata)[names(subset_metadata) == 'Paired_PB'] <- "Post PB"
names(subset_metadata)[names(subset_metadata) == 'Prior.SCT'] <- "Prior SCT"
names(subset_metadata)[names(subset_metadata) == 'Prior.HMA'] <- "Prior HMA"
names(subset_metadata)[names(subset_metadata) == 'sec'] <- "Secondary AML"
names(subset_metadata)[names(subset_metadata) == 'Resp'] <- "Response to Aza/Nivo"
names(subset_metadata)[names(subset_metadata) == 'Karyotype_new'] <- "Karyotype"
names(subset_metadata)[names(subset_metadata) == 'Del5.5q'] <- "Del5/5q"
names(subset_metadata)[names(subset_metadata) == 'Del7.7q'] <- "Del7/7q"
names(subset_metadata)[names(subset_metadata) == 'Del17.17p'] <- "Del17/17p"

subset_metadata <- subset_metadata %>% 
  mutate(Diagnosis = recode(Diagnosis, Mega.="Megakaryocytic AML"))
subset_metadata <- subset_metadata %>% 
  mutate(Diagnosis = recode(Diagnosis, AMOL="AML"))
subset_metadata <- subset_metadata %>% 
  mutate(Diagnosis = recode(Diagnosis, AMML="AML"))
subset_metadata <- subset_metadata %>% 
  mutate(`Response to Aza/Nivo` = recode(`Response to Aza/Nivo`, CR="CR/CRi/HI"))

subset_metadata <- subset_metadata %>% 
  mutate(DonorID = recode(DonorID, PB1B="PT1", PB2B="PT2", PB3B="PT3", PB4B="PT4", PB5B="PT5", PB6B="PT6", PB7B="PT7",
                          PB8B="PT8", PB9B="PT9", PB10B="PT10",PB11B="PT11", PB12B="PT12", PB13B="PT13", PB14B="PT14",
                          PB15B="PT15", PB16B="PT16", PB17B="PT17", PB18B="PT18", PB19B="PT19", PB20B="PT20",
                          BM21B="PT21"))
groups <- c("Baseline_BM", "Baseline_PB", "Post_BM", "Post_PB")

# Oncoprint separated by response
response_group <- c("PT6","PT7","PT8","PT9","PT10","PT11","PT12","PT13","PT14","PT15",
                    "PT1","PT2","PT3","PT4","PT5","PT16","PT17","PT18","PT19","PT20","PT21")
subset_metadata <- subset_metadata %>% 
  slice(match(response_group, DonorID))

rownames(subset_metadata) <- subset_metadata$DonorID

subset_metadata["Secondary AML"][subset_metadata["Secondary AML"] == "2"] <- "Yes"
subset_metadata["Secondary AML"][subset_metadata["Secondary AML"] == "1"] <- "No"

subset_metadata["Karyotype"][subset_metadata["Karyotype"] == "Nondiploid"] <- "Other Intermediate"

subset_metadata["Prior SCT"][subset_metadata["Prior SCT"] == "TRUE"] <- "T"
subset_metadata["Prior SCT"][subset_metadata["Prior SCT"] == "FALSE"] <- "F"
subset_metadata["Prior HMA"][subset_metadata["Prior HMA"] == "TRUE"] <- "T"
subset_metadata["Prior HMA"][subset_metadata["Prior HMA"] == "FALSE"] <- "F"

# Prepare data frame
cyto_4 <- subset_metadata[,c(10:13,65,67,69,25,15,18:20,125:141)]
cyto_4 <- cyto_4[rev(order(cyto_4[,2],cyto_4[,1], cyto_4[,4], cyto_4[,3],cyto_4[,9])),]
cyto_4 <- t(cyto_4)

# Prepare legends
lgd1 = Legend(title = "Sample Present",
              labels = c("True", "False"), 
              legend_gp = gpar(fill = c("steelblue2", "lightgray")))
lgd2= Legend(labels=c("Yes", "No"), 
             legend_gp = gpar(fill=c("#8CBEB1", "lightgrey")), 
             title="Prior Therapy")
lgd4= Legend(labels=c("Yes", "No"), 
             legend_gp = gpar(fill=c("seagreen", "lightgray")), 
             title="Secondary AML")
lgd5= Legend(labels=c("CR/CRi/HI", "NR"), 
             legend_gp = gpar(fill=c("#E79498", "lightgrey")), 
             title="Response to Aza/Nivo")
lgd6= Legend(labels=c("Complex", "Other Intermediate", "Diploid", "Not Done"),
             legend_gp = gpar(fill=c("royalblue3","#2E8289","#91D5DE", "darkgrey")),
             title="Karyotype")
lgd7 = Legend(labels = c("Yes", "No"), 
              legend_gp = gpar(fill = c("plum3","lightgrey")),
              title = "Mutation")
lgd8 = Legend(labels = c("Yes", "No", "Not Tested"), 
              legend_gp = gpar(fill = c("#CEA347","lightgrey","darkgrey")),
              title = "Chromosome Deletion")

pd = packLegend(lgd1, lgd2, lgd4, lgd5, lgd6, lgd7, lgd8, direction = "horizontal", gap = unit("1", "cm"))

# Draw Oncoprint
Heatmap(cyto_4, col = c("TRUE" = "steelblue2", "FALSE"="lightgray",
                        "T"="#8CBEB1", "F"="lightgrey",
                        "Yes" = "seagreen", "No" = "lightgray",
                        "CR/CRi/HI"="#E79498","NR"="lightgrey",
                        "Complex"="royalblue3","Other Intermediate"="#2E8289", "Diploid"="#91D5DE","Not Done"="darkgrey",
                        "Del5/5q"="#CEA347", "Not Done"="darkgrey", "No"="lightgrey",
                        "Del7/7q"="#CEA347", "Not Done"="darkgrey", "No"="lightgrey",
                        "Del17/17p"="#CEA347", "Not Done"="darkgrey", "No"="lightgrey",
                        "0"="lightgray", "1"="plum3"),
        column_title = "", row_title=" ", row_split=c(rep("A", 4), rep("B",2), rep("C",2), rep("E",4), 
                                                      rep("F",17)),
        row_gap = unit(1, "mm"),
        width = ncol(cyto_4)*unit(6, "mm"), height = nrow(cyto_4)*unit(6, "mm"), border=FALSE, 
        rect_gp = gpar(col = "white", lwd = 1),
        show_heatmap_legend = FALSE)

draw(pd, just=c("center","bottom"), y=unit(0.5, "cm"))

# Heatmap
# Vector of patients 1-20: all have peripheral blood samples obtained at baseline 
pb_baseline_cohort <- c("PB1B","PB2B","PB3B","PB4B","PB5B","PB6B","PB7B","PB8B","PB9B","PB10B",
                        "PB11B","PB12B","PB13B","PB14B","PB15B","PB16B","PB17B","PB18B","PB19B","PB20B")

#Subset data
baseline_blood <- filter(df1, DonorID %in% pb_baseline_cohort) # subset for cohort as above
baseline_blood <- baseline_blood %>% 
  mutate(DonorID = recode(DonorID, PB1B="PT1", PB2B="PT2", PB3B="PT3", PB4B="PT4", PB5B="PT5",
                          PB6B="PT6", PB7B="PT7", PB8B="PT8", PB9B="PT9", PB10B="PT10",
                          PB11B="PT11", PB12B="PT12", PB13B="PT13", PB14B="PT14", PB15B="PT15",
                          PB16B="PT16", PB17B="PT17", PB18B="PT18", PB19B="PT19", PB20B="PT20"))

baseline_blood <- baseline_blood %>% arrange(CellSubset)
baseline_blood["Karyotype_new"][baseline_blood["Karyotype_new"] == "Nondiploid"] <- "Other Intermediate"

mat1 = data.matrix(baseline_blood[, 35:66]) # create a matrix with dataframe containing only cytokine values
mat_scaled1 = t(apply(mat1, 2, scale))      # scale cytokine expressions, 1 indicates rows -> t(mat), so use 2 (this will across be rows)
mat_scaled1[is.na(mat_scaled1)] <- 0
colfun4 = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))        # create color function

# create vectors
patients <- c(baseline_blood$DonorID)
cytokines <- c("CCL11", "GMCSF", "GranzymeB", "IFNg", "IL10", "IL12", "IL13", "IL15", 
               "IL17A", "IL17F", "IL1b", "IL2", "IL21", "IL22", "IL4", "IL5", "IL6",
               "IL7", "IL8", "IL9", "IP10", "MCP1", "MCP4", "MIP1a", "MIP1b", "Perforin",
               "RANTES", "sCD137", "sCD40L", "TGFb1", "TNFa", "TNFb")
ckgroups <- c("Chemoattractive","Stimulatory","Effector","Effector","Regulatory","Stimulatory","Regulatory","Stimulatory",
              "Inflammatory","Inflammatory","Inflammatory","Stimulatory","Stimulatory","Regulatory","Regulatory","Stimulatory",
              "Inflammatory","Stimulatory","Stimulatory","Stimulatory","Chemoattractive","Inflammatory","Inflammatory","Effector",
              "Chemoattractive","Effector","Chemoattractive","Regulatory","Regulatory","Regulatory","Effector","Effector")


topannodf = data.frame(cellsubset = baseline_blood$CellSubset,
                       Del77q = baseline_blood$`Del7/7q`,
                       Karyotype = baseline_blood$Karyotype_new,
                       PriorHMA = baseline_blood$`PRIOR HMA`,
                       SecondaryAML = baseline_blood$sec,
                       Response = baseline_blood$StatusR)

rightannodf = data.frame(groups = ckgroups)
rightanno = rowAnnotation(
  df = rightannodf,
  col = list(groups = c("Chemoattractive" = "darkslategrey", "Stimulatory" = "deepskyblue", "Effector" = "darkseagreen3", 
                        "Regulatory" = "gold2", "Inflammatory" = "indianred2")))
topanno = HeatmapAnnotation(
  df = topannodf,
  annotation_name_gp = gpar(fontsize = 10),
  col = list(cellsubset = c("CD4+"="#00539CFF", "CD8+"="#EEA47FFF"),
             Del55q = c("Del5/5q"="#CEA347", "No"="lightgray", "Not Done"="darkgrey"),
             Del77q = c("Del7/7q"="#CEA347", "No"="lightgrey", "Not Done"="darkgrey"),
             Del1717p = c("Del17/17p"="#CEA347", "No"="lightgrey", "Not Done"="darkgrey"),
             Karyotype = c("Complex"="royalblue3", "Diploid"="#91D5DE", "Other Intermediate"="#2E8289", "Not Done"="darkgrey"),
             PriorSCT = c("TRUE"="#8CBEB1", "FALSE"="lightgrey"),
             PriorHMA = c("TRUE"="#8CBEB1", "FALSE"="lightgrey"),
             SecondaryAML = c("1"="lightgrey", "2"="seagreen"),
             Response = c("NR" = "lightgrey", "R" = "#6BBAE5"),
             Sex = c("Male" = "paleturquoise", "Female" = "plum3"),
             Race = c("White"="cyan3", "Black"="coral3", "Latinx"="goldenrod4"),
             TreatmentResponse = c("CR"="#E79498", "NR"="lightgrey"),
             Age = colorRamp2(c(0,100), c("ivory","red")),
             OS = colorRamp2(c(0,30), c("ivory","red")),
             WBC = colorRamp2(c(0,50), c("ivory","red")),
             PLT = colorRamp2(c(0,200), c("ivory","red")),
             Blast = colorRamp2(c(0,75), c("ivory","red"))))

leftannodf = data.frame(groups = ckgroups)
rightanno = rowAnnotation(
  df = rightannodf,
  col = list(groups = c("Chemoattractive" = "darkslategrey", "Stimulatory" = "deepskyblue", 
                        "Effector" = "darkseagreen3", "Regulatory" = "gold2", "Inflammatory" = "indianred2")))

bloodmap = Heatmap((mat_scaled1), name = "Expression",
                   clustering_method_rows = "average",     # cluster rows by avg distance between points
                   colfun4,
                   top_annotation = topanno,
                   right_annotation = rightanno,
                   column_split = baseline_blood$CellSubset,
                   row_order = cytokines,
                   row_split = ckgroups,
                   cluster_column_slices = F,
                   cluster_row_slices = F,
                   column_title = "Peripheral Blood CD4 v. CD8 Baseline Cytokine Dynamics",
                   row_labels = cytokines,
                   column_labels = patients,
                   show_parent_dend_line = FALSE,
                   show_column_dend = T,
                   show_row_dend = F,
                   row_names_gp = grid::gpar(fontsize = 12),
                   column_names_gp = grid::gpar(fontsize = 12),
                   column_names_rot = 90)
draw(bloodmap)

# Boxplots 
BLblood <- baseline_blood %>% select(CCL11:TNFb)
BLblood <- names(BLblood)

for (i in BLblood) {
  sel <- c("PatientID", "CellSubset", "Time", "StatusR", "Sample")
  sel2 <- c(sel, i)
  BL_blood_df <- baseline_blood %>% select(all_of(sel2))
  BL_blood_df$CellSubset <- factor(BL_blood_df$CellSubset, levels = c("CD4+","CD8+"))
  
  blood_plot2 <- ggpaired(BL_blood_df, x = "CellSubset", y = i, 
                          palette = c("#8AAAE5","#EEA47FFF"), fill = "CellSubset",
                          line.color = "black", line.size = 0.1, id = "PatientID") +
    geom_point(size = 2) +
    labs(title = paste("PB CD4 vs. CD8", i), x = NULL, y = "Value") +
    theme(legend.position = "right") 
  
  frm <- reformulate("CellSubset", colnames(BL_blood_df)[6])
  test_bl_blood <- t_test(frm, data = BL_blood_df, paired = TRUE, p.adjust.method = "fdr") %>% 
    adjust_pvalue() %>% 
    add_significance() %>% 
    add_xy_position(x = "CellSubset")
  blood_plot2 <- blood_plot2 + stat_pvalue_manual(test_bl_blood, label = "p.adj.signif", size = 4, y.position = "y.position")
  
  print(blood_plot2)
  }


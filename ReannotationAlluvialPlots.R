#### Create alluvial plots to show changes ####
install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)

#### TCGA changes ####

dim(clean.WHO.CNS5_v2)
# [1] 1099   35


## for TCGA compare TCGA.Histology vs WHO_CNS5_histology
table(clean.WHO.CNS5_v2$TCGA.Histology)
table(clean.WHO.CNS5_v2$WHO_CNS5_histology)

TCGA.hist <- as.data.frame(table(clean.WHO.CNS5_v2$TCGA.Histology, clean.WHO.CNS5_v2$WHO_CNS5_histology))

ggplot(data = TCGA.hist, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) +
  geom_stratum() +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            size=7) +
  scale_x_discrete(limits = c("WHO_CNS4", "WHO_CNS5"), expand = c(0.15, 0.05),
                   labels = c("WHO_CNS4", "WHO_CNS5")) +
  theme_minimal() +  # Remove gray background
  theme(axis.text.y = element_blank(),  # Remove Y-axis labels
        axis.title.y = element_blank(),  # Remove Y-axis title
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(size = 30)  # Set X-axis title font size
  ) 

#ggsave("OP/Reanno_TCGAhist.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 600)
ggsave("OP/Reanno_TCGAhist_v2.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 300)


## for TCGA compare TCGA.Grade vs WHO_CNS5_grade
TCGA.gr <- as.data.frame(table(clean.WHO.CNS5_v2$TCGA.Grade, clean.WHO.CNS5_v2$WHO_CNS5_grade))

ggplot(data = TCGA.gr, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) +
  geom_stratum() +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            size = 7) +
  scale_x_discrete(limits = c("WHO_CNS4", "WHO_CNS5"), expand = c(0.15, 0.05),
                   labels = c("WHO_CNS4", "WHO_CNS5")) +
  theme_minimal() +  # Remove gray background
  theme(axis.text.y = element_blank(),  # Remove Y-axis labels
        axis.title.y = element_blank(),  # Remove Y-axis title
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(size = 30)  # Set X-axis title font size
  ) 

#ggsave("OP/Reanno_TCGAgr.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 600)
ggsave("OP/Reanno_TCGAgr_v2.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 300)

## for TCGA compare IDH.status vs WHO_CNS5_IDHstatus
TCGA.idh <- as.data.frame(table(clean.WHO.CNS5_v2$IDH.status, clean.WHO.CNS5_v2$WHO_CNS5_IDHstatus))

ggplot(data = TCGA.idh, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) +
  geom_stratum() +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            size=7) +
  scale_x_discrete(limits = c("WHO_CNS4", "WHO_CNS5"), expand = c(0.15, 0.05),
                   labels = c("WHO_CNS4", "WHO_CNS5")) +
  theme_minimal() +  # Remove gray background
  theme(axis.text.y = element_blank(),  # Remove Y-axis labels
        axis.title.y = element_blank(),  # Remove Y-axis title
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(size = 30)  # Set X-axis title font size
  ) 

ggsave("OP/Reanno_TCGAidh_v2.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 600)

#### GLASS changes ####
dim(clean.GLASS.CNS5.TP)
# [1] 330  53

## for GLASS compare histology vs WHO_CNS5_histology
GLASS.hist <- as.data.frame(table(clean.GLASS.CNS5.TP$histology, clean.GLASS.CNS5.TP$WHO_CNS5_histology))

ggplot(data = GLASS.hist, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) +
  geom_stratum() +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            size=7) +
  scale_x_discrete(limits = c("WHO_CNS4", "WHO_CNS5"), expand = c(0.15, 0.05),
                   labels = c("WHO_CNS4", "WHO_CNS5")) +
  theme_minimal() +  # Remove gray background
  theme(axis.text.y = element_blank(),  # Remove Y-axis labels
        axis.title.y = element_blank(),  # Remove Y-axis title
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(size = 30)  # Set X-axis title font size
  ) 

ggsave("OP/Reanno_GLASShist_v2.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 600)

## for GLASS compare grade vs WHO_CNS5_grade
GLASS.gr <- as.data.frame(table(clean.GLASS.CNS5.TP$grade, clean.GLASS.CNS5.TP$WHO_CNS5_grade))

ggplot(data = GLASS.gr, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) +
  geom_stratum() +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            size=7) +
  scale_x_discrete(limits = c("WHO_CNS4", "WHO_CNS5"), expand = c(0.15, 0.05),
                   labels = c("WHO_CNS4", "WHO_CNS5")) +
  theme_minimal() +  # Remove gray background
  theme(axis.text.y = element_blank(),  # Remove Y-axis labels
        axis.title.y = element_blank(),  # Remove Y-axis title
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(size = 30)  # Set X-axis title font size
  ) 

ggsave("OP/Reanno_GLASSgr_v2.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 600)

## for GLASS compare idh_status vs WHO_CNS5_IDH.status
GLASS.idh <- as.data.frame(table(clean.GLASS.CNS5.TP$IDH.status, clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status))

ggplot(data = GLASS.idh, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) +
  geom_stratum() +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            size=7) +
  scale_x_discrete(limits = c("WHO_CNS4", "WHO_CNS5"), expand = c(0.15, 0.05),
                   labels = c("WHO_CNS4", "WHO_CNS5")) +
  theme_minimal() +  # Remove gray background
  theme(axis.text.y = element_blank(),  # Remove Y-axis labels
        axis.title.y = element_blank(),  # Remove Y-axis title
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.text.x = element_text(size = 30)  # Set X-axis title font size
  ) 

ggsave("OP/Reanno_GLASSidh_v2.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 600)
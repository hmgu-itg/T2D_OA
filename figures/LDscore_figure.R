library(data.table)
library(ggplot2)
library(RColorBrewer)

ldsc.dt <- data.table(OA.phen=c("AllOA", "KneeOA", "HipOA", "KneeHipOA", "TKR", "THR", "TJR"),
                      rg=c(0.2301,0.2414,0.0783,0.2108,0.256,0.0301,0.1631),
                      se=c(0.026,0.0277,0.0298,0.0278,0.0334,0.0273,0.0291),
                      # p=c(),
                      emp.p=c(0.00654,0.00493,0.1423,0.0001188,0.00348,0.3238,0.0227),
                      cat=c("both", "knee", "hip", "both", "knee", "hip", "both"))

my.colors <- RColorBrewer::brewer.pal(8, "Dark2")[c(1,2,6)]
my.colors <- RColorBrewer::brewer.pal(8, "Paired")[c(3,2,1)]
my.colors <- RColorBrewer::brewer.pal(9, "RdPu")[c(4,6,8)]

ggplot(ldsc.dt, aes(x=rg, y=reorder(OA.phen,rg), xmin=rg-se, xmax=rg+se, fill=cat)) +
# ggplot(ldsc.dt, aes(y=rg, x=reorder(OA.phen,rg), ymin=rg-se, ymax=rg+se, color=cat)) +
  geom_point(size = 3, color="black", pch=21) +
  geom_errorbar(width=0.2,position=position_dodge(0.05), aes(color=cat)) +
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_brewer(palette="Paired") +
  scale_color_brewer(palette="Paired") +
  # scale_fill_manual(values=my.colors) +
  # scale_color_manual(values=my.colors) +
  ylab("") +
  xlab("genetic correlation between T2D and each OA phenotype")
ggsave("/project_data/ldsc_figure.png", width=6, height=3, dpi=400)

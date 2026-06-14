setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot")
library(ggplot2)
library(forcats)
library(readr)
library(SuperExactTest)
library(gridExtra)

clade2_relaxedGO <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO/clade2_relaxed_fdr0.05_p1_GO_BP.csv")
clade5_relaxedGO <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO/clade5_relaxed_fdr0.05_p1_GO_BP.csv")
clade7_relaxedGO <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO/clade7_relaxed_fdr0.05_p1_GO_BP.csv")
clade2_intensifiedGO <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO/clade2_intensified_fdr0.05_p1_GO_BP.csv")
clade5_intensifiedGO <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO/clade5_intensified_fdr0.05_p1_GO_BP.csv")
clade7_intensifiedGO <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO/clade7_intensified_fdr0.05_p1_GO_BP.csv")

clade2_relaxedGO <- subset(clade2_relaxedGO, clade2_relaxedGO$pvalue<0.05)
clade5_relaxedGO <- subset(clade5_relaxedGO, clade5_relaxedGO$pvalue<0.05)
clade7_relaxedGO <- subset(clade7_relaxedGO, clade7_relaxedGO$pvalue<0.05)
clade2_intensifiedGO <- subset(clade2_intensifiedGO, clade2_intensifiedGO$pvalue<0.05)
clade5_intensifiedGO <- subset(clade5_intensifiedGO, clade5_intensifiedGO$pvalue<0.05)
clade7_intensifiedGO <- subset(clade7_intensifiedGO, clade7_intensifiedGO$pvalue<0.05)

clade2_relaxedGODescription <- clade2_relaxedGO$Description
clade5_relaxedGODescription <- clade5_relaxedGO$Description
clade7_relaxedGODescription <- clade7_relaxedGO$Description
clade2_intensifiedGODescription <- clade2_intensifiedGO$Description
clade5_intensifiedGODescription <- clade5_intensifiedGO$Description
clade7_intensifiedGODescription <- clade7_intensifiedGO$Description


set.seed(1234)
#generate random strings
n=29132
input = list(clade2 = arrFamClade2,
             clade5 = arrFamClade5,
             clade7 = arrFamClade7)
Result=supertest(input,n=n)
res_adj <- summary(Result, method="BH")
order2=c("001","010","100","101","011","110","111")
print(order2)
plot(Result, Layout="landscape", sort.by=order2, keep=FALSE, show.elements=FALSE, elements.cex=0.7,
     # 各元素字号独立控制
     set.labels.cex = 2,          # 左侧集合标签（默认1）
     count.cex = 2,               # 顶部数量标签（默认1）
     degree.label.cex = 2,        # 交集组合标签（默认1）
     legend.cex = 2,              # 图例文字（默认1）
     color.scale.name.cex = 2,    # 颜色条标题（默认1）
     elements.list=subset(summary(Result)$Table,Observed.Overlap <= 20),
     show.expected.overlap=TRUE,expected.overlap.style="hatchedBox",color.on = NULL, 
     color.expected.overlap='red',
     elements.list = subset(summary(Result)$Table, Observed.Overlap <= 20)) # 调整图例位置


set.seed(1234)
#generate random strings
n=29132
input = list(clade2 = clade2_contracted,
             clade7 = clade7_contracted)
Result=supertest(input,n=n)
res_adj <- summary(Result, method="BH")
print(order2)
plot(Result, Layout="landscape", keep=FALSE, show.elements=FALSE, elements.cex=0.7,
     # 各元素字号独立控制
     set.labels.cex = 2,          # 左侧集合标签（默认1）
     count.cex = 2,               # 顶部数量标签（默认1）
     degree.label.cex = 2,        # 交集组合标签（默认1）
     legend.cex = 2,              # 图例文字（默认1）
     color.scale.name.cex = 2,    # 颜色条标题（默认1）
     elements.list=subset(summary(Result)$Table,Observed.Overlap <= 20),
     show.expected.overlap=TRUE,expected.overlap.style="hatchedBox",color.on = NULL, 
     color.expected.overlap='red',
     elements.list = subset(summary(Result)$Table, Observed.Overlap <= 20)) # 调整图例位置



fig3GO <- read_csv("fig3GO.csv")
fig3GO$Description <- as.factor(fig3GO$Description)
fig3GO$Description <- fct_inorder(fig3GO$Description)
#首先绘制气泡图，这里把分组信息加上去
plot <- ggplot(fig3GO, aes(Clade, Description)) +
  geom_point(aes(color=pvalue, size=Count))+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  #  scale_color_gradient(low='#8B0000',high='#191970')+
  scale_color_gradient(low='#963634',high='#3E68A6')+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=1))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))
ggsave("fig3GO.tiff", plot, 
width=6, height=8, units="in", dpi=300, 
compression = "lzw")

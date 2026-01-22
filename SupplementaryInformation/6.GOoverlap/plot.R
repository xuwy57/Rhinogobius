setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot")
library(ggplot2)
##"#b72230", "#317cb7"
#####fig2a
#red:intensified;blue:relaxed
#a
data <- data.frame(Clade = c("Clade2", "Clade2", "Clade5", "Clade5", "Clade7", "Clade7"),
                      selective_pressure = rep(c("relaxed", "intensified"), 3),
                      count = c(2299, 236, 1323, 121, 4263, 320))
data$selective_pressure <- factor(data$selective_pressure, levels = c("intensified", "relaxed"))
data <- data[order(data$Clade, data$selective_pressure), ]
ggplot(data, aes(x = Clade, y = count, fill = selective_pressure)) +
  geom_bar(position="stack", stat = "identity") +
  labs(y = "Gene Count", fill = "Selective Pressure") +
  scale_fill_manual(name=NULL,labels = c("Intensified", "Relaxed"),values = c("#dc6d57", "#6dadd1"))+
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), # 在堆积层中央
            color = "black", size = 3.5) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(0,1,1,3),"cm"),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(hjust = 0.5,vjust=0.5, size = 14,color = 'black'),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 12),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA)) 

#bc 
library(ggVennDiagram)
x <- list(Clade2_intensified=arrTargetGenes_clade2,          
          Clade5_intensified=arrTargetGenes_clade5,          
          Clade7_intensified=arrTargetGenes_clade7)  
ggVennDiagram(x, label="count") + scale_fill_gradient(low="#fbe3d5",high = "#6d011f")

x <- list(Clade2_relaxed=arrTargetGenes_clade2,          
          Clade5_relaxed=arrTargetGenes_clade5,      
          Clade7_relaxed=arrTargetGenes_clade7)  
ggVennDiagram(x,, label="count", label_size=8) + scale_fill_gradient(low="#e9f1f4",high = "#104680")

#####fig3 
library(forcats)
library(ggplot2)
library(readr)
fig1GO <- read_csv("fig1GO.csv")
fig1GO <- read_csv("fig1GO_MF.csv")
fig1GO$Description <- as.factor(fig1GO$Description)
fig1GO$Description <- fct_inorder(fig1GO$Description)
#首先绘制气泡图，这里把分组信息加上去
ggplot(fig1GO, aes(Clade, Description)) +
  geom_point(aes(color=pvalue, size=Count))+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
#  scale_color_gradient(low='#8B0000',high='#191970')+
  scale_color_gradient(low='#963634',high='#3E68A6')+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=1))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))


library(ggVennDiagram)
x <- list(Clade2_intensified=arrTargetGenes_clade2_DEG,   
          Clade5_intensified=arrTargetGenes_clade5_DEG,          
          Clade7_intensified=arrTargetGenes_clade7_DEG)  
ggVennDiagram(x) + scale_fill_gradient(low="#F1E4E8",high = "#A35C72")

x <- list(Clade2_relaxed=arrTargetGenes_clade2_DEG,          
          Clade5_relaxed=arrTargetGenes_clade5_DEG,      
          Clade7_relaxed=arrTargetGenes_clade7_DEG)  
ggVennDiagram(x) + scale_fill_gradient(low="#E1E8F1",high = "#3C779F")






#fig3
Sets = c("Clade2 ∩ Clade5 all & relaxed", "Clade2 ∩ Clade7 all & relaxed", "Clade5 ∩ Clade7 all & relaxed", 
         "Clade2 ∩ Clade5 all & intensified", "Clade2 ∩ Clade7 all & intensified", "Clade5 ∩ Clade7 all & intensified",
         "Clade2 ∩ Clade5 overlapped with DEG & relaxed", "Clade2 ∩ Clade7 overlapped with DEG & relaxed", "Clade5 ∩ Clade7 overlapped with DEG & relaxed", 
         "Clade2 ∩ Clade5 overlapped with DEG & intensified", "Clade2 ∩ Clade7 overlapped with DEG & intensified", "Clade5 ∩ Clade7 overlapped with DEG & intensified")
pvalue <- c(3.802264e-214, 1.788926e-266, 4.790301e-197, 5.713664e-17, 5.571169e-14, 6.202359e-14)
p <- round(-log10(pvalue),2)
type <- c("relaxed", "relaxed", "relaxed", "intensified", "intensified", "intensified")
colors <- c("#A35C72","#3C779F")
            #, "#BDB231", "#EEE111", "#706B2C")
data <- data.frame(Sets = Sets, 
                   type = factor(type),
                   p = p)
data$label <- format(data$p, decimal.mark = ",", nsmall = 2)  # 将数字中的小数点替换为逗号                   
#棒棒糖图主体和数值
p1 <- ggplot(data, aes(x = Sets, y = p, color = type)) +  
  geom_point(size = 4) +  
  geom_segment(aes(x = Sets, xend = Sets, y = 0, yend = p), linewidth = 1.5) +  
  geom_text(aes(y = p, label = p),            
            size = 4, angle = 90,            
            nudge_y = 8,            
            hjust = 0,            
            vjust = 0.5) +  
  scale_y_continuous(limits = c(0, 350),                      
                     breaks = seq(0, 100, 20),                     
                     expand = c(0, 0)) +  
  scale_color_manual(values = colors) +  
  labs(y = expression(paste("-log10(pvalue)"))) +  
#  theme_classic2() +   
  theme(legend.position = "none",    
        panel.grid = element_blank(),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.title.y = element_text(size = 16),       
        axis.text.y = element_text(size = 12),        
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),        
#        axis.line.x = element_blank(),        
#        axis.ticks.x = element_blank(),        
#        axis.line.y = element_line(linewidth = 1.5, color = "#808181"),        
#        axis.ticks.y = element_line(linewidth = 1.5, color = "#808181"),        
#        axis.ticks.length.y = unit(0.2, "cm"),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))
p1
#棒棒糖图基底注释条带
p2 <- ggplot(data, aes(x = Sets, y = 0.1, fill = type)) +  
  geom_tile(width = 0.95) + 
  scale_x_discrete(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colors) +  
  theme_bw() +  
  theme(legend.position = "none",        
        panel.border = element_rect(linewidth = 2),        
        panel.grid = element_blank(),        
        axis.title = element_blank(),        
        axis.ticks = element_blank(),        
        axis.text.y = element_blank(),        
        axis.text.x = element_text(size = 15,                                    
                                   color = "#808181",                                    
                                   family = "sans",                                   
                                   hjust = 1,                                   
                                   vjust = 1,                                   
                                   angle = 45))
lollipops <- p1 / p2 + plot_layout(heights = c(1, 0.05))
lollipops
ggsave(file = "lollipops.pdf", plot = lollipops, width = 8, height = 4.7)



# 安装并加载 SuperExactTest 包（如果尚未安装）
if (!requireNamespace("SuperExactTest", quietly = TRUE)) {
  install.packages("SuperExactTest")
}
library(SuperExactTest)

# 构造示例集合
setA <- 1:100          # 集合 A 包含数字 1 到 100
setB <- 50:150         # 集合 B 包含数字 50 到 150
setC <- 80:180         # 集合 C 包含数字 80 到 180

# 将三个集合放入列表中
set_list <- list(A = setA, B = setB, C = setC)

# 进行多集合重叠显著性检验
result <- supertest(set_list)

# 查看检验结果
print(result)

# 绘制结果图（可选）
plot(result)

# 1. 加载R包
library(phytools)
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/raxml/")

# 2. 加载示例数据
tree_file <- "FigTree_pop.tre"  
tree <- read.tree(tree_file)
trait_file <- "get_clade_inf.txt"   # 替换为你的性状文件路径
trait_data <- read.csv(trait_file, row.names = 1, stringsAsFactors = FALSE)
x <- setNames(trait_data[,1], rownames(trait_data))

# 3. 查看数据基本信息
tree                    # 查看树的基本信息（节点数、物种数等）
x                       # 查看性状向量，检查状态类别（如"TG", "Tw", "Tr"等）
# 确保树是二叉的
if (!is.binary(tree)) {
  warning("树包含多分叉，正在将其随机拆分为二叉（添加零长度分支）...")
  tree <- multi2di(tree)
}

# 确保树已定根
if (!is.rooted(tree)) {
  warning("树未定根，正在用第一个物种作为外群进行定根...")
  # 可以选择一个合适的外群，这里用第一个物种作为外群，请根据实际替换
  tree <- root(tree, outgroup = "Glossogobius_giuris", resolve.root = TRUE)
}

# 检查物种名是否与树中的 tip.label 一致
# 不一致的物种会被 ace 函数自动排除，并给出警告
if (!all(names(x) %in% tree$tip.label)) {
  warning("部分物种名在树中未找到，请检查拼写！")
  # 可选择保留交集
  common_species <- intersect(names(x), tree$tip.label)
  x <- x[common_species]
  tree <- keep.tip(tree, common_species)
}


# 假设你的原始树是 tree，性状向量是 x
# 假设你知道某些内部节点的状态，并将它们存储在一个名为 'known_states' 的向量中
# 这个向量的名称是节点编号（node IDs），值是该节点的状态。
# 例如：
# known_states <- c("47" = "landlocked", "53" = "amphidromous")

# --- 步骤1：创建一个新树，用于添加分支 ---
tree_with_known <- tree
# 先简单绘制树（不显示编号）
plotTree(tree, type = "phylogram", fsize = 0.7)
# 在内部节点处添加编号（节点编号为灰色小字）
nodelabels(cex = 0.6, frame = "none", col = "red")
known_states <- c("59" = "Amphidromous", "115" = "Amphidromous", "60" = "Amphidromous"
                  , "61" = "Amphidromous", "63" = "Amphidromous", "64" = "Amphidromous"
                  , "62" = "Amphidromous")
# --- 步骤2：循环为每个已知节点添加零长度分支 ---
# 重要：每添加一个新分支，树的节点编号都会改变，所以需要动态匹配节点。
for (i in 1:length(known_states)) {
  original_node_id <- as.numeric(names(known_states)[i])
  # 使用 matchNodes 将原始树的节点编号，映射到当前新树 (tree_with_known) 的节点编号上
  # method="distances" 通过比较节点到根的距离来匹配，比较稳健
  node_match <- matchNodes(tree, tree_with_known, method = "distances")
  # 获取在当前树中的正确节点编号
  current_node_id <- node_match[which(node_match[,1] == original_node_id), 2]
  # 添加一个长度为 0 的分支到这个节点上，新“种”的名字可以自定义，比如用节点编号
  new_tip_name <- paste0("Node_", original_node_id, "_fixed")
  tree_with_known <- bind.tip(tree_with_known, 
                              tip.label = new_tip_name,
                              edge.length = 0,       # 关键：长度为 0
                              where = current_node_id)
}
# --- 步骤3：构建新的性状向量，包含原始物种和这些“固定”的节点 ---
# 将已知节点的状态赋给新添加的 tip
x_known <- known_states
names(x_known) <- paste0("Node_", names(known_states), "_fixed")
# 合并原始物种的性状和这些新“节点种”的性状
x_with_known <- c(x, x_known)
# 步骤4：使用 fitMk 拟合模型（你可以指定任何模型，如 "ARD" 或自定义矩阵）
# 注意：新树包含了零长度分支，不会影响模型对进化速率的估计，但会提供信息。

# 获取状态顺序
states <- sort(unique(x_with_known))
k <- length(states)
# 创建一个空的 design matrix，先全部填 0
model_custom <- matrix(0, k, k, dimnames = list(states, states))
# 按照我们的生物学假设填入数字索引
# 规则：相同数字的速率被约束相等
# A→B 
model_custom["Amphidromous", "Landlocked"] <- 1
# B→A
model_custom["Landlocked", "Amphidromous"] <- 0.8
# 查看我们构建的 design matrix
print(model_custom)
# 拟合自定义约束模型
fit_custom <- fitMk(tree_with_known, x_with_known, model = model_custom, pi = c(1, 0))

anc_states <- ancr(fit_custom)

# 设置颜色（根据实际性状状态数量调整）
states <- sort(unique(x_with_known))
cols <- setNames(RColorBrewer::brewer.pal(length(states), "Set1"), states)

# 调整画布边距（下、左、上、右）
par(mar = c(1, 1, 1, 1))

# 一个简单的策略是直接用 tree_with_known 绘图，但将用于固定的 tip 的标签设置为空或不显示。
# 获取所有 tip 的坐标
last_pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
# 找到固定用的 tip 的索引
fixed_tips <- grep("_fixed", tree_with_known$tip.label)
# 在这些 tip 的位置上不显示标签（设置颜色为透明或直接不绘制）
# 但这需要更精细的绘图控制，一个更简单的做法是，在绘图前就重命名它们的标签为空字符串
tip_labels_for_plot <- tree_with_known$tip.label
tip_labels_for_plot[fixed_tips] <- ""
plot(tree_with_known, tip.label = tip_labels_for_plot, type = "phylogram", cex=0.8, label.offset=0.1,         # 物种名与分支末端的距离
     use.edge.length = TRUE,   # 显示分支长度
     direction = "rightwards")

# 在内部节点添加饼图（使用 ancr 输出的祖先状态概率）
nodelabels(node = 1:tree_with_known$Nnode + Ntip(tree_with_known),
           pie = anc_states$ace,
           piecol = cols,
           cex = 0.2)

# 在尖端添加饼图（当前状态）
tiplabels(pie = to.matrix(x_with_known[tree_with_known$tip.label], states),
          piecol = cols,
          cex = 0.2,
          offset = 0.1)

# 添加图例
add.simmap.legend(colors = cols,
                  prompt = FALSE,
                  x = 0.1,
                  y = 10,
                  fsize = 0.8)









# 1. 加载R包
library(phytools)
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/raxml/")

# 2. 加载示例数据
tree_file <- "FigTree_pop.tre"  
tree <- read.tree(tree_file)
trait_file <- "get_clade_inf.txt"   # 替换为你的性状文件路径
trait_data <- read.csv(trait_file, row.names = 1, stringsAsFactors = FALSE)
x <- setNames(trait_data[,1], rownames(trait_data))

# 3. 查看数据基本信息
tree                    # 查看树的基本信息（节点数、物种数等）
x                       # 查看性状向量，检查状态类别（如"TG", "Tw", "Tr"等）
# 确保树是二叉的
if (!is.binary(tree)) {
  warning("树包含多分叉，正在将其随机拆分为二叉（添加零长度分支）...")
  tree <- multi2di(tree)
}

# 确保树已定根
if (!is.rooted(tree)) {
  warning("树未定根，正在用第一个物种作为外群进行定根...")
  # 可以选择一个合适的外群，这里用第一个物种作为外群，请根据实际替换
  tree <- root(tree, outgroup = "Glossogobius_giuris", resolve.root = TRUE)
}

# 检查物种名是否与树中的 tip.label 一致
# 不一致的物种会被 ace 函数自动排除，并给出警告
if (!all(names(x) %in% tree$tip.label)) {
  warning("部分物种名在树中未找到，请检查拼写！")
  # 可选择保留交集
  common_species <- intersect(names(x), tree$tip.label)
  x <- x[common_species]
  tree <- keep.tip(tree, common_species)
}
states <- sort(unique(x))
# 步骤4：使用 fitMk 拟合模型（你可以指定任何模型，如 "ARD" 或自定义矩阵）
k <- length(states)
# 创建一个空的 design matrix，先全部填 0
model_custom <- matrix(0, k, k, dimnames = list(states, states))
# 按照我们的生物学假设填入数字索引
# 规则：相同数字的速率被约束相等
# A→B 
model_custom["Amphidromous", "Landlocked"] <- 1
# B→A
model_custom["Landlocked", "Amphidromous"] <- 0.8
# 查看我们构建的 design matrix
print(model_custom)
# 拟合自定义约束模型
fit_custom <- fitMk(tree, x, model = model_custom)

anc <- ancr(fit_custom)

# 设置颜色（根据实际性状状态数量调整）
states <- sort(unique(x))
cols <- setNames(RColorBrewer::brewer.pal(length(states), "Set1"), states)

# 绘制树（直角分支）
plot(tree, type = "phylogram", cex = 0.8, label.offset = 0.1, use.edge.length = TRUE)

# 在内部节点添加饼图（祖先状态概率）
nodelabels(node = 1:tree$Nnode + Ntip(tree),
           pie = anc$ace,
           piecol = cols,
           cex = 0.3)

# 在尖端添加饼图（当前物种的状态）
tiplabels(pie = to.matrix(x[tree$tip.label], states),
          piecol = cols,
          cex = 0.2,
          offset = 0.05)

# 添加图例
add.simmap.legend(colors = cols,
                  prompt = FALSE,
                  x = 0.1,
                  y = 10,
                  fsize = 0.8)





# ====================== 1. 加载包 ======================
library(phytools)
library(RColorBrewer)

setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/raxml/")

# ====================== 2. 加载数据 ======================
tree <- read.tree("FigTree_pop.tre")
trait_data <- read.csv("get_clade_inf.txt", row.names = 1, stringsAsFactors = FALSE)
x <- setNames(trait_data[,1], rownames(trait_data))

# ====================== 3. 数据预处理 ======================
if (!is.binary(tree)) tree <- multi2di(tree)
if (!is.rooted(tree)) {
  tree <- root(tree, outgroup = "Glossogobius_giuris", resolve.root = TRUE)
}

common_species <- intersect(names(x), tree$tip.label)
x <- x[common_species]
tree <- keep.tip(tree, common_species)

states <- sort(unique(x))
cat("States detected:", paste(states, collapse = ", "), "\n\n")

# ====================== 定义三种模型 ======================
k <- length(states)

model_A <- matrix(0, k, k, dimnames = list(states, states))
model_A["Amphidromous", "Landlocked"] <- 1
model_A["Landlocked", "Amphidromous"] <- 0.8

model_B <- matrix(1, k, k, dimnames = list(states, states))

model_C <- matrix(0, k, k, dimnames = list(states, states))
model_C["Amphidromous", "Landlocked"] <- 0.8
model_C["Landlocked", "Amphidromous"] <- 1

models <- list(
  "Asymmetric (A→L=1, L→A=0.8)"     = model_A,
  "Symmetric (both=1)"              = model_B,
  "Reverse Asymmetric (A→L=0.8, L→A=1)" = model_C
)

# ====================== 拟合模型并比较（更稳健的写法） ======================
results <- data.frame(
  Model  = character(0),
  logLik = numeric(0),
  npar   = numeric(0),
  AIC    = numeric(0),
  stringsAsFactors = FALSE
)

anc_list <- list()

cat("=== Model Fitting and Comparison ===\n\n")

for (name in names(models)) {
  cat("Fitting model:", name, "...\n")
  
  fit <- fitMk(tree, x, model = models[[name]])
  anc_list[[name]] <- ancr(fit)
  
  cat("  log-likelihood :", round(fit$logLik, 3), "\n")
  cat("  Number of parameters :", fit$numpar, "\n")
  cat("  AIC                :", round(AIC(fit), 2), "\n\n")
}
# ====================== 保存最佳模型的结果 ======================
best_model_name <- "Symmetric (both=1)"
anc_best <- anc_list[[best_model_name]]

cat("\nBest model selected:", best_model_name, "\n")

# ====================== 6. 报告最佳模型的关键祖先概率 ======================
best_model_name <- "Asymmetric (A→L=1, L→A=0.8)"
anc_best <- anc_list[[best_model_name]]

cat("\n=== Ancestral Probabilities under Best Model (Asymmetric) ===\n")

# MRCA (根节点)
root_prob <- anc_best$ace[1, "Amphidromous"]
cat("MRCA (Root) P(Amphidromous) =", round(root_prob*100, 2), "%\n")

# 请你先运行 plot(tree) 找到以下三个 Clade 的祖先节点编号，然后替换下面的数字
clade2_node <- 59   # ←←← 替换为 Clade 2 的祖先节点编号
clade5_node <- XXX   # ←←← 替换为 Clade 5 的祖先节点编号
clade7_node <- XXX   # ←←← 替换为 Clade 7 的祖先节点编号

cat("Clade 2  ancestor P(Amphidromous) =", round(anc_best$ace[clade2_node - Ntip(tree), "Amphidromous"]*100, 2), "%\n")
cat("Clade 5  ancestor P(Amphidromous) =", round(anc_best$ace[clade5_node - Ntip(tree), "Amphidromous"]*100, 2), "%\n")
cat("Clade 7  ancestor P(Amphidromous) =", round(anc_best$ace[clade7_node - Ntip(tree), "Amphidromous"]*100, 2), "%\n")

# ====================== 7. 绘图（最佳模型） ======================
cols <- setNames(brewer.pal(length(states), "Set1"), states)

png("Ancestral_Reconstruction_BestModel.png", width = 5000, height = 3800, res = 300)

plot(tree, type = "phylogram", cex = 0.75, label.offset = 0.20, 
     use.edge.length = TRUE, no.margin = TRUE)

nodelabels(node = 1:tree$Nnode + Ntip(tree),
           pie = anc_best$ace,
           piecol = cols,
           cex = 0.35)

tiplabels(pie = to.matrix(x[tree$tip.label], states),
          piecol = cols,
          cex = 0.25, offset = 0.08)

add.simmap.legend(colors = cols, prompt = FALSE, x = 0.05, y = 15, fsize = 0.9)

dev.off()

cat("\nBest model plot saved as: Ancestral_Reconstruction_BestModel.png\n")


# 绘制树并显示节点编号
plot(tree, type = "phylogram", cex = 0.8, label.offset = 0.1)
nodelabels(text = 1:tree$Nnode + Ntip(tree), cex = 0.5, frame = "none", col = "blue")













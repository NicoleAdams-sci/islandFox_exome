# load in libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

library(tidyverse)
library(ggtree)
library(phangorn)


# load in metadata
#mydata <- read.table("/Users/adamsn23/Downloads/ExomeSamples_forMap_4-26-2018.csv", header=TRUE, sep=",")
#mydata <- read.table("/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/ExomeSamples_forMap_4-26-2018.csv", header=TRUE, sep=",")
mydata <- read.table("~/Documents/NicoleAdams/ulit/ExomeSamples_forMap_4-26-2018.csv", header=TRUE, sep=",")
mydata.1 <- subset(mydata, subset = mydata$Time=="Historical")
mydata.2 <- subset(mydata, subset = mydata$Time=="Modern")

mydata$ISLAND <- factor(mydata$ISLAND, levels=c("SMI", "SRI", "SCZ", "CAT", "SCL", "SNI"))
island.cols <- c("midnightblue", "royalblue2", "cadetblue2", "darkred", "firebrick1", "darksalmon")

# load in tree data
tree <- read.tree("/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/outtree")

# manipulate metadata to match tree
mydata$treeIDaa <- gsub("_", "-", mydata$Sample_ID) 
mydata$treeIDa <- gsub("LANHM", "LA", mydata$treeIDaa)
mydata$treeIDb <- gsub("SBMNH", "SB", mydata$treeIDa)
mydata$treeIDc <- gsub("CAT-51-20852", "CAT-51", mydata$treeIDb)
mydata$treeIDd <- gsub("CAT-55-45F62", "CAT-55", mydata$treeIDc)
mydata$treeIDe <- gsub("CAT-56-F7507", "CAT-56", mydata$treeIDd)
mydata$treeIDf <- gsub("CAT-58-B0F74", "CAT-58", mydata$treeIDe)
mydata$treeIDg <- gsub("CAT-59-1477A", "CAT-59", mydata$treeIDf)
mydata$treeIDh <- gsub("CAT-60-35E15", "CAT-60", mydata$treeIDg)
mydata$treeID <- gsub("CAT-63-97338", "CAT-63", mydata$treeIDh)

tree$tip.label <- gsub("SCL-40749", "SNI-SCL-40749", tree$tip.label)


# remove samples from meta not in tree
samples2rm <- c("LA_006843", "LA_085730", "SB_1074", "SB_2063")
mydata4tree <- mydata %>% filter(!treeID %in% samples2rm) %>% dplyr::select(treeID, ISLAND, Time)

# get internal node numbers
# ggtree(tree) + geom_text2(aes(label=node))
# get close up on tight clade
# viewClade(tt, MRCA(tt, tip=c("LA-008354", "SB-2937"))) + geom_text2(aes(label=node))

tt <- ggtree(tree) #+ geom_treescale()
#MRCA(tree.p, "SNI-40530", "LA-008354")

tree.p <- tt %<+% mydata4tree + 
  geom_tiplab(aes(color = factor(Time)), # color for label font
              geom = "text", size=3, show.legend=F) +
  scale_color_grey(start=0.1, end=0.6) +
  geom_cladelabel(node=156, label = "CAT", color = "darkred", barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=142, label = "SNI", color='darksalmon', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=140, label = "CAT", color = "darkred", barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=126, label = "SCL", color='firebrick1', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=110, label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_strip(taxa1 = 'SB-1246', taxa2 = 'SB-1245', label = "SRI", color='royalblue2', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE) +
 # geom_cladelabel(node=36, label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_strip(taxa1 = 'LA-006391', taxa2 = 'SB-1246', label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE) +
 # geom_cladelabel(node=90, label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) 
  geom_strip(taxa1 = 'SB-1931', taxa2 = 'LA-007729', label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE)

  # geom_strip('68', '82', label = "CAT", color='darkred', barsize = 2, offset = 0.01, fontsize = 6) +
  # geom_strip('64', '55', label = "SNI", color='darksalmon', barsize = 2, offset = 0.01, fontsize = 6) +
  # geom_strip('42', '49', label = "SCL", color='firebrick1', barsize = 2, offset = 0.01, fontsize = 6) +
  # geom_strip('33', '20', label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6) +
  # geom_strip('108', '36', label = "SRI", color='royalblue2', barsize = 2, offset = 0.01, fontsize = 6) +
  # geom_strip('5', '88', label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.01, fontsize = 6) 

#ggsave(tree.p, filename = "/Users/Nicole/Documents/Fox_lab/MHC_Sequencing/ms_drafts/tree_6-8-2021", width = 11, height =8 )
#ggsave(tree.p, filename = "/Users/adamsn23/Downloads/tree_2-1-2022.png", width = 11, height =8)




## bootstrap 
# https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html#standard-scripts-for-nucleotide-analysis
library(phangorn)
library(tidyverse)
library(ggtree)

#fasta.f <- read.phyDat("/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/snphylo.output.fasta", format = "fasta")
fasta.f <- read.phyDat("/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/snphylo.output_bwaAln.fasta", format = "fasta")


dm  <- dist.ml(fasta.f)
treeUPGMA  <- upgma(dm)

plot(treeUPGMA, main="UPGMA")

fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(fasta.f,  fun)

plotBS(treeUPGMA, bs_upgma, main="UPGMA")


parsimony(treeUPGMA, fasta.f)

treeRatchet  <- pratchet(fasta.f, trace = 0, minit=100)
parsimony(treeRatchet, fasta.f)

treeRatchet  <- acctran(treeRatchet, fasta.f)
treeRatchet  <- di2multi(treeRatchet)

if(inherits(treeRatchet, "multiPhylo")){
  treeRatchet <- unique(treeRatchet)
}

plotBS(midpoint(treeRatchet), type="phylogram")
add.scale.bar()

#mt <- modelTest(fasta.f)

mt <- modelTest(fasta.f, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
                control = pml.control(trace = 0))



fit_mt <- pml_bb(mt, control = pml.control(trace = 0))
fit_mt


fit <- pml(treeUPGMA, fasta.f)
fitHKY <- optim.pml(fit, model = "HKY")


logLik(fitHKY)
bs <- bootstrap.pml(fitHKY, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
#saveRDS(bs, file = "/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/bootstrapTree_1-10-2023.RDS")

plotBS(midpoint(fitHKY$tree), bs, p=50, type="p") #shows bootsraps >50 support

#write.tree(plotBS(midpoint(fitHKY$tree), bs, type="p"), file = "/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/bootstrapTree.txt")
write.tree(plotBS(midpoint(fitHKY$tree), bs, type="p"), file = "/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/bootstrapTree_bwaAln.txt")

#bs.tr <- read.tree("/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/bootstrapTree.txt")
#bs.tr <- read.tree("/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/bootstrapTree_bwaAln.txt")
bs.tr <- read.tree("~/Documents/NicoleAdams/ulit/bootstrapTree_bwaAln.txt") #umich imac
bs.tt <- ggtree(bs.tr)

bs.tt$data$label <- gsub("SCL-40749", "SNI-SCL-40749", bs.tt$data$label)



# color bootstrap nodes code from doi: 10.1371/journal.pone.0243927

tree.p2 <- bs.tt %<+% mydata4tree + 
  geom_tiplab(aes(color = factor(Time)), # color for label font
              geom = "text", size=6, show.legend=F, hjust = -.03) +
  scale_color_grey(start=0.1, end=0.6) +
  #geom_nodelab(geom="text", aes(subset=!is.na(as.numeric(label)) & as.numeric(label) > 50)) +
  geom_nodepoint( shape=21, size=4, aes(subset=!is.na(as.numeric(label)),
                                        fill=cut(as.numeric(label),c(0,50,70,85,100)))) +
  scale_fill_manual(values = c("black","mediumpurple1","lavender","white"), guide = 'legend',
                    name = 'Bootstrap Percentage',
                    breaks = c('(85,100]', '(70,85]','(50,70]', '(0,50]'),
                    labels=c("> 85", "70-85", "50-70", "< 50")) + 
  geom_cladelabel(node=92, label = "CAT", color = "darkred", barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=113, label = "SNI", color='darksalmon', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=127, label = "SCL", color='firebrick1', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=148, label = "SMI", color='midnightblue', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
 # geom_cladelabel(node=163, label = "SRI", color='royalblue2', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
 # geom_cladelabel(node=141, label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  geom_strip(taxa1 = 'SRI-38956', taxa2 = 'SB-1369', label = "SRI", color='royalblue2', barsize = 2, offset = 0.02, fontsize = 6, hjust =-0.3, align = TRUE) +
  geom_strip(taxa1 = 'LA-006391', taxa2 = 'SB-1369', label = "SMI", color='midnightblue', barsize = 2, offset = 0.02, fontsize = 6, hjust =-0.3, align = TRUE) + 
  geom_strip(taxa1 = 'SB-1931', taxa2 = 'LA-007729', label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.02, fontsize = 6, hjust =-0.3, align = TRUE)


tree.bs.p <- tree.p2 + theme(legend.position = c(0,0.95), legend.justification = c(0,0.85),
                             legend.text = element_text(size=14),
                             legend.title = element_text(size=16))

#ggsave(tree.bs.p, filename = "/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/tree_bootstrap_1-11-2023.png", width = 11, height =16)









########## Tree without labels

tree.p3 <- bs.tt %<+% mydata4tree + 
  #geom_tiplab(aes(color = factor(Time)), # color for label font
              #geom = "text", size=6, show.legend=F, hjust = -.03) +
  #scale_color_grey(start=0.1, end=0.6) +
  geom_tippoint(aes(color=ISLAND, shape=Time), size=4, show.legend = F) +
  scale_color_manual(values=island.cols, guide="none") +
  #geom_nodelab(geom="text", aes(subset=!is.na(as.numeric(label)) & as.numeric(label) > 50)) +
  geom_nodepoint( shape=23, size=3, aes(subset=!is.na(as.numeric(label)),
                                        fill=cut(as.numeric(label),c(0,50,70,85,100)))) +
  scale_fill_manual(values = c("black","mediumpurple1","lavender","white"), 
                    name = 'Bootstrap Percentage',
                    breaks = c('(85,100]', '(70,85]','(50,70]', '(0,50]'),
                    labels=c("> 85", "70-85", "50-70", "< 50"),
                    guide = guide_legend(override.aes = list(shape=23))) + 
    geom_cladelabel(node=92, label = "CAT", color = "darkred", barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=113, label = "SNI", color='darksalmon', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=127, label = "SCL", color='firebrick1', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=148, label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  # geom_cladelabel(node=163, label = "SRI", color='royalblue2', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  # geom_cladelabel(node=141, label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_strip(taxa1 = 'SRI-38956', taxa2 = 'SB-1369', label = "SRI", color='royalblue2', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE) +
  geom_strip(taxa1 = 'LA-006391', taxa2 = 'SB-1369', label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE) + 
  geom_strip(taxa1 = 'SB-1931', taxa2 = 'LA-007729', label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE)


tree.bs.p2 <- tree.p3 + theme(legend.position = c(0,0.95), legend.justification = c(0,0.85),
                             legend.text = element_text(size=14),
                             legend.title = element_text(size=16))

#ggsave(tree.bs.p2, filename = "/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/tree_bootstrap_nolabs_1-12-2023.png", width = 11, height =16)








##################### FOR BWA ALN ##################
# define bootstrap colors
colfunc <- colorRampPalette(c("gray30", "white"))
colfunc <- colorRampPalette(c("white", "gray30"))
boot.cols <- colfunc(4)

tree.p2 <- bs.tt %<+% mydata4tree + 
  geom_tiplab(aes(color = factor(Time)), # color for label font
              geom = "text", size=6, show.legend=F, hjust = -.03) +
  scale_color_grey(start=0.1, end=0.6) +
  #geom_nodelab(geom="text", aes(subset=!is.na(as.numeric(label)) & as.numeric(label) > 50)) +
  geom_nodepoint( shape=21, size=4, aes(subset=!is.na(as.numeric(label)),
                                        fill=cut(as.numeric(label),c(0,50,70,85,100)))) +
  #scale_fill_manual(values = c("black","mediumpurple1","lavender","white"), guide = 'legend',
  #scale_fill_grey(guide = 'legend',
  scale_fill_manual(values = c(boot.cols), guide = 'legend',
                    name = 'Bootstrap Percentage',
                    breaks = c('(85,100]', '(70,85]','(50,70]', '(0,50]'),
                    labels=c("> 85", "70-85", "50-70", "< 50")) + 
  geom_cladelabel(node=87, label = "CAT", color = "darkred", barsize = 2, offset = 0.04, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=106, label = "SNI", color='darksalmon', barsize = 2, offset = 0.04, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=119, label = "SCL", color='firebrick1', barsize = 2, offset = 0.04, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=138, label = "SMI", color='midnightblue', barsize = 2, offset = 0.04, fontsize = 6, align = TRUE) +
  # geom_cladelabel(node=163, label = "SRI", color='royalblue2', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  # geom_cladelabel(node=141, label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  geom_strip(taxa1 = 'SRI-38956', taxa2 = 'SB-1369', label = "SRI", color='royalblue2', barsize = 2, offset = 0.04, fontsize = 6, hjust =-0.3, align = TRUE) +
  geom_strip(taxa1 = 'LA-006391', taxa2 = 'SB-1369', label = "SMI", color='midnightblue', barsize = 2, offset = 0.04, fontsize = 6, hjust =-0.3, align = TRUE) + 
  geom_strip(taxa1 = 'SCZ-12939', taxa2 = 'LA-007729', label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.04, fontsize = 6, hjust =-0.3, align = TRUE)


tree.bs.p <- tree.p2 + theme(legend.position = c(0,0.95), legend.justification = c(0,0.85),
                             legend.text = element_text(size=14),
                             legend.title = element_text(size=16))

#ggsave(tree.bs.p, filename = "/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/bwaAln_tree_bootstrap_1-29-2023.png", width = 11, height =16)


tree.p3 <-  bs.tt %<+% mydata4tree + 
  #geom_tiplab(aes(color = factor(Time)), # color for label font
     #geom = "text", size=6, show.legend=F, hjust = -.03) +
  #scale_color_grey(start=0.1, end=0.6) +
  geom_tippoint(aes(color=ISLAND, shape=Time), size=4, show.legend = F) +
  scale_color_manual(values=island.cols, guide="none") +
  #geom_nodelab(geom="text", aes(subset=!is.na(as.numeric(label)) & as.numeric(label) > 50)) +
  geom_nodepoint( shape=22, size=4, aes(subset=!is.na(as.numeric(label)),
                                        fill=cut(as.numeric(label),c(0,50,70,85,100)))) +
  #scale_fill_manual(values = c("black","mediumpurple1","lavender","white"), guide = 'legend',
  #scale_fill_grey(guide = 'legend',
  scale_fill_manual(values = c(boot.cols), guide = 'legend',
                    name = 'Bootstrap Percentage',
                    breaks = c('(85,100]', '(70,85]','(50,70]', '(0,50]'),
                    labels=c("> 85", "70-85", "50-70", "< 50")) + 
  geom_cladelabel(node=87, label = "CAT", color = "darkred", barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=106, label = "SNI", color='darksalmon', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=119, label = "SCL", color='firebrick1', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  geom_cladelabel(node=138, label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6, align = TRUE) +
  # geom_cladelabel(node=163, label = "SRI", color='royalblue2', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  # geom_cladelabel(node=141, label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.02, fontsize = 6, align = TRUE) +
  geom_strip(taxa1 = 'SRI-38956', taxa2 = 'SB-1369', label = "SRI", color='royalblue2', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE) +
  geom_strip(taxa1 = 'LA-006391', taxa2 = 'SB-1369', label = "SMI", color='midnightblue', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE) + 
  geom_strip(taxa1 = 'SCZ-12939', taxa2 = 'LA-007729', label = "SCZ", color='cadetblue2', barsize = 2, offset = 0.01, fontsize = 6, hjust =-0.3, align = TRUE)


tree.bs.p3 <- tree.p3 + theme(legend.position = c(0,0.95), legend.justification = c(0,0.85),
                             legend.text = element_text(size=14),
                             legend.title = element_text(size=16))

#ggsave(tree.bs.p3, filename = "/Users/adamsn23/OneDrive - Michigan State University/Documents/Ulit/exome/bwaAln_tree_bootstrap_nolabs_1-29-2023.png", width = 11, height =16)
ggsave(tree.bs.p3, filename = "~/Documents/NicoleAdams/ulit/bwaAln_tree_bootstrap_nolabs_gray_2-14-2023.png", width = 11, height =16)

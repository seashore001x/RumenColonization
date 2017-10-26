library(dada2)
library(phyloseq)
library(plyr)
library(reshape2)
library(ggplot2)
library(DESeq2)
tax <- readRDS('../dada2/tax_final.rds')
seq <- readRDS('../dada2/seqtab_final.rds')
tre <- import_qiime(treefilename = 'tree.nwk')
rep <- import_qiime(refseqfilename = 'uniqueSeqs.fasta')

# Handoff to phyloseq ----
# Make a data.frame holding the sample data
samples.out <- rownames(seq)
subject <- substr(samples.out, 1, 2)
time <- factor(unlist(regmatches(samples.out, gregexpr('[0-9.]+(?=.raw)',samples.out, perl = T))), levels = c('0.5', '2', '4', '8', '16', '24', '48'))
treatment <- substr(samples.out, 3, 4)
samdf <- data.frame(Subject=subject, Time=time, Treatment=treatment)

rownames(samdf) <- samples.out

# Construct phyloseq object (straightforward from dada2 outputs) ----
seq <- t(as.matrix(seq))

rownames(seq) <- paste0("Seq", seq(nrow(seq)))
rownames(tax) <- paste0("Seq", seq(nrow(seq)))

# correct the sample order
colnames(seq)[grep('[12]A[AC]48', colnames(seq))] <- c("S2AA48.raw", "S2AC48.raw", "S1AA48.raw", "S1AC48.raw")
cache <- colnames(seq)[grep('2CC', colnames(seq))] 
cache2 <- colnames(seq)[grep('2CA', colnames(seq))]
colnames(seq)[c(grep('2CC', colnames(seq)), grep('2CA', colnames(seq)))] <- c(cache2, cache)

ps <- phyloseq(otu_table(seq, taxa_are_rows=T), 
               sample_data(samdf), 
               tax_table(tax),
               phy_tree(tre),
               refseq(rep))


# rarefying----
set.seed(20171012)
ps_rare <- rarefy_even_depth(ps)

# Taxa filtering
# Filter archeae & Chloroplast data
ps1 <- subset_taxa(ps_rare, Kingdom != 'Archaea') # exclude archeae data
ps1 <- subset_taxa(ps1, Phylum != 'Cyanobacteria/Chloroplast') # exclude Chloroplast data

# the lactobacilus is too high in sheep 2 at 8 h and 16 h under CA treatment; 
# the strepotococcus is too high in sheep 
ps1 <- prune_samples(!((sample_data(ps1)$Time == 8|sample_data(ps1)$Time == 16)&sample_data(ps1)$Subject == 'S2'&sample_data(ps1)$Treatment == 'CA'), ps1)



# Prevalence filtering ----
# Define prevalence of each taxa
# (in how many samples did each taxa appear at least once)
prev0 <- apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(ps1),
                    tax_table(ps1))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 3 # THe number can be defined by yoursef

# filter taxas appered lower than 3 times
ps2 = prune_taxa((prev0 >= prevalenceThreshold), ps1)

# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps2, Phylum %in% names(keepPhyla))


# Taxa prevalence v. total counts.
ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)

# Abundance value transformation ----
plot_abundance = function(physeq, ylabn = "",
                          Facet = "Phylum",
                          Color = "Phylum"){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "Time", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + ylab(ylabn) +
    scale_y_log10()
}

ps2ra = transform_sample_counts(ps2, function(x){x / sum(x)})

ps2ra_AA <- subset_samples(ps2ra, Treatment == 'AA')
ps2ra_AC <- subset_samples(ps2ra, treatment == 'AC')
ps2ra_CA <- subset_samples(ps2ra, treatment == 'CA')
ps2ra_CC <- subset_samples(ps2ra, treatment == 'CC')
ps2_AA <- subset_samples(ps2, Treatment == 'AA')
ps2_AC <- subset_samples(ps2, Treatment == 'AC')
ps2_CA <- subset_samples(ps2, Treatment == 'CA')
ps2_CC <- subset_samples(ps2, Treatment == 'CC')

plot_abundance(ps2ra_AA,"Relative Abundances")
plot_abundance(ps2ra_AC,"Relative Abundances")

# bar plot
library(wesanderson)
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

plot_bar(ps2ra, "Time", fill="Phylum", facet_grid=~Treatment) + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  scale_fill_manual(values=phylum_colors)+
  scale_color_manual(values=phylum_colors)+
  theme_bw()

detach(package:wesanderson, unload = T)


# Alpha diversity ----
alpha_diversity <- estimate_richness(ps, measures = c('Observed', 'Chao1', 'Shannon', 'Simpson'))
alpha_diversity <- merge(alpha_diversity, samdf, by = 0, sort = F)
alpha_diversity <-
  transform(alpha_diversity, row.names = Row.names, Row.names = NULL)

alpha_diversity_melt <- melt(alpha_diversity[,-3], id = c('Subject', 'Time', 'Treatment'))

alpha_diversity_mean <- ddply(alpha_diversity, .(Time, Treatment), numcolwise(mean))[,-5]
alpha_diversity_se <- ddply(alpha_diversity, .(time, Treatment), numcolwise(function(x)sd(x)/sqrt(length(x))))[,-5]

ggplot(alpha_diversity_melt, aes(x = as.numeric(Time), y = value, color = Treatment)) + geom_smooth(method = "glm", formula = y ~ splines::bs(x, 3), se = T) + theme_bw() + scale_x_continuous(
  breaks = c(1, 2, 3, 4, 5, 6, 7),
  label = c(0.5, 2, 4, 8, 16, 24, 28)
) + labs(x = "Incubated Time, /h", y = '') + facet_wrap(~variable, scales = "free_y")


alpha_diversity_plot <- plot_richness(ps, x = 'Time' ,measures = c('Observed', 'Chao1', 'Shannon', 'Simpson'), color="Treatment") + geom_point()


# Beta diversity ----
ps_ord <- ordinate(ps_rare, "MDS", "bray")
plot_ordination(ps2, ps_ord, type="samples", color="Subject", shape = 'Treatment', axes = c(1,2))+geom_point(size = 3) + theme_bw()

# PCoA after log transform
pslog <- transform_sample_counts(ps2ra, function(x) log(1 + x))
out.bc.log <- ordinate(pslog, method = "MDS", distance = "bray")
evals <- out.bc.log$values$Eigenvalues
plot_ordination(ps2ra, out.bc.log, color = "Time", shape = 'Treatment') +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Subject") +
  geom_point(size = 2.5)


# convert phyloseq data into genus level
ps_genus <- tax_glom(ps2, taxrank = 'Genus', NArm = T) # taxa that not assigned to genus are not wanted
ps_genus_AA <- tax_glom(ps2_AA, taxrank = 'Genus', NArm = T)
ps_genus_AC <- tax_glom(ps2_AC, taxrank = 'Genus', NArm = T)
ps_genus_CA <- tax_glom(ps2_CA, taxrank = 'Genus', NArm = T)
ps_genus_CC <- tax_glom(ps2_CC, taxrank = 'Genus', NArm = T)

ps_genus_AA_merge <- merge_samples(ps_genus_AA, "Time")
genus_table_AA_merge <- t(otu_table(ps_genus_AA_merge))
rownames(genus_table_AA_merge) <- tax_table(ps_genus_AA_merge)[,'Genus']
write.table(genus_table_AA_merge, 'genus_table_AA_merge.txt', sep = '\t', col.names = NA, quote = F)

ps_genus_AC_merge <- merge_samples(ps_genus_AC, "Time")
genus_table_AC_merge <- t(otu_table(ps_genus_AC_merge))
rownames(genus_table_AC_merge) <- tax_table(ps_genus_AC_merge)[,'Genus']
write.table(genus_table_AC_merge, 'genus_table_AC_merge.txt', sep = '\t', col.names = NA, quote = F)

ps_genus_CA_merge <- merge_samples(ps_genus_CA, "Time")
genus_table_CA_merge <- t(otu_table(ps_genus_CA_merge))
rownames(genus_table_CA_merge) <- tax_table(ps_genus_CA_merge)[,'Genus']
write.table(genus_table_CA_merge, 'genus_table_CA_merge.txt', sep = '\t', col.names = NA, quote = F)

ps_genus_CC_merge <- merge_samples(ps_genus_CC, "Time")
genus_table_CC_merge <- t(otu_table(ps_genus_CC_merge))
rownames(genus_table_CC_merge) <- tax_table(ps_genus_CC_merge)[,'Genus']
write.table(genus_table_CC_merge, 'genus_table_CC_merge.txt', sep = '\t', col.names = NA, quote = F)


# Differential abundance testing with ANCOM-----
library(ancom.R)

# Ancom analysis
# convert to ancom format
anc <- data.frame(merge(t(otu_table(ps_genus)), sample_data(ps_genus), by = 0))
rownames(anc) <- anc[,1]
anc <- anc[,-1]

colnames(anc)[na.omit(match(colnames(anc), rownames(tax_table(ps_genus))))] <- as.character((as.data.frame(tax_table(ps_genus))$Genus))

# compare each treatment within different time points
anc_sub <- list()
anc_res <- list()
anc_tax <- list()
anc_plot <- list()
for(i in 1:length(levels(anc$Time))){
  j = levels(anc$Time)[i]
  anc_sub[[i]] <- anc[anc$Time ==j,] 
  anc_res[[i]] <- ANCOM(anc_sub[[i]][,-c(158,159)], multcorr = 3, repeated=F)
  anc_tax[[i]] <- anc_res[[i]]$detected
  anc_plot[[i]] <- plot_ancom(anc_res[[i]])
}

# compare each time points within each treatment
anc_time <- list()
anc_time_res <- list()
anc_time_tax <- list()
anc_time_plot <- list()
for(i in 1:length(levels(anc$Treatment))){
  j = levels(anc$Treatment)[i]
  anc_time[[i]] <- anc[anc$Treatment == j,]
  anc_time_res[[i]] <- ANCOM(anc_time[[i]][,-c(158,160)], multcorr = 3, repeated=F)
  anc_time_tax[[i]] <- anc_time_res[[i]]$detected
  anc_time_plot[[i]] <- plot_ancom(anc_time_res[[i]])
}

  
detach(ancom.R)

# ANOSIM
anosim <- anc[,-c(158,159,160)]
anosim(anosim[c(1:14,29:42,57:70),], anc$Treatment[c(1:14,29:42,57:70)], permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))
anosim(anosim[c(15:28,43:56,71:84),], anc$Treatment[c(15:28,43:56,71:84)], permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))


# write convert phyloseq object into BIOM file
write.table(otu_table(ps_rare), 'OTU_table_rarefied.txt', sep = '\t', quote = F, col.names=NA)
write.table(sample_data(ps_rare), 'sample_data.txt', sep = '\t', quote =F, col.names = NA)
write.table(otu_table(ps_genus), 'Genus_table.txt', sep = '\t', quote = F, col.names=NA)




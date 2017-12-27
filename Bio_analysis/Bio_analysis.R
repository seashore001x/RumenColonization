library(phyloseq)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tibble)
source('../../Rrumen/phyloseq_to_lefse.R')

# import data
tax_gg <- readRDS('../dada2/tax_gg.rds')
tax_rdp <- readRDS('../dada2/tax_rdp.rds')
tax_silva <- readRDS('../dada2/tax_silva.rds')
seq <- readRDS('../dada2/seqtab_final.rds')
tre <- import_qiime(treefilename = 'tree.nwk')
rep <- import_qiime(refseqfilename = 'uniqueSeqs.fasta')

# choose one of the assign taxa
tax <- tax_rdp

# Handoff to phyloseq ----
# Make a data.frame holding the sample data
samples.out <- rownames(seq)
subject <- substr(samples.out, 1, 2)
time <- factor(unlist(regmatches(samples.out, gregexpr('[0-9.]+(?=.raw)',samples.out, perl = T))), levels = c('0.5', '2', '4', '8', '16', '24', '48'))
treatment <- substr(samples.out, 3, 4)
rumenenv <- substr(treatment, 1, 1)
forage <- substr(treatment, 2, 2)
samdf <- data.frame(Subject=subject, Time=time, Treatment=treatment, RumenEnv=rumenenv, Forage=forage)  # sample data
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
# change 16h in AC treatment between S2 and S3
colnames(seq)[match(c('S2AC16.raw', 'S3AC16.raw'), colnames(seq))] <- c('S3AC16.raw', 'S2AC16.raw')

# construct phyloseq file
ps <- phyloseq(otu_table(seq, taxa_are_rows=T), 
               sample_data(samdf), 
               tax_table(tax),
               phy_tree(tre),
               refseq(rep))


# rarefying----
set.seed(20171012)
ps_rare <- rarefy_even_depth(ps)

# Taxa filtering----
# Filter archeae & Chloroplast data
ps1 <- subset_taxa(ps_rare, Kingdom != 'Archaea') # exclude archeae data

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
prevalenceThreshold = 3 # THe number can be defined by yourself

# filter taxas appered lower than 3 times
ps2 = prune_taxa((prev0 >= prevalenceThreshold), ps1)

# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps2, Phylum %in% names(keepPhyla))


# # Taxa prevalence v. total counts.
# ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
#   geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
#   geom_point(size = 2, alpha = 0.7) +
#   scale_y_log10() + scale_x_log10() +
#   xlab("Total Abundance") +
#   facet_wrap(~Phylum)
# 
# # Abundance value transformation ----
# plot_abundance = function(physeq, ylabn = "",
#                           Facet = "Phylum",
#                           Color = "Phylum"){
#   mphyseq = psmelt(physeq)
#   mphyseq <- subset(mphyseq, Abundance > 0)
#   ggplot(data = mphyseq,
#          mapping = aes_string(x = "Time", y = "Abundance",
#                               color = Color, fill = Color)) +
#     geom_violin(fill = NA) +
#     geom_point(size = 1, alpha = 0.3,
#                position = position_jitter(width = 0.3)) +
#     facet_wrap(facets = Facet) + ylab(ylabn) +
#     scale_y_log10()
# }

ps2ra = transform_sample_counts(ps2, function(x){x / sum(x)})

ps2ra_AA <- subset_samples(ps2ra, Treatment == 'AA')
ps2ra_AC <- subset_samples(ps2ra, treatment == 'AC')
ps2ra_CA <- subset_samples(ps2ra, treatment == 'CA')
ps2ra_CC <- subset_samples(ps2ra, treatment == 'CC')
ps2_AA <- subset_samples(ps2, Treatment == 'AA')
ps2_AC <- subset_samples(ps2, Treatment == 'AC')
ps2_CA <- subset_samples(ps2, Treatment == 'CA')
ps2_CC <- subset_samples(ps2, Treatment == 'CC')

# phylum bar plot
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
ps_ord <- ordinate(ps2, "MDS", "bray")
plot_ordination(ps2, ps_ord, type="samples", color='Subject', shape = 'Treatment', axes = c(1,2))+geom_point(size = 3) + theme_bw()

# Beta diversity for each invidiual
ps_ord <- ordinate(subset_samples(ps2, Subject == 'S2'), 'MDS', 'bray')
plot_ordination(ps2, ps_ord, type="samples", color="Time", shape = 'Treatment', axes = c(1,2)) + theme_bw() + geom_point(size = 3)

# PCoA after log transform
pslog <- transform_sample_counts(ps2ra, function(x) log(1 + x))
out.bc.log <- ordinate(pslog, method = "MDS", distance = "bray")
evals <- out.bc.log$values$Eigenvalues
plot_ordination(ps2ra, out.bc.log, color = "Subject", shape = 'Treatment') +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Subject") +
  geom_point(size = 2.5)


# convert phyloseq data into different levels
ps_phylum <- tax_glom(ps2, taxrank = 'Phylum', NArm = T) # taxa that not assigned to genus are not wanted
ps_class <- tax_glom(ps2, taxrank = 'Class', NArm = T)
ps_order <- tax_glom(ps2, taxrank = 'Order', NArm = T)
ps_family <- tax_glom(ps, taxrank = 'Family', NArm = T)
ps_genus <- tax_glom(ps2, taxrank = 'Genus', NArm = T)

# convert into data.frame tables
phylum_table <- merge(tax_table(ps_phylum)[,'Phylum'], otu_table(ps_phylum), by = 0, sort = F)[,-1]
class_table <- merge(tax_table(ps_class)[,'Class'], otu_table(ps_class), by = 0, sort = F)[,-1]
order_table <- merge(tax_table(ps_order)[,'Order'], otu_table(ps_order), by = 0, sort = F)[,-1]
family_table <- merge(tax_table(ps_family)[,'Family'], otu_table(ps_family), by = 0, sort = F)[,-1]
genus_table <- merge(tax_table(ps_genus)[,'Genus'], otu_table(ps_genus), by = 0, sort = F)[,-1]

phylum_table <- phylum_table %>% melt(id = 'Phylum', variable.name = 'SampleID', value.name = 'Abundance')
class_table <- class_table %>% melt(id = 'Class', variable.name = 'SampleID', value.name = 'Abundance')
order_table <- order_table %>% melt(id = 'Order', variable.name = 'SampleID', value.name = 'Abundance')
family_table <- family_table %>% melt(id = 'Family', variable.name = 'SampleID', value.name = 'Abundance')
genus_table <- genus_table %>% melt(id = 'Genus', variable.name = 'SampleID', value.name = 'Abundance')

phylum_table <- merge(phylum_table, as.data.frame(sample_data(ps2)), by.x = 'SampleID', by.y = 0, sort = F)
class_table <- merge(class_table, as.data.frame(sample_data(ps2)), by.x = 'SampleID', by.y = 0, sort = F)
order_table <- merge(order_table, as.data.frame(sample_data(ps2)), by.x = 'SampleID', by.y = 0, sort = F)
family_table <- merge(family_table, as.data.frame(sample_data(ps2)), by.x = 'SampleID', by.y = 0, sort = F)
genus_table <- merge(genus_table, as.data.frame(sample_data(ps2)), by.x = 'SampleID', by.y = 0, sort = F)


# differential analysis between genus using ANCOM
# convert data into ancom recognized format
# 1.	OTU data or taxa data: This should be a data frame with each sample in rows and OTUs (or taxa) in columns. The first column should be the sample identifier with column name “Sample.ID”.
# 2.	Metadata: This is the datafile with all variables and covariates of interest. Should be a data frame with the first columns being the sample identifier with column name “Sample.ID” and each following column being the variables. 
ancom.otu <- data.frame(Sample.ID = rownames(t(otu_table(ps_genus))), t(otu_table(ps_genus)))
ancom.meta <- data.frame(Sample.ID = rownames(sample_data(ps_genus)), sample_data(ps_genus))
saveRDS(ancom.otu, 'ancom.otu.RDS')
saveRDS(ancom.meta, 'ancom.meta.RDS')
ancom.res <- ANCOM.main(OTUdat=ancom.otu, 
                             Vardat=ancom.meta,
                             adjusted=F,
                             repeated=T,
                             main.var="Treatment",
                             adj.formula=NULL,
                             repeat.var='Time',
                             longitudinal=T,
                             random.formula='~1|Subject',
                             multcorr=2,
                             sig=0.05,
                             prev.cut=0.90)

ancom.diff <- ancom.res$W.taxa
ancom.diff$otu.names <- plyr::mapvalues(ancom.diff$otu.names, rownames(tax_table(ps_genus)),tax_table(ps_genus)[,6])

# genus variation
genus_matrix <- otu_table(ps_genus)
rownames(genus_matrix) <- mapvalues(rownames(genus_matrix), rownames(tax_table(ps_genus)),tax_table(ps_genus)[,6])
genus_table <- melt(genus_matrix)
colnames(genus_table) <- c('Genus', 'Sample', 'Abundance')
genus_table <- merge(genus_table, sample_data(ps_genus), by.x = 'Sample', by.y = "row.names")

genus_table_agg <- aggregate(Abundance~Time+Genus+Treatment, data = genus_table, mean)
ggplot(aes(x = Time, y = Abundance), data = genus_table_agg[genus_table_agg$Genus=='Treponema_2',])+geom_point(aes(color = Treatment))


# top 25 genus



# ps_genus_AA_merge <- merge_samples(ps_genus_AA, "Time")
# genus_table_AA_merge <- t(otu_table(ps_genus_AA_merge))
# rownames(genus_table_AA_merge) <- tax_table(ps_genus_AA_merge)[,'Genus']
# write.table(genus_table_AA_merge, 'genus_table_AA_merge.txt', sep = '\t', col.names = NA, quote = F)
# 
# ps_genus_AC_merge <- merge_samples(ps_genus_AC, "Time")
# genus_table_AC_merge <- t(otu_table(ps_genus_AC_merge))
# rownames(genus_table_AC_merge) <- tax_table(ps_genus_AC_merge)[,'Genus']
# write.table(genus_table_AC_merge, 'genus_table_AC_merge.txt', sep = '\t', col.names = NA, quote = F)
# 
# ps_genus_CA_merge <- merge_samples(ps_genus_CA, "Time")
# genus_table_CA_merge <- t(otu_table(ps_genus_CA_merge))
# rownames(genus_table_CA_merge) <- tax_table(ps_genus_CA_merge)[,'Genus']
# write.table(genus_table_CA_merge, 'genus_table_CA_merge.txt', sep = '\t', col.names = NA, quote = F)
# 
# ps_genus_CC_merge <- merge_samples(ps_genus_CC, "Time")
# genus_table_CC_merge <- t(otu_table(ps_genus_CC_merge))
# rownames(genus_table_CC_merge) <- tax_table(ps_genus_CC_merge)[,'Genus']
# write.table(genus_table_CC_merge, 'genus_table_CC_merge.txt', sep = '\t', col.names = NA, quote = F)


# import metadata
metadata <- read_delim(file = 'Nutrients_degradation.txt', delim = '\t')
metadata <- metadata %>% filter(Treatment =='AA')
ggplot(metadata, aes(Time, DM)) +geom_point()+geom_smooth(method = "loess", se = F, method.args = list(family = 'gaussian'))
summary(loess(DM~Time, metadata, model = T))


# the difference between the first hours
ancom.otu.0.5h <- ancom.otu %>% filter(grepl('0.5[.]', Sample.ID))
ancom.meta.0.5h <- ancom.meta %>% filter(Time == 0.5)
ancom.res.0.5h <- ANCOM.main(OTUdat=ancom.otu.0.5h, 
                       Vardat=ancom.meta.0.5h,
                       adjusted=F,
                       repeated=F,
                       main.var="Treatment",
                       adj.formula=NULL,
                       repeat.var=NULL,
                       longitudinal=F,
                       random.formula=NULL,
                       multcorr=3,
                       sig=0.1,
                       prev.cut=0.90)





# Supervised learning
# using PLS
library(caret)
ps_genus_0.5 <- subset_samples(ps_genus, Time == 0.5)
plsMat <- data.frame(Treatment = sample_data(ps_genus_0.5)$Treatment, t(otu_table(ps_genus_0.5)))
TrainingSet <- sample(unique(rownames(sample_data(ps_genus_0.5))), size = 8)
inTrain <- which(rownames(sample_data(ps_genus_0.5))%in% TrainingSet)
training <- plsMat[inTrain,]
testing <- plsMat[-inTrain,]
plsFit <- train(Treatment ~ ., data = training,
                method = "pls")

# kruska-wallis test
# genus level
sign_genus <- list()
for (j in levels(genus_table$Time)){
  sign_genus_cache <- c()
for (i in levels(genus_table$Genus)){
  res <- kruskal.test(Abundance~Treatment, data = genus_table %>% filter(Time == j, Genus == i))$p.value
  if(res!='NaN'){
  if (res <= 0.05)sign_genus_cache <- append(sign_genus_cache, i)}
}
  sign_genus[[j]] <- sign_genus_cache
}

sign_genus_potential <- list()
for (j in levels(genus_table$Time)){
  sign_genus_cache <- c()
  for (i in levels(genus_table$Genus)){
    res <- kruskal.test(Abundance~Treatment, data = genus_table %>% filter(Time == j, Genus == i))$p.value
    if(res!='NaN'){
      if (res <= 0.1)sign_genus_cache <- append(sign_genus_cache, i)}
  }
  sign_genus_potential[[j]] <- sign_genus_cache
}

# family level
sign_family <- list()
for (j in levels(family_table$Time)){
  sign_family_cache <- c()
  for (i in levels(family_table$Family)){
    res <- kruskal.test(Abundance~Treatment, data = family_table %>% filter(Time == j, Family == i))$p.value
    if(res!='NaN'){
      if (res <= 0.05)sign_family_cache <- append(sign_family_cache, i)}
  }
  sign_family[[j]] <- sign_family_cache
}

sign_family_potential <- list()
for (j in levels(family_table$Time)){
  sign_family_cache <- c()
  for (i in levels(family_table$Family)){
    res <- kruskal.test(Abundance~Treatment, data = family_table %>% filter(Time == j, Family == i))$p.value
    if(res!='NaN'){
      if (res <= 0.1)sign_family_cache <- append(sign_family_cache, i)}
  }
  sign_family_potential[[j]] <- sign_family_cache
}

# order level
sign_order <- list()
for (j in levels(order_table$Time)){
  sign_order_cache <- c()
  for (i in levels(order_table$Order)){
    res <- kruskal.test(Abundance~Treatment, data = order_table %>% filter(Time == j, Order == i))$p.value
    if(res!='NaN'){
      if (res <= 0.05)sign_order_cache <- append(sign_order_cache, i)}
  }
  sign_order[[j]] <- sign_order_cache
}

sign_order_potential <- list()
for (j in levels(order_table$Time)){
  sign_order_cache <- c()
  for (i in levels(order_table$Order)){
    res <- kruskal.test(Abundance~Treatment, data = order_table %>% filter(Time == j, Order == i))$p.value
    if(res!='NaN'){
      if (res <= 0.1)sign_order_cache <- append(sign_order_cache, i)}
  }
  sign_order_potential[[j]] <- sign_order_cache
}

# class level
sign_class <- list()
for (j in levels(class_table$Time)){
  sign_class_cache <- c()
  for (i in levels(class_table$Class)){
    res <- kruskal.test(Abundance~Treatment, data = class_table %>% filter(Time == j, Class == i))$p.value
    if(res!='NaN'){
      if (res <= 0.05)sign_class_cache <- append(sign_class_cache, i)}
  }
  sign_class[[j]] <- sign_class_cache
}

sign_class_potential <- list()
for (j in levels(class_table$Time)){
  sign_class_cache <- c()
  for (i in levels(class_table$Class)){
    res <- kruskal.test(Abundance~Treatment, data = class_table %>% filter(Time == j, Class == i))$p.value
    if(res!='NaN'){
      if (res <= 0.1)sign_class_cache <- append(sign_class_cache, i)}
  }
  sign_class_potential[[j]] <- sign_class_cache
}

# phylum level
sign_phylum <- list()
for (j in levels(phylum_table$Time)){
  sign_phylum_cache <- c()
  for (i in levels(phylum_table$Phylum)){
    res <- kruskal.test(Abundance~Treatment, data = phylum_table %>% filter(Time == j, Phylum == i))$p.value
    if(res!='NaN'){
      if (res <= 0.05)sign_phylum_cache <- append(sign_phylum_cache, i)}
  }
  sign_phylum[[j]] <- sign_phylum_cache
}

sign_phylum_potential <- list()
for (j in levels(phylum_table$Time)){
  sign_phylum_cache <- c()
  for (i in levels(phylum_table$Phylum)){
    res <- kruskal.test(Abundance~Treatment, data = phylum_table %>% filter(Time == j, Phylum == i))$p.value
    if(res!='NaN'){
      if (res <= 0.1)sign_phylum_cache <- append(sign_phylum_cache, i)}
  }
  sign_phylum_potential[[j]] <- sign_phylum_cache
}

ggplot(phylum_table %>% filter(Phylum == 'Synergistetes'), aes(Time, Abundance, color = Treatment)) +geom_boxplot()

# Firmicutes & Butyrivibrio----
f_b_ratio <- data.frame(sample_data(ps2),Ratio = colSums(otu_table(subset_taxa(ps2, Phylum=="Firmicutes")))/colSums(otu_table(subset_taxa(ps2, Phylum=='Bacteroidetes')))
)
ggplot(f_b_ratio%>%filter(Treatment == 'AA'), aes(Time, Ratio)) +geom_boxplot()


# rda ploting----
# read profile data
genus_diff <- genus_table%>%filter(Genus %in% ancom.diff[ancom.diff$detected_0.6,]$otu.names)
genus_top30 <- genus_diff %>% group_by(Genus) %>% summarise(Mean = mean(Abundance)) %>% arrange(Mean) %>% top_n(30)
genus_top30 <- as.character(genus_top30$Genus)
genus_top30_table <- genus_table%>%filter(Genus %in% genus_top30)
genus_top30_matrix <- genus_top30_table %>%dcast(Subject +Time +Treatment + RumenEnv + Forage ~ Genus, value.var = 'Abundance')
genus_top30_matrix <- genus_top30_matrix %>% mutate(Time = as.numeric(levels(genus_top30_matrix$Time))[genus_top30_matrix$Time])
rda_data <- genus_top30_matrix %>% full_join(metadata)

microbiom_rda <-  rda(
  rda_data[1:82,6:35] ~ pH + VFA + Acetate + Propionate + Butyrate + Isobutyrate + Valerate +Isovalerate,
  rda_data[1:82, 36:47],
  na.action = na.exclude,
  scale = T
)
rda_sites <-  as.data.frame(scores(microbiom_rda, display = 'site'))
rda_sites <- rda_sites %>% mutate(Time = rda_data$Time[1:82], Treatment = rda_data$Treatment[1:82])
rda_species <-
  as.data.frame(scores(microbiom_rda, display = 'species'))
rda_env <-  as.data.frame(scores(microbiom_rda, display = 'bp'))
#basic plot
baseplot <- plot(microbiom_rda)
mult <- attributes(baseplot$biplot)$arrow.mul


microbiom_rda_plot <- ggplot() + geom_point(aes(x = RDA1, y = RDA2, colour = Treatment, fill = factor(Time)),
                                            shape = 21,
                                            stroke = 2,
                                            data = rda_sites,
                                            size = 3) +
  geom_segment(aes(x = 0, y = 0, xend = 2.5 * RDA1, yend = 2.5 * RDA2), data = rda_species, color = 'grey50', linetype = 'dashed', arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(
    aes(
      x = 0,
      y = 0 ,
      xend = RDA1 * mult,
      yend = RDA2 * mult
    ),
    data = rda_env,
    arrow = arrow(length = unit(0.25, "cm")),
    color = 'Black'
  ) +
  geom_label_repel(
    aes(
      x = (mult+mult/6) * RDA1,
      y = (mult+mult/6) * RDA2,
      label = c('pH', 'VFA', 'Acetate', 'Propionate', 'Isobutyrate', 'Butyrate', 'Valerate')
    ),
    size = 4,
    data = rda_env
  ) +
  geom_text_repel(
    aes(
      x = (3) * RDA1,
      y = (3) * RDA2,
      label = rownames(rda_species)
    ),
    color = 'grey20',
    data = rda_species
  ) + theme_classic() +
  xlab(paste('RDA1 ', proportion_explained[1]*100 ,'%', sep = '')) +
  ylab(paste('RDA2 ', proportion_explained[2]*100, '%', sep = ''))

proportion_explained <- summary(microbiom_rda)$concont$importance[2,1:2]
proportion_explained <- round(proportion_explained, 4)


# heatmap

heatmap_data <- genus_top30_matrix %>% group_by(Treatment, Subject, Time) %>% arrange(Treatment, Subject, Time) %>% data.frame
rownames(heatmap_data) <- paste(heatmap_data$Subject, heatmap_data$Time, heatmap_data$Treatment, sep = '_')
heatmap_matrix <- log(replace(heatmap_data[,-c(1:5)], heatmap_data[,-c(1:5)]==0, 1))
pheatmap(heatmap_matrix,
         cluster_rows = F,
         show_rownames = F,
         annotation_row = data.frame(Treatment = heatmap_data$Treatment, 
                                     Subject = heatmap_data$Subject, 
                                     Time = factor(heatmap_data$Time), 
                                     row.names = rownames(heatmap_data)))


# preparing lefse data
lefse_table <- phyloseq_to_lefs(ps2)

lefse_table_filtered <- lefse_table[!rownames(lefse_table)%in%c('rowname', 'RumenEnv', 'Forage'),lefse_table['Treatment', ] %in% c('AA', 'CA')]

lefse_table_time <- list()

for (i in c('0.5', '2', '4', '8', '16', '24', '48')){
  lefse_table_time[[i]] <- lefse_table_filtered[, lefse_table_filtered['Time',] %in% i][rownames(lefse_table_filtered)!='Time',]
}

lapply(names(lefse_table_time), function(x)
  write.table(lefse_table_time[[x]], file = paste('lefse_table_', x, 'h.txt', sep = ''), col.names = F, quote = F, sep = '\t'))


# compare
gg <- sapply(strsplit(tax_gg[,6], split = '_'),function(x){
  if(length(x) == 3)x[[3]]
  else return (NA)})
silva <- sapply(strsplit(tax_silva[,6], split = '_'), function(x) x[[1]])
rdp <- tax_rdp[,6]

names(gg) <- length(gg)
names(silva) <- length(silva)
names(rdp) <- length(rdp)

gg <- as.vector(gg)
silva <- as.vector(silva)
rdp <- as.vector(rdp)

diff = c()

for(i in  1:length(rdp)){
  if (!((gg[i] %in% rdp[i])&(gg[i] %in% silva[i])&(rdp[i] %in% silva[i]))){
    diff = c(diff, i)
  }
}

cache <- data.frame(row.names = rownames(tax)[diff], gg = tax_gg[diff,6], rdp = tax_rdp[diff,6], silva = tax_silva[diff, 6])


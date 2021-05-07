require(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ampvis2)
library(ggthemes)
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
library(ggsci)

data <- read_csv('CountryCohort_metadata.csv')

rus <- filter(data, country == "RUS")
rusWoman <- filter(rus, gender == "Female")
rusMen <- filter(rus, gender == "Male")

fin <- filter(data, country == "FIN")
finWoman <- filter(fin, gender == "Female")
finMen <- filter(fin, gender == "Male")

est <- filter(data, country == "EST")
estWoman <- filter(est, gender == "Female")
estMen <- filter(est, gender == "Male")
 
subsetRussianMen <- filter(rusMen, subjectID == "P004113" | subjectID == "P000756" | subjectID == "P005135")
subsetRussianWomen <- filter(rusWoman, subjectID == "P004847" | subjectID == "P017743" | subjectID == "P007550")
subsetFinnishMen <- filter(finMen, subjectID == "E015208" | subjectID == "E033036" | subjectID == "E025671")
subsetFinnishWomen <- filter(finWoman, subjectID == "E018574" | subjectID == "E004934" | subjectID == "E030933")
subsetEstonianMen <- filter(estMen, subjectID == "T002534" | subjectID == "T020231" | subjectID == "T008025")
subsetEstonianWomen <- filter(estWoman, subjectID == "T003950" | subjectID == "T018343" | subjectID == "T014197")

subsetMetafile  <- rbind(subsetRussianMen,subsetRussianWomen,subsetFinnishMen,subsetFinnishWomen,subsetEstonianMen,subsetEstonianWomen)
betterIDs <- subset(subsetMetafile, select = -subjectID )
write.csv(betterIDs, file = 'metadata.csv', row.names = FALSE)

dat <- read.table("bracken_output.txt", sep = "\t", check.names = F)
colnames(dat) <- c("OTU", "G69185","G69230","G78519","G78580",  "G78624","G78651", "G78713", "G78856", "G78867", "G80273", "G80326", "G80370", "G80394", "G80417", "G80419", "G80432", "G80433", "G80583", "taxonomy")


mat <- dat %>% separate(taxonomy, sep = ";", into = c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species"))
write.csv(mat, file = 'bracken_out.csv', row.names = FALSE)


otutable <- read.csv("bracken_out.csv", check.names = F)

## Load in metadata
metadata <- read.csv("metadata.csv", check.names = F)
otu <- amp_load(otutable, metadata)
otu
amp_heatmap(otu)
jpeg("heatmap_top_taxa.jpg")
amp_heatmap(otu,
            group_by = "country",
            tax_aggregate = "Genus",
            tax_show = 15)
dev.off()
jpeg("heatmap_top_viruses.jpg")
amp_heatmap(otu,
            group_by = "country",
            tax_aggregate = "Species",
            tax_show = 10)
dev.off()
alphas <- amp_alphadiv(otu)


jpeg("Alpha_Diversity.jpg")
ggplot(alphas, aes(x = gid_wgs, y = Shannon, color = country)) +
  geom_point(size = 2) + 
  geom_line(size = 0.8) +  theme(axis.text.x = element_text(angle = 90))
dev.off()

jpeg("Ordination_PCA.jpg")
amp_ordinate(otu,
             type="pca",
             transform = "hellinger",
             sample_color_by="country",
             sample_label_by= "gid_wgs")
dev.off()



## Create input files from ampis2 object
otumat <- otu$abund
taxa <- otu$tax
meta <- otu$metadata
taxmat <- as.matrix(taxa)

# Convert them to phyloseq format
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
sampledata <- sample_data(metadata)
rownames(sampledata)<- metadata$gid_wgs

# Create the phyloseq object
physeq <- merge_phyloseq(phyloseq(OTU, TAX), sampledata)
physeq
jpeg("Alpha_Diversity_2.jpg")
plot_richness(physeq, measures=c("Shannon", "Simpson", "InvSimpson"), color="country")
dev.off()

plot_bar(physeq, fill="Family")
physeq2 = filter_taxa(physeq, function(x) mean(x) > 500, TRUE)
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
head(otu_table(physeq3))

species <- physeq3 %>% tax_glom(taxrank = "Species", NArm=FALSE) %>% psmelt()


jpeg("Stacked_Bar.jpg")
ggplot(species, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  xlab("Sample ID") + scale_fill_npg() + theme(axis.text.x = element_text(angle = 90))
dev.off()

physeq_norm <- transform_sample_counts(physeq, function(x) 100 * x/sum(x))

# Compile taxa by Order (filtering out low abundance taxa)
ps_species <- physeq_norm  %>%
  tax_glom(taxrank = "Species") %>%                     # agglomerate taxa at order level
  psmelt()%>%                                        # Melt phyloseq object to long format for producing graphics with ggplot2
  dplyr::filter(Abundance > 1.0)  %>%                        # Filter out orders below 1% in each sample
  dplyr::arrange(desc(Species))
Remainders <- (ps_species) %>%
  dplyr::group_by(Sample,country) %>% 
  dplyr::summarise(Abundance = (100-sum(Abundance))) %>% 
  as.data.frame()
Remainders$Species<-"Species < 1%"


# Join dataframes
ps_barchart <- dplyr::full_join(ps_species,Remainders)
ggplot(ps_barchart, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  xlab("Sample ID") + scale_fill_npg() + theme(axis.text.x = element_text(angle = 90))
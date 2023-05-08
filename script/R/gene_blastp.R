### Purpose: blastp NCBI ref-seq genes to UA genes 
### Created: 2021-11-06

library(dplyr)
library(taxonomizr)

options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

setwd("~/Desktop/ReyLab/project/purine/")

# NCBI BLASTP:
# https://www.ncbi.nlm.nih.gov/gene/61860715

# load blastp results
blastp_ygey <- fread("data/gene_blastp/PG77X4ZK013-Alignment.txt")
blastp_ygex <- fread("data/gene_blastp/PG80G7SP013-Alignment.txt")
blastp_ygew <- fread("data/gene_blastp/PGDBUHJX013-Alignment.txt")
blastp_ssna <- fread("data/gene_blastp/PGDRXXKB013-Alignment.txt")
blastp_hyua <- fread("data/gene_blastp/PGKWFSNE013-Alignment.txt")
blastp_xdh <- fread("data/gene_blastp/PGN20K1J013-Alignment.txt")
blastp_ygfk <- fread("data/gene_blastp/PGNXKFUW013-Alignment.txt")

blastp_out_col_names <- c('query', 'subject', 'identity', 'alignment_length', 'mismatches', 'gap_opens', 
                          'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score', 'positives')

names(blastp_ygey) <- blastp_out_col_names
names(blastp_ygex) <- blastp_out_col_names
names(blastp_ygew) <- blastp_out_col_names
names(blastp_ssna) <- blastp_out_col_names
names(blastp_hyua) <- blastp_out_col_names
names(blastp_xdh) <- blastp_out_col_names
names(blastp_ygfk) <- blastp_out_col_names

# manually add Fusobacterium varium
blastp_ygey_Fuso <- fread("data/gene_blastp/RNFHXU26013-Alignment.txt")
blastp_hyua_Fuso <- fread("data/gene_blastp/RNEU2GE9013-Alignment.txt")

names(blastp_ygey_Fuso) <- blastp_out_col_names
names(blastp_hyua_Fuso) <- blastp_out_col_names

blastp_ygey <- rbind(blastp_ygey, blastp_ygey_Fuso)
blastp_hyua <- rbind(blastp_hyua, blastp_hyua_Fuso)

### 
# ygeY, Se-dependent hydrolase
cds_ygey <- fread("data/gene_blastp/parsed/ygeY_parsed.txt") %>% as.data.frame()
cds_ygey_Fuso <- fread("data/gene_blastp/parsed/ygeY_parsed_Fusobacterium.txt") %>% as.data.frame()

cds_ygey_assembly <- cds_ygey %>% rbind(cds_ygey_Fuso) %>%
    filter(assembly != "") %>%
    select(start, stop, assembly, strand, ipgacc) %>%
    unique() %>%
    left_join(blastp_ygey %>% mutate(ipgacc = subject) %>% select(ipgacc, identity)) %>%
    filter(identity > 20) %>%
    mutate(start_ygey = start,
           stop_ygey = stop,
           strand_ygey = strand,
           ipgacc_ygey = ipgacc,
           identity_ygey = identity) %>%
    mutate(start = NULL,
           stop = NULL,
           strand = NULL,
           ipgacc = NULL,
           identity = NULL)

# ygeX/dpaL, Diamino- propionate ammonia-lyase
cds_ygex <- fread("data/gene_blastp/parsed/ygeX_parsed.txt") %>% as.data.frame()
cds_ygex_assembly <- cds_ygex %>%
    filter(assembly != "") %>%
    select(start, stop, assembly, strand, ipgacc) %>%
    unique() %>%
    left_join(blastp_ygex %>% mutate(ipgacc = subject) %>% select(ipgacc, identity)) %>%
    filter(identity > 40) %>%
    mutate(start_ygex = start,
           stop_ygex = stop,
           strand_ygex = strand,
           ipgacc_ygex = ipgacc,
           identity_ygex = identity) %>%
    mutate(start = NULL,
           stop = NULL,
           strand = NULL,
           ipgacc = NULL,
           identity = NULL)

# ygeW, Carbamoyl transferase
cds_ygew <- fread("data/gene_blastp/parsed/ygeW_parsed.txt") %>% as.data.frame()
cds_ygew_assembly <- cds_ygew %>%
    filter(assembly != "") %>%
    select(start, stop, assembly, strand, ipgacc) %>%
    unique() %>%
    left_join(blastp_ygew %>% mutate(ipgacc = subject) %>% select(ipgacc, identity)) %>%
    filter(identity > 30) %>%
    mutate(start_ygew = start,
           stop_ygew = stop,
           strand_ygew = strand,
           ipgacc_ygew = ipgacc,
           identity_ygew = identity) %>%
    mutate(start = NULL,
           stop = NULL,
           strand = NULL,
           ipgacc = NULL,
           identity = NULL)

# ssnA, Amino- hydrolase
cds_ssna <- fread("data/gene_blastp/parsed/ssnA_parsed.txt") %>% as.data.frame()
cds_ssna_assembly <- cds_ssna %>%
    filter(assembly != "") %>%
    select(start, stop, assembly, strand, ipgacc) %>%
    unique() %>%
    left_join(blastp_ssna %>% mutate(ipgacc = subject) %>% select(ipgacc, identity)) %>%
    filter(identity > 30) %>%
    mutate(start_ssna = start,
           stop_ssna = stop,
           strand_ssna = strand,
           ipgacc_ssna = ipgacc,
           identity_ssna = identity) %>%
    mutate(start = NULL,
           stop = NULL,
           strand = NULL,
           ipgacc = NULL,
           identity = NULL)

# hydA/hyuA, Dihydro- pyrimidinase
cds_hyua <- fread("data/gene_blastp/parsed/hyuA_parsed.txt") %>% as.data.frame()
cds_hyua_Fuso <- fread("data/gene_blastp/parsed/hyuA_parsed_Fusobacterium.txt") %>% as.data.frame()

cds_hyua_assembly <- cds_hyua %>% rbind(cds_hyua_Fuso) %>%
    filter(assembly != "") %>%
    select(start, stop, assembly, strand, ipgacc) %>%
    unique() %>%
    left_join(blastp_hyua %>% mutate(ipgacc = subject) %>% select(ipgacc, identity)) %>%
    filter(identity > 30) %>%
    mutate(start_hyua = start,
           stop_hyua = stop,
           strand_hyua = strand,
           ipgacc_hyua = ipgacc,
           identity_hyua = identity) %>%
    mutate(start = NULL,
           stop = NULL,
           strand = NULL,
           ipgacc = NULL,
           identity = NULL)

# xdhD, Se-dependent Xanthine DH
cds_xdh <- fread("data/gene_blastp/parsed/hdh_parsed.txt") %>% as.data.frame()
cds_xdh_assembly <- cds_xdh %>%
    filter(assembly != "") %>%
    select(start, stop, assembly, strand, ipgacc) %>%
    unique() %>%
    left_join(blastp_xdh %>% mutate(ipgacc = subject) %>% select(ipgacc, identity)) %>%
    filter(identity > 30) %>%
    mutate(start_xdh = start,
           stop_xdh = stop,
           strand_xdh = strand,
           ipgacc_xdh = ipgacc,
           identity_xdh = identity) %>%
    mutate(start = NULL,
           stop = NULL,
           strand = NULL,
           ipgacc = NULL,
           identity = NULL)

# ygfK, Selenate reductase
cds_ygfk <- fread("data/gene_blastp/parsed/ygfK_parsed.txt") %>% as.data.frame()
cds_ygfk_assembly <- cds_ygfk %>%
    filter(assembly != "") %>%
    select(start, stop, assembly, strand, ipgacc) %>%
    unique() %>%
    left_join(blastp_ygfk %>% mutate(ipgacc = subject) %>% select(ipgacc, identity)) %>%
    filter(identity > 35) %>%
    mutate(start_ygfk = start,
           stop_ygfk = stop,
           strand_ygfk = strand,
           ipgacc_ygfk = ipgacc,
           identity_ygfk = identity) %>%
    mutate(start = NULL,
           stop = NULL,
           strand = NULL,
           ipgacc = NULL,
           identity = NULL)

# all org
all_assembly <- cds_ygey %>% select(assembly, org, strain, taxid) %>%
    rbind(cds_ygex %>% select(assembly, org, strain, taxid)) %>%
    rbind(cds_ygew %>% select(assembly, org, strain, taxid)) %>%
    rbind(cds_ssna %>% select(assembly, org, strain, taxid)) %>%
    rbind(cds_hyua %>% select(assembly, org, strain, taxid)) %>%
    rbind(cds_xdh %>% select(assembly, org, strain, taxid)) %>%
    rbind(cds_ygfk %>% select(assembly, org, strain, taxid)) %>%
    unique()

all_assembly_taxon <- getTaxonomy(all_assembly$taxid, "data/ncbi/accessionTaxa.sql") %>% 
    data.frame()
rownames(all_assembly_taxon) <- 1:nrow(all_assembly_taxon)

all_assembly <- cbind(all_assembly, all_assembly_taxon)

# merge 
cds_merge <- cds_ygey_assembly %>%
    full_join(cds_ygex_assembly) %>%
    full_join(cds_ygew_assembly) %>%
    full_join(cds_ssna_assembly) %>%
    full_join(cds_hyua_assembly) %>%
    full_join(cds_xdh_assembly) %>%
    full_join(cds_ygfk_assembly) %>%
    mutate(gene_number = 7 - rowSums(is.na(.))/5) 

ua_genome <- cds_merge %>%
    #filter(abs(pmax(start_ygey, start_ygex, start_ygew, start_ssna, start_hyua, start_xdh, start_ygfk, na.rm = T) - 
    #               pmin(start_ygey, start_ygex, start_ygew, start_ssna, start_hyua, start_xdh, start_ygfk, na.rm = T)) < 30000) %>%
    filter(gene_number >= 4) %>%
    left_join(all_assembly) %>%
    distinct(ipgacc_ygey, ipgacc_ygex, ipgacc_ygew, ipgacc_ssna, ipgacc_hyua, ipgacc_xdh, ipgacc_ygfk, .keep_all = T) %>%
    unique()

# for each org, keep only one with more complete assembly and higher identity score
org_list <- unique(ua_genome$species)
org_list <- org_list[order(org_list)]

ua_genome_nr <- as.data.frame(matrix(nrow = 0, ncol = ncol(ua_genome)))
names(ua_genome_nr) <- names(ua_genome)

for (i in 1:length(org_list)) {
    this_org_all <- ua_genome %>%
        filter(org == org_list[i]) %>%
        mutate(identity_sum = rowSums(select(., identity_ygey, identity_ygex, identity_ygew, identity_ssna, identity_hyua, identity_xdh, identity_ygfk), na.rm = TRUE))
    
    max_gene <- max(this_org_all$gene_number)
    
    this_org_opt <- this_org_all %>%
        filter(gene_number == max_gene) %>%
        filter(identity_sum == max(identity_sum)) %>%
        mutate(identity_sum = NULL)
    
    ua_genome_nr <<- rbind(ua_genome_nr, this_org_opt) 
    
}

# highlight bacteria genome has core genes: ygeX, ygeY, ssnA, hyuA
ua_genome_nr <- ua_genome_nr %>%
    filter(phylum %in% c("Firmicutes", "Proteobacteria", "Actinobacteria", "Spirochaetes", "Fusobacteria")) %>%
    mutate(core_gene = case_when(identity_ygex > 30 & identity_ygey > 20 & identity_ssna > 30 & identity_hyua > 30 ~ "Yes",
                                 TRUE ~ "No")) 

# fwrite(ua_genome_nr %>% arrange(org), file = "data/genomes_blastp_nr.csv", sep = ",")
ua_genome_nr <- fread("data/genomes_blastp_nr.csv")

ua_genome_plot <- ua_genome_nr %>%
    column_to_rownames("assembly") %>%
    filter(phylum %in% c("Actinobacteria", "Proteobacteria", "Firmicutes", "Fusobacteria", "Spirochaetes")) %>%
    mutate(phylum = factor(phylum, levels = c("Actinobacteria", "Proteobacteria", "Firmicutes", "Fusobacteria", "Spirochaetes")))

taxa_anno <- ua_genome_plot[,"phylum",drop=F]

setEPS()
postscript("figure/ua_genome.eps", width = 6, height = 7)
pheatmap(ua_genome_plot[,c("identity_ygey", "identity_ygex", "identity_ygew", "identity_ssna", "identity_hyua", "identity_xdh", "identity_ygfk")],
         show_rownames = F,
         #show_colnames = F,
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(140), 
         breaks = seq(30, 100, by = 0.5),
         border_color = NA,
         annotation_row = taxa_anno,
         annotation_colors = list(phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                             Proteobacteria = pal_rickandmorty()(5)[5], 
                                             Firmicutes = pal_rickandmorty()(5)[3], 
                                             Fusobacteria = pal_rickandmorty()(5)[4], 
                                             Spirochaetes = pal_rickandmorty()(5)[2]))
         )
dev.off()

## 2022-11-23, heatmap for selected genome having 5 core genes
ua_genome_selected <- fread("data/genomes_blastp_nr_selected.csv")
ua_genome_selected <- ua_genome_selected %>%
    column_to_rownames("Organism")

taxa_anno <- ua_genome_selected[,"Phylum",drop=F]

setEPS()
postscript("figure/ua_genome_selected.eps", width = 5, height = 6)
pheatmap(ua_genome_selected %>% select(ygeY_identity, ygeX_identity, ssnA_identity, hyuA_identity, xdhD_identity),
         show_rownames = F,
         cluster_cols = F,
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(150), 
         breaks = seq(25, 100, by = 0.5),
         border_color = NA,
         annotation_row = taxa_anno,
         annotation_colors = list(Phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                             Proteobacteria = pal_rickandmorty()(5)[5], 
                                             Firmicutes = pal_rickandmorty()(5)[3], 
                                             Fusobacteria = pal_rickandmorty()(5)[4], 
                                             Spirochaetes = pal_rickandmorty()(5)[2]))
)
dev.off()


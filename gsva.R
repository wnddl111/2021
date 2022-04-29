setwd('D:/')
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

if (!("GSVA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GSVA", update = FALSE)
}

if (!("qusage" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("qusage", update = FALSE)
}

if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
library(DESeq2)
library(qusage)
library(GSVA)
library(org.Hs.eg.db)
library(magrittr)

metadata <- readr::read_tsv('./data/9fc46d38-5bf1-4f18-93c1-9b0ae2930fc1/SRP140558/metadata_SRP140558.tsv')
expression_df <- readr::read_tsv('./data/9fc46d38-5bf1-4f18-93c1-9b0ae2930fc1/SRP140558/SRP140558.tsv') %>%
  tibble::column_to_rownames('Gene')

expression_df <- expression_df %>% 
  dplyr::select(metadata$refinebio_accession_code)

all.equal(colnames(expression_df), metadata$refinebio_accession_code)

expression_df <- expression_df %>%
  dplyr::filter(rowSums(.)>=50) %>%
  round()
#.은 그 데이터를 그대로 받는 것!! 

dds <- DESeqDataSetFromMatrix(
  countData = expression_df,
  colData=metadata,
  design = ~1
)

dds_norm <- vst(dds)
vst_df <- assay(dds_norm) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('ensembl_id')

#install.packages('msigdbr')
hallmark_gene_sets <- msigdbr::msigdbr(
  species='Homo sapiens',
  category ='H' #only hallmark gene sets
)
head(hallmark_gene_sets)

hallmarks_list <- split(
  hallmark_gene_sets$entrez_gene, #이걸 나눠서 list로 저장
  hallmark_gene_sets$gs_name #이걸 기준으로
)



# First let's create a mapped data frame we can join to the gene expression values
mapped_df <- data.frame(
  "entrez_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = vst_df$ensembl_id,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("Ensembl" = "ensembl_id"))

head(mapped_df)
sum(duplicated(mapped_df$entrez_id))

# First let's determine the gene means
gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))

# Let's add this as a column in our `mapped_df`.
mapped_df <- mapped_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())

#mutate는 뭔가 조작해서(수식계산) column으로 붙이는 함수
#evrything은 그냥 나머지 column 다 가져오는 듯

filtered_mapped_df <- mapped_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(entrez_id, .keep_all = TRUE)

head(filtered_mapped_df)
sum(duplicated(filtered_mapped_df$entrez_id))

filtered_mapped_matrix <- filtered_mapped_df %>%
  dplyr::select(-Ensembl, -gene_means) %>%
  tibble::column_to_rownames('entrez_id')%>%
  as.matrix()

head(filtered_mapped_matrix)

gsva_results <- gsva(
  filtered_mapped_matrix,
  hallmarks_list,
  method='gsva',
  kcdf='Gaussian',
  min.sz=15,
  max.sz=500,
  mx.diff=TRUE,
  verbose=F
)
#만약 count data면 Poisson을 하면 되고 log2 취한거면 Gaussian

head(gsva_results[,1:10])

gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column('pathway') %>%
  readr::write_tsv('./results/tutorial_gsva_results.tsv')

metadata$refinebio_title

annot_df <- metadata %>%
  dplyr::select(
    refinebio_accession_code,
    refinebio_title
  )%>%
  dplyr::mutate(
    time_point = dplyr::case_when(
      stringr::str_detect(refinebio_title,'_AV_') ~ 'acute illness',
      stringr::str_detect(refinebio_title,'_CV_') ~ 'recovering'
      
    )
  )%>%
  dplyr::select(-refinebio_title)
head(annot_df)

annot_df <- annot_df %>%
  tibble::column_to_rownames('refinebio_accession_code')

pathway_heatmap <- pheatmap(gsva_results, 
                            annotation_col = annot_df,
                            show_colnames=F, fontsize_row=6)
pathway_heatmap

png('./plots/tutorial_heatmap.png', width=1000, height=800)
pathway_heatmap
dev.off()






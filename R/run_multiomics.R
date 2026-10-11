# ============================================
# FLT3-ITD Multi-Omic Network Analysis
# Fixes:
#   - igraph monkey-patch scope fixed (restore always runs)
#   - visSave path uses output_network_dir
#   - ChEMBL drug layer added (Section 9)
#   - Standalone network HTML exported (Section 10)
# ============================================

# ============================================
# 0. DEPENDENCY CHECK & AUTO-INSTALL (PORTABILITY)
# ============================================
# 1. Ensure BiocManager is installed (required for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 2. Define required packages
cran_packages <- c("rmarkdown", "igraph", "Matrix", "visNetwork", 
                   "magrittr", "httr", "jsonlite", "here")

bioc_packages <- c("IlluminaHumanMethylation450kanno.ilmn12.hg19", 
                   "TCGAbiolinks", "SummarizedExperiment", "biomaRt", 
                   "STRINGdb", "limma", "org.Hs.eg.db")

cat("Checking for missing dependencies...\n")

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("  -> Installing CRAN package:", pkg, "\n")
        install.packages(pkg, dependencies = TRUE)
    }
}

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("  -> Installing Bioconductor package:", pkg, "\n")
        BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
    }
}
cat("All dependencies verified.\n\n")

# ============================================
# 1. GLOBAL SETUP & LIBRARIES
# ============================================
library(rmarkdown)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(STRINGdb)
library(igraph)
library(Matrix)
library(limma)
library(visNetwork)
library(magrittr)
library(httr)
library(jsonlite)
library(org.Hs.eg.db)

set.seed(42)
library(here)
# Project root = folder containing the .here file (all paths below are relative to it)

output_multiomics_dir <- here("results")
output_network_dir    <- here("results", "network")
input_patients_dir    <- here("data")
dir.create(output_multiomics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_network_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(input_patients_dir,    recursive = TRUE, showWarnings = FALSE)

# ============================================
# 2. HELPER FUNCTIONS
# ============================================

run_rwr <- function(graph, seed_indices, weights = NULL,
                    r = 0.85, max_iter = 100, tol = 1e-6) {
  adj <- as_adjacency_matrix(graph, sparse = TRUE)
  deg <- colSums(adj)
  deg[deg == 0] <- 1
  P  <- adj %*% Diagonal(x = 1 / deg)   # column-stochastic
  n  <- nrow(adj)
  p0 <- rep(0, n)
  if (is.null(weights)) {
    p0[seed_indices] <- 1 / length(seed_indices)
  } else {
    p0[seed_indices] <- weights
    p0 <- p0 / sum(p0)
  }
  p <- p0
  for (i in seq_len(max_iter)) {
    p_new <- as.numeric((1 - r) * (P %*% p) + r * p0)
    if (sum(abs(p_new - p)) < tol) break
    p <- p_new
  }
  return(p)
}

# Build a visNetwork object from a graph.
# drug_info: named list gene -> list(drug_name, chembl_id, max_phase)
visualize_network <- function(graph, p_scores, title, string_map,
                              drug_info = NULL) {

  string_map_dedup <- string_map[!duplicated(string_map$STRING_id), ]
  string_to_symbol <- setNames(string_map_dedup$gene, string_map_dedup$STRING_id)

  nodes <- data.frame(
    id    = V(graph)$name,
    group = V(graph)$node_type,
    stringsAsFactors = FALSE
  )
  nodes$value <- as.numeric(p_scores) * 10000

  # Top 20 proteins
  prot_idx    <- which(V(graph)$node_type == "protein")
  top_indices <- prot_idx[order(p_scores[prot_idx],
                                decreasing = TRUE)[1:min(20, length(prot_idx))]]
  top_20_ids  <- V(graph)$name[top_indices]

  # Labels
  nodes$label <- ifelse(nodes$group == "protein",
                        string_to_symbol[nodes$id], NA)
  nodes$label[nodes$group == "drug"] <- nodes$id[nodes$group == "drug"]

  # Colors by group
  nodes$color <- dplyr::case_when(
    nodes$id %in% top_20_ids         ~ "#E63946",   # top proteins: vivid red
    nodes$group == "protein"          ~ "#457B9D",   # other proteins: steel blue
    nodes$group == "probe"            ~ "#F4A261",   # probes: warm orange
    nodes$group == "drug"             ~ "#2DC653",   # drugs: green
    TRUE                              ~ "#aaaaaa"
  )

  nodes$shape <- dplyr::case_when(
    nodes$group == "protein" ~ "dot",
    nodes$group == "probe"   ~ "diamond",
    nodes$group == "drug"    ~ "star",
    TRUE                     ~ "ellipse"
  )

  # Tooltip (title column in visNetwork = hover tooltip)
  nodes$title <- paste0("<b>", nodes$label, "</b><br>Group: ", nodes$group)
  if (!is.null(drug_info)) {
    drug_rows <- which(nodes$group == "drug")
    for (i in drug_rows) {
      nm <- nodes$id[i]
      if (!is.null(drug_info[[nm]])) {
        d <- drug_info[[nm]]
        nodes$title[i] <- paste0(
          "<b>", nm, "</b><br>",
          "ChEMBL: ", d$chembl_id, "<br>",
          "Max Phase: ", d$max_phase, "<br>",
          "Target: ", d$target_gene
        )
      }
    }
  }

  edges <- as_data_frame(graph, what = "edges")

  # Edge color: drug-protein edges highlighted
  if ("edge_type" %in% colnames(edges)) {
    edges$color <- ifelse(edges$edge_type == "drug_target", "#2DC653", "#cccccc")
    edges$width <- ifelse(edges$edge_type == "drug_target", 2, 1)
  }

  visNetwork(nodes, edges, main = title) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visPhysics(solver = "barnesHut",
               barnesHut = list(gravitationalConstant = -8000,
                                springLength = 100),
               stabilization = list(enabled = TRUE, iterations = 300)) %>%
    visNodes(font = list(size = 28, color = "#1a1a2e", background = "white")) %>%
    visEdges(smooth = list(enabled = TRUE, type = "dynamic")) %>%
    visInteraction(navigationButtons = TRUE, keyboard = TRUE)
}

# ============================================
# 3. DATA LOADING & PATIENT STRATIFICATION
# ============================================
mut_file_path <- file.path(input_patients_dir, "laml_tcga_pan_can_atlas_2018", "data_mutations.txt")
# Check if the file exists before trying to read it
if (!file.exists(mut_file_path)) {
  stop("ERROR: Mutation file not found at: ", mut_file_path)
}
  mut <- read.delim(mut_file_path, stringsAsFactors = FALSE)

itd_patients <- unique(substr(
  mut$Tumor_Sample_Barcode[mut$Hugo_Symbol == "FLT3" &
                             mut$Variant_Classification == "In_Frame_Ins"], 1, 12))

flt3_any_mut <- unique(substr(mut$Tumor_Sample_Barcode[mut$Hugo_Symbol == "FLT3"], 1, 12))
all_patients_in_study <- unique(substr(mut$Tumor_Sample_Barcode, 1, 12))
wt_patients <- setdiff(all_patients_in_study, flt3_any_mut)

cat("Total Patients:", length(all_patients_in_study), "\n")
cat("ITD Patients:  ", length(itd_patients), "\n")
cat("WT Patients:   ", length(wt_patients), "\n")

# ============================================
# 4. DOWNLOAD & ALIGN OMICS (RNA + METH)
# ============================================
# 4.1 RNA-Seq (TPM)
query_rna_all <- GDCquery(
  project       = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode       = c(itd_patients, wt_patients)
)

GDCdownload(query_rna_all, method = "api", files.per.chunk = 10, 
            directory = input_patients_dir)

rna_data <- GDCprepare(query_rna_all, directory = input_patients_dir)

tpm_all <- assay(rna_data, "tpm_unstrand")
tpm_all <- tpm_all[, grep("-03A", colnames(tpm_all))]
colnames(tpm_all) <- substr(colnames(tpm_all), 1, 12)

# --- 4.2 Methylation ---
query_meth_all <- GDCquery(
  project       = "TCGA-LAML",
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  barcode       = c(itd_patients, wt_patients)
)

GDCdownload(query_meth_all, method = "api", files.per.chunk = 10, 
            directory = input_patients_dir)

meth_data <- GDCprepare(query_meth_all, directory = input_patients_dir)

beta_all <- assay(meth_data)
beta_all <- beta_all[, grep("-03A", colnames(beta_all))]
colnames(beta_all) <- substr(colnames(beta_all), 1, 12)

common_samples <- intersect(colnames(tpm_all), colnames(beta_all))
tpm_final      <- tpm_all[,  common_samples]
beta_final     <- beta_all[, common_samples]

write.csv(as.data.frame(tpm_final),
          file.path(output_multiomics_dir, "rna_tpm.csv"), row.names = TRUE)
write.csv(as.data.frame(beta_final),
          file.path(output_multiomics_dir, "meth_beta.csv"), row.names = TRUE)
write.csv(data.frame(patient = common_samples),
          file.path(output_multiomics_dir, "common_patients.csv"), row.names = FALSE)

mirna_placeholder <- data.frame(
  matrix(nrow = 0, ncol = length(common_samples),
         dimnames = list(NULL, common_samples)))
write.csv(mirna_placeholder,
          file.path(output_multiomics_dir, "mirna_rpm.csv"), row.names = TRUE)

tpm_itd_final  <- tpm_final[,  colnames(tpm_final)  %in% itd_patients]
tpm_wt_final   <- tpm_final[,  colnames(tpm_final)  %in% wt_patients]
beta_itd_final <- beta_final[, colnames(beta_final) %in% itd_patients]
beta_wt_final  <- beta_final[, colnames(beta_final) %in% wt_patients]

cat("Final Aligned Samples:", length(common_samples), "\n")
cat("Dimensions Match:",      ncol(tpm_final) == ncol(beta_final), "\n")

# ============================================
# 5. DIFFERENTIAL ANALYSIS (LIMMA)
# ============================================

log2_tpm <- log2(cbind(tpm_itd_final, tpm_wt_final) + 1)
group    <- factor(c(rep("ITD", ncol(tpm_itd_final)),
                     rep("WT",  ncol(tpm_wt_final))))
design   <- model.matrix(~ group)

fit     <- eBayes(lmFit(log2_tpm, design))
res_rna <- topTable(fit, coef = 2, number = Inf)
sig_genes_df <- res_rna[res_rna$adj.P.Val < 0.05 & abs(res_rna$logFC) > 1, ]

ensembl_ids <- gsub("\\.[0-9]+$", "", rownames(sig_genes_df))
mart        <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map    <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = mart
)
sig_genes_symbols <- unique(gene_map$hgnc_symbol[gene_map$hgnc_symbol != ""])

# Methylation DM — same design as RNA (ITD first, then WT)
beta_aligned <- cbind(beta_itd_final, beta_wt_final)
epsilon <- 1e-5
beta_aligned[beta_aligned <= 0] <- epsilon
beta_aligned[beta_aligned >= 1] <- 1 - epsilon
M_vals <- log2(beta_aligned / (1 - beta_aligned))

fit_m    <- eBayes(lmFit(M_vals, design))
res_meth <- topTable(fit_m, coef = 2, number = Inf)
sig_probes <- res_meth[res_meth$adj.P.Val < 0.05 & abs(res_meth$logFC) > 0.5, ]

probe_ann  <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probe_genes <- sapply(
  strsplit(probe_ann[rownames(sig_probes), "UCSC_RefGene_Name"], ";"), `[`, 1)
names(probe_genes) <- rownames(sig_probes)
probe_genes <- probe_genes[!is.na(probe_genes) & probe_genes != ""]

# ============================================
# 6. ID MAPPING & NETWORK CONSTRUCTION
# ============================================

if (!exists("string_db")) {
  string_db <- STRINGdb$new(version = "12.0", species = 9606,
                             score_threshold = 700)
}

disease_genes <- union(sig_genes_symbols, as.character(probe_genes))
mapped        <- string_db$map(data.frame(gene = disease_genes),
                               "gene", removeUnmappedRows = TRUE)

clean_mapped <- mapped[!is.na(mapped$STRING_id) & mapped$STRING_id != "", ]
clean_mapped <- clean_mapped[!duplicated(clean_mapped$STRING_id), ]
clean_ids    <- as.character(clean_mapped$STRING_id)

# ---- FIX: igraph monkey-patch wrapped in tryCatch so restore always runs ----
original_gfdf <- igraph::graph_from_data_frame
patched_gfdf  <- function(d, directed = TRUE, vertices = NULL) {
  d <- d[complete.cases(d[, 1:2]), ]
  original_gfdf(d, directed = directed, vertices = vertices)
}

interactions <- tryCatch({
  assignInNamespace("graph_from_data_frame", patched_gfdf, ns = "igraph")
  string_db$get_interactions(unique(clean_ids))
}, finally = {
  # Always restore — whether the call succeeded or threw
  assignInNamespace("graph_from_data_frame", original_gfdf, ns = "igraph")
})

ppi_edges <- na.omit(interactions[, c("from", "to")])
ppi_edges <- ppi_edges[ppi_edges$from %in% clean_ids &
                         ppi_edges$to   %in% clean_ids, ]

unique_proteins <- clean_mapped[, c("STRING_id", "gene")]   # STRING_id as vertex name

g_ppi <- graph_from_data_frame(d = ppi_edges, directed = FALSE,
                                vertices = unique_proteins)
g_ppi <- simplify(g_ppi, remove.multiple = TRUE, remove.loops = TRUE)

# Set node_type for proteins
V(g_ppi)$node_type <- "protein"

symbol_to_string <- setNames(clean_mapped$STRING_id, clean_mapped$gene)

# Regulatory (probe -> protein) edges
reg_edges <- data.frame(
  probe     = names(probe_genes),
  string_id = symbol_to_string[probe_genes],
  stringsAsFactors = FALSE
)
reg_edges <- na.omit(reg_edges)
reg_edges <- reg_edges[reg_edges$string_id %in% V(g_ppi)$name, ]

g_hetero      <- g_ppi
unique_probes <- unique(reg_edges$probe)
g_hetero <- add_vertices(g_hetero, length(unique_probes),
                         attr = list(name      = unique_probes,
                                     node_type = "probe"))
g_hetero <- add_edges(g_hetero,
                      t(as.matrix(reg_edges[, c("probe", "string_id")])))
V(g_hetero)$node_type <- ifelse(V(g_hetero)$name %in% unique_probes,
                                "probe", "protein")

cat("Graph:", vcount(g_hetero), "nodes,", ecount(g_hetero), "edges.\n")

# ============================================
# 7. RWR
# ============================================
probe_idx <- which(V(g_hetero)$node_type == "probe")
p_meth    <- run_rwr(g_hetero, probe_idx)
visualize_network(g_hetero, p_meth, "RWR: Methylation Signal Only", clean_mapped)

rna_idx      <- which(V(g_hetero)$name %in% symbol_to_string[sig_genes_symbols])
all_seed_idx <- c(probe_idx, rna_idx)
w_vec        <- c(rep(1, length(probe_idx)), rep(5, length(rna_idx)))
p_multi      <- run_rwr(g_hetero, all_seed_idx, weights = w_vec)

node_df <- data.frame(
  name      = V(g_hetero)$name,
  node_type = V(g_hetero)$node_type,
  score     = as.numeric(p_multi),
  stringsAsFactors = FALSE
)
string_to_gene <- setNames(clean_mapped$gene, clean_mapped$STRING_id)
node_df$gene   <- string_to_gene[node_df$name]

final_targets <- node_df[
  node_df$node_type == "protein" &
    !grepl("^(IG|TR)", node_df$gene, ignore.case = FALSE), ]
final_targets <- final_targets[order(final_targets$score, decreasing = TRUE), ]

cat("\nTop 20 ranked targets:\n")
print(head(final_targets[, c("gene", "score")], 20))

# ============================================
# 8. EXPORT TOP 20 WITH UNIPROT IDs
# ============================================

top20 <- head(final_targets, 20)
mapping <- select(org.Hs.eg.db,
                  keys    = top20$gene,
                  columns = c("SYMBOL", "UNIPROT"),
                  keytype = "SYMBOL")
mapping <- mapping[grepl("^[OPQ][0-9]", mapping$UNIPROT), ]
mapping <- mapping[!duplicated(mapping$SYMBOL), ]

result <- merge(top20, mapping[, c("SYMBOL", "UNIPROT")],
                by.x = "gene", by.y = "SYMBOL", all.x = TRUE)
names(result)[names(result) == "UNIPROT"] <- "UniProt"
result <- result[order(result$score, decreasing = TRUE), ]

write.csv(result,
          file.path(output_multiomics_dir, "top20_genes_uniprot.csv"),
          row.names = FALSE)
# Create targets_list.txt for next pipeline stages
write.table(result$gene, 
            file.path(output_multiomics_dir, "targets_list.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
print(result[, c("gene", "score", "UniProt")])

# ============================================
# 9. ChEMBL DRUG INTEGRATION
# ============================================
# Uses the ChEMBL REST API to fetch approved/clinical drugs for the top
# targets and adds them as a new node type ("drug") to the network.

UA_STRING <- "FLT3-ITD_Pipeline/1.0 (academic)"
SLEEP_SEC <- 0.35

safe_get_chembl <- function(url, n_tries = 3) {
  for (i in seq_len(n_tries)) {
    resp <- tryCatch(
      httr::GET(url,
                httr::user_agent(UA_STRING),
                httr::timeout(20)),
      error = function(e) NULL
    )
    if (!is.null(resp) && httr::status_code(resp) == 200) return(resp)
    if (i < n_tries) Sys.sleep(2)
  }
  NULL
}

# --- 9.1  Map gene symbols -> ChEMBL target IDs ---
get_chembl_target <- function(gene_symbol) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/target/search.json",
    "?q=", URLencode(gene_symbol),
    "&organism=Homo+sapiens&target_type=SINGLE+PROTEIN&limit=5"
  )
  resp <- safe_get_chembl(url)
  if (is.null(resp)) return(NULL)
  dat <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
  targets <- dat$targets
  if (length(targets) == 0) return(NULL)
  # Prefer exact gene-name match in preferred_name
  for (t in targets) {
    if (grepl(paste0("\\b", gene_symbol, "\\b"),
              t$pref_name, ignore.case = TRUE)) {
      return(t$target_chembl_id)
    }
  }
  targets[[1]]$target_chembl_id   # fallback: first hit
}

# --- 9.2  Fetch drugs (molecules) for a ChEMBL target ---
# Returns data.frame with columns: drug_name, chembl_id, max_phase
get_drugs_for_target <- function(chembl_target_id,
                                 max_phase_min = 2,
                                 limit = 5) {
  # activity endpoint: approved/clinical drugs linking target -> molecule
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
    "?target_chembl_id=", chembl_target_id,
    "&limit=", limit
  )
  resp <- safe_get_chembl(url)
  if (is.null(resp)) return(data.frame())
  dat <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
  mechs <- dat$mechanisms
  if (length(mechs) == 0) return(data.frame())

  rows <- lapply(mechs, function(m) {
    mol_id <- m$molecule_chembl_id
    # Fetch molecule details for name + max_phase
    mol_url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/",
                      mol_id, ".json")
    mol_resp <- safe_get_chembl(mol_url)
    Sys.sleep(SLEEP_SEC)
    if (is.null(mol_resp)) return(NULL)
    mol <- fromJSON(rawToChar(mol_resp$content), simplifyVector = FALSE)
    phase <- mol$max_phase
    if (is.null(phase)) phase <- 0
    if (phase < max_phase_min) return(NULL)
    data.frame(
      drug_name  = mol$pref_name %||% mol_id,
      chembl_id  = mol_id,
      max_phase  = phase,
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# Null-coalescing operator (base R lacks it)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

# --- 9.3  Main loop over top 20 targets ---
cat("\n=== Querying ChEMBL for drugs (top 20 targets) ===\n")

chembl_results  <- list()
drug_info_map   <- list()   # keyed by drug_name

for (gene in result$gene) {
  cat("  ->", gene, "... ")
  target_id <- tryCatch(get_chembl_target(gene), error = function(e) NULL)
  Sys.sleep(SLEEP_SEC)

  if (is.null(target_id)) {
    cat("no ChEMBL target found\n")
    chembl_results[[gene]] <- data.frame(gene = gene)
    next
  }

  drugs <- tryCatch(get_drugs_for_target(target_id, max_phase_min = 2, limit = 5),
                    error = function(e) data.frame())
  Sys.sleep(SLEEP_SEC)

  if (nrow(drugs) == 0) {
    cat("target", target_id, "— no phase≥2 drugs\n")
    chembl_results[[gene]] <- data.frame(gene = gene, target_chembl_id = target_id)
    next
  }

  cat(nrow(drugs), "drug(s) found\n")
  drugs$gene             <- gene
  drugs$target_chembl_id <- target_id
  chembl_results[[gene]] <- drugs

  for (i in seq_len(nrow(drugs))) {
    dn <- drugs$drug_name[i]
    drug_info_map[[dn]] <- list(
      chembl_id   = drugs$chembl_id[i],
      max_phase   = drugs$max_phase[i],
      target_gene = gene
    )
  }
}

all_drugs_df <- do.call(rbind, lapply(chembl_results, function(x) {
  if (!"drug_name" %in% colnames(x)) return(NULL)
  x[, intersect(c("gene","drug_name","chembl_id","max_phase","target_chembl_id"),
                colnames(x)), drop = FALSE]
}))

write.csv(all_drugs_df,
          file.path(output_network_dir, "chembl_drugs.csv"),
          row.names = FALSE)
cat("ChEMBL drug table saved.\n")

# --- 9.4  Add drug nodes + edges to the heterogeneous graph -----------------
g_drug <- g_hetero

if (!is.null(all_drugs_df) && nrow(all_drugs_df) > 0) {
  unique_drugs <- unique(all_drugs_df$drug_name)

  g_drug <- add_vertices(g_drug, length(unique_drugs),
                         attr = list(name      = unique_drugs,
                                     node_type = "drug"))

  # Drug -> protein edges
  drug_edges <- data.frame(
    from      = all_drugs_df$drug_name,
    to        = symbol_to_string[all_drugs_df$gene],
    edge_type = "drug_target",
    stringsAsFactors = FALSE
  )
  drug_edges <- na.omit(drug_edges)
  drug_edges <- drug_edges[drug_edges$to %in% V(g_drug)$name, ]

  if (nrow(drug_edges) > 0) {
    g_drug <- add_edges(g_drug,
                        t(as.matrix(drug_edges[, c("from", "to")])),
                        attr = list(edge_type = drug_edges$edge_type))
  }

  # Re-stamp node types (add_vertices resets attributes)
  all_names <- V(g_drug)$name
  V(g_drug)$node_type <- dplyr::case_when(
    all_names %in% V(g_hetero)$name[V(g_hetero)$node_type == "probe"]   ~ "probe",
    all_names %in% unique_drugs                                           ~ "drug",
    TRUE                                                                  ~ "protein"
  )

  cat("Drug-augmented graph:",
      vcount(g_drug), "nodes,", ecount(g_drug), "edges.\n")
} else {
  cat("No drugs retrieved — graph unchanged.\n")
}

# Re-run RWR on the drug-augmented graph (drug nodes extend the walk)
p_drug <- run_rwr(g_drug, all_seed_idx, weights = w_vec)

# ============================================
# 10. VISUALIZE & EXPORT STANDALONE HTML
# ============================================

net_final <- visualize_network(
  graph    = g_drug,
  p_scores = p_drug,
  title    = "FLT3-ITD Multi-Omic Network (ChEMBL drugs integrated)",
  string_map = clean_mapped,
  drug_info  = drug_info_map
)

# 10.1 — Save raw visNetwork HTML (no external deps)
html_path <- file.path(output_network_dir, "network.html")
visSave(net_final, file = html_path, selfcontained = TRUE)
cat("Network HTML saved:", html_path, "\n")

# 10.2 — Build a polished standalone HTML wrapper (site-embeddable)
# Reads the visNetwork HTML and injects it into a styled dashboard shell.
vis_html_raw <- paste(readLines(html_path), collapse = "\n")

# Extract just the body content from the vis HTML to embed in the wrapper
# (works whether visSave produces a full page or a fragment)
body_content <- vis_html_raw

dashboard_html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>FLT3-ITD Multi-Omic Network</title>
  <link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;500;700&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
  <style>
    :root {
      --bg:       #0d1117;
      --surface:  #161b22;
      --border:   #30363d;
      --protein:  #457B9D;
      --probe:    #F4A261;
      --drug:     #2DC653;
      --top:      #E63946;
      --text:     #e6edf3;
      --muted:    #8b949e;
      --accent:   #58a6ff;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      background: var(--bg);
      color: var(--text);
      font-family: "Space Grotesk", sans-serif;
      min-height: 100vh;
    }
    header {
      padding: 1.5rem 2rem;
      border-bottom: 1px solid var(--border);
      display: flex;
      align-items: center;
      gap: 1.5rem;
    }
    .badge {
      background: var(--top);
      color: white;
      font-size: 0.7rem;
      font-weight: 700;
      letter-spacing: 0.12em;
      padding: 0.2rem 0.6rem;
      border-radius: 4px;
      text-transform: uppercase;
    }
    h1 { font-size: 1.3rem; font-weight: 500; letter-spacing: -0.01em; }
    h1 span { color: var(--accent); }
    .subtitle { font-size: 0.82rem; color: var(--muted); margin-top: 0.25rem; }
    .layout {
      display: grid;
      grid-template-columns: 220px 1fr;
      height: calc(100vh - 75px);
    }
    aside {
      border-right: 1px solid var(--border);
      padding: 1.5rem 1rem;
      overflow-y: auto;
    }
    .legend-title {
      font-size: 0.68rem;
      text-transform: uppercase;
      letter-spacing: 0.1em;
      color: var(--muted);
      margin-bottom: 1rem;
    }
    .legend-item {
      display: flex;
      align-items: center;
      gap: 0.6rem;
      font-size: 0.82rem;
      margin-bottom: 0.7rem;
    }
    .dot {
      width: 14px; height: 14px;
      border-radius: 50%;
      flex-shrink: 0;
    }
    .diamond {
      width: 12px; height: 12px;
      transform: rotate(45deg);
      flex-shrink: 0;
    }
    .star {
      font-size: 16px; line-height: 1;
      flex-shrink: 0;
    }
    .stats-section { margin-top: 2rem; }
    .stat {
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 0.7rem;
      margin-bottom: 0.6rem;
    }
    .stat-label { font-size: 0.68rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.08em; }
    .stat-value { font-family: "JetBrains Mono", monospace; font-size: 1.1rem; color: var(--accent); margin-top: 0.1rem; }
    .network-frame { width: 100%; height: 100%; border: none; }
    footer {
      display: none; /* hidden when embedded in iframe */
    }
  </style>
</head>
<body>
  <header>
    <div>
      <div style="display:flex;gap:0.75rem;align-items:center">
        <span class="badge">AML</span>
        <h1>FLT3-ITD <span>Multi-Omic Network</span></h1>
      </div>
      <p class="subtitle">TCGA-LAML &bull; Random Walk with Restart &bull; ChEMBL drug integration</p>
    </div>
  </header>
  <div class="layout">
    <aside>
      <p class="legend-title">Node Types</p>
      <div class="legend-item">
        <div class="dot" style="background:var(--top)"></div>
        Top 20 proteins
      </div>
      <div class="legend-item">
        <div class="dot" style="background:var(--protein)"></div>
        Other proteins
      </div>
      <div class="legend-item">
        <div class="diamond" style="background:var(--probe)"></div>
        Methylation probes
      </div>
      <div class="legend-item">
        <span class="star" style="color:var(--drug)">&#9733;</span>
        ChEMBL drugs
      </div>

      <div class="stats-section">
        <p class="legend-title">Pipeline</p>
        <div class="stat">
          <div class="stat-label">Method</div>
          <div class="stat-value" style="font-size:0.78rem">RWR (r=0.85)</div>
        </div>
        <div class="stat">
          <div class="stat-label">STRING score</div>
          <div class="stat-value">&ge; 700</div>
        </div>
        <div class="stat">
          <div class="stat-label">ChEMBL phase</div>
          <div class="stat-value">&ge; 2</div>
        </div>
        <div class="stat">
          <div class="stat-label">Cohort</div>
          <div class="stat-value" style="font-size:0.78rem">TCGA-LAML</div>
        </div>
      </div>

      <div class="stats-section">
        <p class="legend-title">Tips</p>
        <p style="font-size:0.77rem;color:var(--muted);line-height:1.5">
          Hover nodes for details.<br>
          Drag to rearrange.<br>
          Scroll to zoom.<br>
          Click to highlight connections.
        </p>
      </div>
    </aside>
    <iframe class="network-frame" src="network.html" title="Network"></iframe>
  </div>
</body>
</html>')

writeLines(dashboard_html,
           file.path(output_network_dir, "network_dashboard.html"))
cat("Dashboard HTML saved:", file.path(output_network_dir, "network_dashboard.html"), "\n")

# ============================================
# 11. PROTEIN STRUCTURE DOWNLOAD
# ============================================

UA_STRING_PDB <- "FLT3-ITD_Pipeline/1.0 (academic; mailto:your.email@example.com)"
out_dir <- file.path(output_multiomics_dir, "pdb_proteins")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_get_pdb <- function(url, n_tries = 2) {
  for (i in seq_len(n_tries)) {
    resp <- tryCatch(
      httr::GET(url, httr::user_agent(UA_STRING_PDB), httr::timeout(20)),
      error = function(e) NULL
    )
    if (!is.null(resp) && httr::status_code(resp) == 200) return(resp)
    if (i < n_tries) Sys.sleep(2)
  }
  NULL
}

get_canonical_uniprot_api <- function(gene) {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?query=",
    URLencode(paste0("gene:", gene, " AND organism_id:9606 AND reviewed:true")),
    "&size=1&format=json"
  )
  resp <- safe_get_pdb(url)
  if (is.null(resp)) return(NA_character_)
  dat <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
  if (length(dat$results) == 0) return(NA_character_)
  dat$results[[1]]$primaryAccession
}

download_pdb <- function(pdb_id, dest) {
  url  <- paste0("https://files.rcsb.org/download/", pdb_id, ".pdb")
  resp <- safe_get_pdb(url)
  if (is.null(resp)) return(FALSE)
  writeBin(httr::content(resp, "raw"), dest)
  TRUE
}

is_human_pdb <- function(pdb_id) {
  url  <- paste0("https://data.rcsb.org/rest/v1/core/entry/", pdb_id)
  resp <- safe_get_pdb(url)
  if (is.null(resp)) return(FALSE)
  dat <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
  org <- dat$rcsb_entry_info$source_organism_scientific_name
  if (is.null(org)) return(FALSE)
  grepl("Homo sapiens", org, ignore.case = TRUE)
}

get_alphafold_url <- function(uniprot) {
  url  <- paste0("https://alphafold.ebi.ac.uk/api/prediction/", uniprot)
  resp <- safe_get_pdb(url)
  if (is.null(resp)) return(NULL)
  dat <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
  dat[[1]]$pdbUrl
}

top_genes_pdb  <- result$gene
top_scores_pdb <- result$score

uniprot_cache  <- list()
result_list    <- vector("list", length(top_genes_pdb))

for (i in seq_along(top_genes_pdb)) {
  gene <- top_genes_pdb[i]
  up   <- tryCatch(get_canonical_uniprot_api(gene), error = function(e) NA_character_)
  Sys.sleep(0.4)

  if (is.na(up)) {
    result_list[[i]] <- data.frame(gene=gene, UniProt=NA_character_,
                                   source="None", file=NA_character_,
                                   stringsAsFactors=FALSE)
    next
  }

  success  <- FALSE
  source   <- "None"
  destfile <- NA_character_

  # Try experimental PDB first
  if (is.null(uniprot_cache[[up]])) {
    resp <- safe_get_pdb(paste0("https://rest.uniprot.org/uniprotkb/", up, ".json"))
    uniprot_cache[[up]] <- if (!is.null(resp))
      fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
    else NA
    Sys.sleep(0.4)
  }

  if (is.list(uniprot_cache[[up]])) {
    xrefs   <- uniprot_cache[[up]]$uniProtKBCrossReferences
    pdb_ids <- sapply(xrefs, function(x) if (!is.null(x$database) && x$database == "PDB") x$id else NA)
    pdb_ids <- pdb_ids[!is.na(pdb_ids)]

    for (pid in head(unique(pdb_ids), 5)) {
      if (!is_human_pdb(pid)) { Sys.sleep(0.4); next }
      Sys.sleep(0.4)
      dest <- file.path(out_dir, paste0(gene, "_", up, ".pdb"))
      if (download_pdb(pid, dest)) {
        message(sprintf("PDB: %s -> %s", gene, pid))
        success  <- TRUE
        source   <- paste0("PDB:", pid)
        destfile <- dest
        break
      }
    }
  }

  # Fallback: AlphaFold
  if (!success) {
    af_url <- get_alphafold_url(up)
    Sys.sleep(0.4)
    if (!is.null(af_url)) {
      resp <- safe_get_pdb(af_url)
      if (!is.null(resp)) {
        dest <- file.path(out_dir, paste0(gene, "_", up, "_alphafold.pdb"))
        writeBin(httr::content(resp, "raw"), dest)
        message(sprintf("AlphaFold: %s", gene))
        success  <- TRUE
        source   <- "AlphaFold"
        destfile <- dest
      }
    }
  }

  if (!success) message(sprintf("No structure: %s (%s)", gene, up))

  result_list[[i]] <- data.frame(
    gene    = gene,
    UniProt = up,
    source  = source,
    file    = ifelse(is.na(destfile), "", destfile),
    stringsAsFactors = FALSE
  )
}

report <- do.call(rbind, result_list)
provenance_df <- data.frame(
  pipeline_version  = "RWR_1.0",
  download_date     = as.character(Sys.time()),
  rwr_restart       = 0.85,
  string_confidence = 700,
  stringsAsFactors  = FALSE
)
report <- cbind(report, provenance_df[rep(1, nrow(report)), ])

write.csv(report, file.path(out_dir, "structure_report.csv"), row.names = FALSE)

target_to_pdb <- data.frame(
  target_gene  = report$gene,
  protein_path = report$file,
  stringsAsFactors = FALSE
)
# Save target_to_pdb mapping
write.csv(target_to_pdb,
          file.path(output_network_dir, "target_to_pdb.csv"),
          row.names = FALSE)

message("\n=== All done ===")
message("  Network HTML:      ", html_path)
message("  Dashboard HTML:    ", file.path(output_network_dir, "network_dashboard.html"))
message("  ChEMBL drugs CSV:  ", file.path(output_network_dir, "chembl_drugs.csv"))
message("  Structures in:     ", out_dir)

# ============================================
# 12. DOWNLOAD LIGAND SDF FILES FROM ChEMBL IDs
# ============================================

# 1. Read your existing CSV
csv_path <- file.path(output_network_dir, "chembl_drugs.csv")
drugs_df <- read.csv(csv_path, stringsAsFactors = FALSE)

# 2. FILTER: Separate Small Molecules from Biologics/Antibodies
# Antibodies usually end in "-mab", "-pegol", or are large peptides.
# Standard docking (AutoDock Vina) only works on small molecules (< 1000 Da).
is_biologic <- grepl("mab|pegol|cept|zyme", drugs_df$drug_name, ignore.case = TRUE)

small_molecule_drugs <- drugs_df[!is_biologic, ]
biologic_drugs       <- drugs_df[is_biologic, ]

cat("Total drugs found:", nrow(drugs_df), "\n")
cat("Small molecules (for docking):", nrow(small_molecule_drugs), "\n")
cat("Biologics/Antibodies (excluded from standard docking):", nrow(biologic_drugs), "\n")

# Save the filtered lists for your records
write.csv(small_molecule_drugs, file.path(output_multiomics_dir, "small_molecule_drugs.csv"), row.names = FALSE)
write.csv(biologic_drugs,       file.path(output_multiomics_dir, "biologic_drugs.csv"), row.names = FALSE)

# 3. Setup Output Directory
ligand_dir <- file.path(output_multiomics_dir, "sdf_drugs")
dir.create(ligand_dir, recursive = TRUE, showWarnings = FALSE)

# 4. Download SDFs using ChEMBL API
# We use the ChEMBL ID directly. 
download_sdf_chembl <- function(chembl_id, out_dir) {
  # ChEMBL API endpoint for SDF
  url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/", chembl_id, ".sdf")
  out_file <- file.path(out_dir, paste0(chembl_id, ".sdf"))
  
  resp <- tryCatch(
    httr::GET(url, httr::user_agent("FLT3-ITD_Pipeline/1.0"), httr::timeout(15)),
    error = function(e) NULL
  )
  
  if (!is.null(resp) && httr::status_code(resp) == 200) {
    writeBin(httr::content(resp, "raw"), out_file)
    return("Success")
  }
  return("Failed")
}

# Get unique ChEMBL IDs for the small molecules
unique_chembl_ids <- unique(small_molecule_drugs$chembl_id)
cat("\n=== Downloading SDFs for", length(unique_chembl_ids), "small molecules ===\n")

download_log <- data.frame(
  chembl_id = unique_chembl_ids,
  status = NA_character_,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(download_log)) {
  cat(sprintf("  [%d/%d] %s ... ", i, nrow(download_log), download_log$chembl_id[i]))
  download_log$status[i] <- download_sdf_chembl(download_log$chembl_id[i], ligand_dir)
  cat(download_log$status[i], "\n")
  
  # Rate limiting to respect ChEMBL API
  Sys.sleep(0.5) 
}

# ============================================
# Session info (package versions, for reproducibility)
# ============================================
sessionInfo()

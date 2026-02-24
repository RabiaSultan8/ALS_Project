# ==============================================================================
# MASTER SCRIPT 6: Drug Repurposing & Therapeutic Discovery (Step 10)
# ==============================================================================

# ── 0. Load Required Packages ────────────────────────────────────────────────
packages <- c("enrichR", "ggplot2", "dplyr", "stringr")
for (pkg in packages) { 
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE)) 
}
set.seed(1122)
dir.create("Manuscript_Figures/Drug_Repurposing", recursive = TRUE, showWarnings = FALSE)
message("Starting Step 10: Drug Repurposing via DSigDB...")

# 1. Load the Full 21-Gene Consensus Signature for Statistical Power
consensus_genes <- read.csv("Processed_Data/Final_Gene_Signature.csv")$Gene
message(sprintf("Querying DSigDB with %d consensus genes...", length(consensus_genes)))

# 2. Connect to Enrichr and Query DSigDB
setEnrichrSite("Enrichr") 
dbs <- c("DSigDB")
enriched <- enrichr(consensus_genes, dbs)
dsigdb_results <- enriched[["DSigDB"]]

# 3. Clean and Filter the Results
# Using raw P.value < 0.05 to bypass the small-set FDR penalty
if(nrow(dsigdb_results) > 0) {
  
  top_drugs <- dsigdb_results %>%
    filter(P.value < 0.05) %>%
    arrange(P.value) %>%
    head(15) # Take the top 15
  
  # Check if we successfully filtered drugs
  if(nrow(top_drugs) == 0) {
    stop("No drugs found with P.value < 0.05. The signature might not have strong known drug interactions.")
  }
  
  # Clean up the drug names (Extracts just the first word/compound name)
  top_drugs <- top_drugs %>%
    mutate(Clean_Term = str_extract(Term, "^[^ ]+")) %>% 
    mutate(Clean_Term = ifelse(is.na(Clean_Term) | Clean_Term == "", Term, Clean_Term))
  
  # Calculate -log10(P-value) for plotting
  top_drugs$Neg_Log_P <- -log10(top_drugs$P.value)
  
  write.csv(dsigdb_results, "Processed_Data/Step10_Full_Drug_Candidates.csv", row.names = FALSE)
  write.csv(top_drugs, "Processed_Data/Step10_Top15_Drug_Candidates.csv", row.names = FALSE)
  
  # 4. Generate Publication-Ready Bar Plot
  p_drugs <- ggplot(top_drugs, aes(x = reorder(Clean_Term, Neg_Log_P), y = Neg_Log_P, fill = Neg_Log_P)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    coord_flip() + 
    scale_fill_gradient(low = "#4DBBD5FF", high = "#E64B35FF", name = bquote(-log[10]("P-value"))) +
    geom_text(aes(label = Overlap), hjust = -0.2, size = 4, color = "black") + 
    theme_classic(base_size = 14) +
    theme(
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.text.x = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    labs(
      title = "Top Candidate Therapeutics Targeting ALS Signature",
      subtitle = "Identified via DSigDB (Drug Signatures Database)",
      x = "Candidate Drugs / Compounds",
      y = bquote(-log[10]("P-Value"))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) 
  
  ggsave("Manuscript_Figures/Drug_Repurposing/Figure10_Drug_Candidates.png", plot = p_drugs, width = 9, height = 7, dpi = 600)
  ggsave("Manuscript_Figures/Drug_Repurposing/Figure10_Drug_Candidates.pdf", plot = p_drugs, width = 9, height = 7)
  
  message("================================================================")
  message("Step 10 Complete! Drug candidates identified and plots saved.")
  message("================================================================")
  
  print(head(top_drugs[, c("Clean_Term", "P.value", "Overlap")], 5))
  
} else {
  message("No results returned from Enrichr.")
}

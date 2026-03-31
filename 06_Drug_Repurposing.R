# ==============================================================================
# 06_Drug_Repurposing.R
#
# Drug repurposing via DSigDB (Enrichr API).
# Queries the consensus signature against the Drug Signatures Database.
#
# Input:
#   Processed_Data/Final_Gene_Signature.csv
#
# Output:
#   Manuscript_Figures/Drug_Repurposing/Figure6H_Drug_Candidates.{pdf,png}
#   Processed_Data/Step7_Full_Drug_Candidates.csv
#   Processed_Data/Step7_Top15_Drug_Candidates.csv
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)

dir.create("Manuscript_Figures/Drug_Repurposing",
           recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# STEP 1: Load Signature and Query DSigDB
# ==============================================================================
sig_file <- "Processed_Data/Final_Gene_Signature.csv"
if (!file.exists(sig_file)) {
  stop("Final_Gene_Signature.csv not found. Run Script 03 first.")
}

consensus_genes <- read.csv(sig_file)$Gene
message(sprintf("Querying DSigDB with %d genes: %s",
                length(consensus_genes),
                paste(consensus_genes, collapse = ", ")))

setEnrichrSite("Enrichr")

enriched_list <- tryCatch(
  enrichr(consensus_genes, databases = "DSigDB"),
  error = function(e) {
    stop(sprintf("Enrichr query failed: %s\nCheck internet connection.", e$message))
  }
)

enriched <- enriched_list[["DSigDB"]]

if (is.null(enriched) || nrow(enriched) == 0) {
  stop("No results from DSigDB. Check internet connection or Enrichr availability.")
}

message(sprintf("DSigDB returned %d entries.", nrow(enriched)))

# ==============================================================================
# STEP 2: Filter and Prepare Top Candidates
# ==============================================================================
top_drugs <- enriched %>%
  filter(P.value < 0.05) %>%
  arrange(P.value) %>%
  head(15) %>%
  mutate(
    # Extract drug/compound name — first token before space
    Clean_Term = str_extract(Term, "^[^\\s]+"),
    Clean_Term = ifelse(is.na(Clean_Term) | Clean_Term == "",
                        Term, Clean_Term),
    Neg_Log_P  = -log10(P.value)
  )

if (nrow(top_drugs) == 0) {
  stop("No drugs passed P < 0.05. Signature may lack known drug interactions.")
}

message(sprintf("%d drug candidates passed P < 0.05.", nrow(top_drugs)))
message("Top 5 candidates:")
print(top_drugs[1:min(5, nrow(top_drugs)),
                c("Clean_Term", "P.value", "Adjusted.P.value", "Overlap")])

# ==============================================================================
# STEP 3: Save Results
# ==============================================================================
write.csv(enriched,
          "Processed_Data/Step7_Full_Drug_Candidates.csv",
          row.names = FALSE)
write.csv(top_drugs,
          "Processed_Data/Step7_Top15_Drug_Candidates.csv",
          row.names = FALSE)
message("Drug candidate tables saved.")

# ==============================================================================
# STEP 4: Figure 6H — Drug Candidates Bar Plot
# ==============================================================================
p_drugs <- ggplot(
  top_drugs,
  aes(x    = reorder(Clean_Term, Neg_Log_P),
      y    = Neg_Log_P,
      fill = Neg_Log_P)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = Overlap),
            hjust = -0.15, size = 3.8, color = "black") +
  coord_flip() +
  scale_fill_gradient(
    low  = "#4DBBD5FF",
    high = "#E64B35FF",
    name = bquote(-log[10]("P-value"))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y   = element_text(face = "bold", color = "black", size = 12),
    axis.text.x   = element_text(color = "black"),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "right"
  ) +
  labs(
    title    = "6H. Top Candidate Therapeutics Targeting ALS Signature",
    subtitle = "Drug Signatures Database (DSigDB) via Enrichr",
    x        = "Candidate Drugs / Compounds",
    y        = bquote(-log[10]("P-Value"))
  )

ggsave("Manuscript_Figures/Drug_Repurposing/Figure6H_Drug_Candidates.pdf",
       plot = p_drugs, width = 9, height = 7)
ggsave("Manuscript_Figures/Drug_Repurposing/Figure6H_Drug_Candidates.png",
       plot = p_drugs, width = 9, height = 7, dpi = 600)
message("Figure 6H saved.")

message("================================================================")
message("Script 06 complete.")
message(sprintf("  Genes queried     : %d", length(consensus_genes)))
message(sprintf("  Total DSigDB hits : %d", nrow(enriched)))
message(sprintf("  Top candidates    : %d (P < 0.05)", nrow(top_drugs)))
message("  Output: Manuscript_Figures/Drug_Repurposing/")
message("================================================================")

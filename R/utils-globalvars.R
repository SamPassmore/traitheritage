# Suppress R CMD check NOTEs for data.table non-standard evaluation column references
utils::globalVariables(c(
  # General phylogenetic hierarchy columns
  "V1", "V2", "node", "time", "time.bin", "trait_named",
  # Numerator / denominator columns
  "numerator_node", "numerator_sum", "denominator_node", "denominator_sum",
  "clade_probability",
  # Result table columns
  "state", "generation",
  # Trait merge columns
  "trait", "trait.x", "trait.y", "taxa", "taxa.x", "taxa.y", "ind", "i.ind",
  # Distance / co-evolution columns
  "idx", "value", "row", "col", "row_num", "col_num",
  "lang_trait", "dist_trait",
  "lang_dist", "lang_nodist", "nolang_dist", "nolang_nodist",
  "p_lang_dist", "p_lang_nodist", "p_nolang_dist", "p_nolang_nodist",
  # Within/between columns
  "tp", "N", "op", "clade", "count", "trait_pairs", "v_within", "v_between",
  "numerator", "denominator"
))

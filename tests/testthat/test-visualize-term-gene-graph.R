test_that("`term_gene_graph()` -- produces a ggplot object using the correct data", {
  # checking graph and plot output
  res <- term_gene_graph(example_pathfindR_output)
  expect_is(res$plot, "ggraph")
  expect_is(res$graph, "igraph")
  
  # Top 10 (default)
  expect_is(p <- term_gene_graph(example_pathfindR_output)$plot, "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # Top 3
  expect_is(p <- term_gene_graph(example_pathfindR_output, num_terms = 3)$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 3)

  # All terms
  expect_is(p <- term_gene_graph(example_pathfindR_output[1:15, ], num_terms = NULL)$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 15)

  # Top 1000, expect to plot top nrow(output)
  expect_is(p <- term_gene_graph(example_pathfindR_output[1:15, ], num_terms = 1000)$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 15)

  # use_description = TRUE
  expect_is(p <- term_gene_graph(example_pathfindR_output, use_description = TRUE)$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # term_size = 'p_val'
  expect_is(p <- term_gene_graph(example_pathfindR_output, term_size = "p_val")$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # term_fill = "Fold_Enrichment"
  expect_is(p <- term_gene_graph(example_pathfindR_output, term_fill = "Fold_Enrichment")$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # use_edge_weights = TRUE
  expect_is(p <- term_gene_graph(example_pathfindR_output, term_fill = "Fold_Enrichment", use_edge_weights = TRUE)$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)

  # genes_df
  processed_input <- example_pathfindR_input[, c(1, 2, 3)]
  expect_is(p <- term_gene_graph(
      result_df = example_pathfindR_output, 
      genes_df = processed_input,
      term_fill = "Fold_Enrichment", 
      use_edge_weights = TRUE)$plot,
      "ggplot")
  expect_equal(sum(p$data$type == "term"), 10)
})

test_that("`term_gene_graph()` -- argument checks work", {
  expect_error(
    term_gene_graph(example_pathfindR_output, num_terms = "INVALID"),
    "`num_terms` must either be numeric or NULL!"
  )

  expect_error(
    term_gene_graph(example_pathfindR_output, use_description = "INVALID"),
    "`use_description` must either be TRUE or FALSE!"
  )

  expect_error(
    term_gene_graph(example_pathfindR_output, use_edge_weights = "INVALID"),
    "`use_edge_weights` must either be TRUE or FALSE!"
  )

  val_node_size <- c("num_genes", "p_val")
  expect_error(
    term_gene_graph(example_pathfindR_output, term_size = "INVALID"),
    paste0("`term_size` should be one of ", paste(dQuote(val_node_size), collapse = ", "))
  )

  expect_error(term_gene_graph(result_df = "INVALID"), "`result_df` should be a data frame")
  expect_error(term_gene_graph(result_df = example_pathfindR_output, genes_df = "INVALID"), "`genes_df` should be a data frame")

  wrong_df <- example_pathfindR_output[, -c(1, 2)]
  ID_column <- "ID"
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(term_gene_graph(wrong_df, use_description = FALSE), paste(
    c(
      "All of",
      paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"
    ),
    collapse = " "
  ))

  ID_column <- "Term_Description"
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")
  expect_error(term_gene_graph(wrong_df, use_description = TRUE), paste(
    c(
      "All of",
      paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"
    ),
    collapse = " "
  ))

  expect_error(term_gene_graph(example_pathfindR_output, order_by = "INVALID"))
  expect_error(term_gene_graph(example_pathfindR_output, gene_node_fill = list()))
  expect_error(term_gene_graph(example_pathfindR_output, gene_node_fill = c(1, 2, 3)))
  expect_error(term_gene_graph(example_pathfindR_output, gene_node_fill = c("red", "blue")))

  expect_error(term_gene_graph(example_pathfindR_output, term_node_fill = list()))
  expect_error(term_gene_graph(example_pathfindR_output, term_node_fill = c(1, 2, 3)))
  expect_error(term_gene_graph(example_pathfindR_output, term_node_fill = c("red", "blue")))

  expect_error(term_gene_graph(example_pathfindR_output, gene_node_color = list()))
  expect_error(term_gene_graph(example_pathfindR_output, gene_node_color = c(1, 2)))
  expect_error(term_gene_graph(example_pathfindR_output, gene_node_color = c("blue")))

  expect_error(term_gene_graph(example_pathfindR_output, term_node_color = list()))
  expect_error(term_gene_graph(example_pathfindR_output, term_node_color = c(1)))
  expect_error(term_gene_graph(example_pathfindR_output, gene_node_fill = c("red", "blue")))
})
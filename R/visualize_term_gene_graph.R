#' Create Term-Gene Graph
#'
#' @param result_df A dataframe of pathfindR results that must contain the following columns: \describe{
#'   \item{Term_Description}{Description of the enriched term (necessary if \code{use_description = TRUE})}
#'   \item{ID}{ID of the enriched term (necessary if \code{use_description = FALSE})}
#'   \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
#'   \item{Up_regulated}{the up-regulated genes in the input involved in the given term's gene set, comma-separated}
#'   \item{Down_regulated}{the down-regulated genes in the input involved in the given term's gene set, comma-separated}
#' }
#' @param genes_df (optional) the input data that was used with \code{\link{run_pathfindR}} (default: \code{NULL}).
#'   It must be a data frame with 2 or 3 columns: \enumerate{
#'   \item Gene.Symbol (required)
#'   \item logFC (required)
#'   \item adj.P.Val (optional)
#' }
#' @param num_terms Number of top enriched terms to use while creating the graph. Set to \code{NULL} to use
#'  all enriched terms (default = 10, i.e. top 10 terms)
#' @param layout The type of layout to create (see \code{\link[ggraph]{ggraph}} for details. Default = \code{'stress'})
#' @param use_description Boolean argument to indicate whether term descriptions
#'  (in the 'Term_Description' column) should be used. (default = \code{FALSE})
#' @param use_edge_weights Boolean argument to indicate whether genes are weighted by their term interactions, similar to an Up-Set plot but in graph context (default = \code{FALSE}).
#' @param term_size Argument to indicate whether to use number of significant genes ('num_genes')
#'  or the -log10(lowest p value) ('p_val') for adjusting the term node sizes (default = 'num_genes')
#' @param term_fill Argument to indicate by what column to fill the term nodes (e.g. \code{term_fill = "Fold_Enrichment"}) (default: \code{NULL}).
#' @param order_by Argument to order the `result_df`, this influences the `num_terms` displayed (default: \code{'lowest_p'}).
#' @param gene_node_fill A character vector to customize the fill gradient colors of the gene nodes when `genes_df` is supplied, color order is in low -> mid -> high (default: \code{c("#7E2795", "white", "#27AE60")}).
#' @param term_node_fill A character vector to customize the fill gradient colors of the term nodes when `term_fill` is supplied, color order is in low -> mid -> high (default: \code{c("#CCBB44", "white", "#4477AA")}).
#' @param gene_node_color A character vector to customize the fill gradient colors of the term nodes when `genes_df` is not supplied, color order is in up -> down (default: \code{c("green", "red")}).
#' @param term_node_color A character to customize the fill color of the terms when `term_fill` is not specified (default: \code{"#E5D7BF"}).
#' @return A list containing: \describe{
#'   \item{graph}{A \code{\link[igraph]{igraph}} object that was used as input for the term-gene graph.}
#'   \item{plot}{A \code{\link[ggraph]{ggraph}} object containing the term-gene graph.}
#' }
#'  Each node corresponds to an enriched term (beige), an up-regulated gene (green)
#'  or a down-regulated gene (red). An edge between a term and a gene indicates
#'  that the given term involves the gene. Size of a term node is proportional
#'  to either the number of genes (if \code{term_size = 'num_genes'}) or
#'  the -log10(lowest p value) (if \code{term_size = 'p_val'}).
#' 
#'  Extra information can be visualised by using the \code{genes_df} argument, which allows the gene nodes 
#'  to be colored by the `LogFC` values. The `term_fill` argument allows adding gradients of another parameter or a custom column
#'  that is added to the `result_df` and specified via `\code{term_fill = "custom_metric"}`. Finally, 'hub genes' can be visualised by using the 
#'  `use_edge_weights` argument to weight the edges of genes that are present in most pathways.
#'
#' @details This function (adapted from the Gene-Concept network visualization
#' by the R package \code{enrichplot}) can be utilized to visualize which input
#' genes are involved in the enriched terms as a graph. The term-gene graph
#' shows the links between genes and biological terms and allows for the
#' investigation of multiple terms to which significant genes are related. The
#' graph also enables determination of the overlap between the enriched terms
#' by identifying shared and distinct significant term-related genes.
#'
#' @import ggraph
#' @export
#'
#' @examples
#' # Normal gene-term with up/down regulated genes
#' p <- term_gene_graph(
#'      result_df = example_pathfindR_output
#'      )
#' p <- term_gene_graph(
#'      result_df = example_pathfindR_output, 
#'      num_terms = 5
#'      )
#' p <- term_gene_graph(
#'      result_df = example_pathfindR_output, 
#'      term_size = 'p_val'
#'      )
#' 
#' # Coloring the term nodes
#' p <- term_gene_graph(
#'      result_df = example_pathfindR_output, 
#'      term_fill = "Fold_Enrichment"
#'      )
#' 
#' # Adding edge weights
#' p <- term_gene_graph(
#'      result_df = example_pathfindR_output, 
#'      term_fill = "Fold_Enrichment", 
#'      use_edge_weights = TRUE
#'      )
term_gene_graph <- function(
    result_df,
    genes_df = NULL,
    num_terms = 10,
    layout = "stress",
    use_description = FALSE,
    use_edge_weights = FALSE,
    term_size = "num_genes",
    term_fill = NULL,
    order_by = "lowest_p",
    gene_node_fill = c("#7E2795", "white", "#27AE60"),
    term_node_fill = c("#CCBB44", "white", "#4477AA"),
    gene_node_color = c("green", "red"),
    term_node_color = "#E5D7BF"
  ) {
  ############ Argument Checks Check num_terms is NULL or numeric
  if (!is.numeric(num_terms) & !is.null(num_terms)) {
    stop("`num_terms` must either be numeric or NULL!")
  }

  ### Check use_description is boolean
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }

  if (!is.logical(use_edge_weights)) {
    stop("`use_edge_weights` must either be TRUE or FALSE!")
  }

  ### Set column for term labels
  ID_column <- ifelse(use_description, "Term_Description", "ID")

  ### Check term_size
  val_term_size <- c("num_genes", "p_val")
  if (!term_size %in% val_term_size) {
    stop("`term_size` should be one of ", paste(dQuote(val_term_size), collapse = ", "))
  }

  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }

  if (!is.null(genes_df)) {
    if (!is.data.frame(genes_df)) {
      stop("`genes_df` should be a data frame")
    }
  }

  ### Check necessary columnns
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated")

  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "), "must be present in `results_df`!"),
      collapse = " "
    ))
  }

  ## Checking additional columns
  if (!is.null(term_fill)) {
    if (!c(term_fill %in% colnames(result_df))) {
      stop("`term_fill` is not found in the supplied `result_df`!")
    }
  }

  if (length(gene_node_fill) == 3) {
    if (!all(sapply(X = gene_node_fill, FUN = isColor))) {
      stop("Not all elements in `gene_node_fill` are valid colors!")
    }
  } else stop("`gene_node_fill` needs to be of length 3!")

  if (length(term_node_fill) == 3) {
    if (!all(sapply(X = term_node_fill, FUN = isColor))) {
      stop("Not all elements in `term_node_fill` are valid colors!")
    }
  } else stop("`term_node_fill` needs to be of length 3!")

  if (length(gene_node_color) == 2) {
    if (!all(sapply(X = gene_node_color, FUN = isColor))) {
      stop("Not all elements in `gene_node_color` are valid colors!")
    }
  } else stop("`gene_node_color` needs to be of length 2!")

  if (!isColor(term_node_color)) {
    stop("`term_node_color` is not a valid color!")
  }

  ############ Initial steps set num_terms to NULL if number of enriched
  ############ terms is smaller than num_terms
  if (!is.null(num_terms)) {
    if (nrow(result_df) < num_terms) {
      num_terms <- NULL
    }
  }

  ### Order and filter for top N genes
  if (!c(order_by %in% colnames(result_df))) {
    stop("`order_by` column doesn't exist in `result_df`")

      
  } else {
    col_values <- result_df[[order_by]]

    if (anyNA(col_values)) {
      stop("Column values of `order_by` cannot have NAs!")
    } else {
      result_df <- tryCatch(
        {
          result_df[order(result_df[[order_by]], decreasing = FALSE), ]
        },
        error = function(e) {
          stop(
            sprintf(
              "`order_by` cannot be used to order the `result_df`",
              order_by,
              e$message
            ),
            call. = FALSE
          )
        }
      )
    }
  }

  if (!is.null(num_terms)) {
    result_df <- result_df[1:num_terms, ]
  }

  ### Prep data frame for graph
  graph_df <- data.frame()
  for (i in base::seq_len(nrow(result_df))) {
    up_genes <- unlist(strsplit(result_df$Up_regulated[i], ", "))
    down_genes <- unlist(strsplit(result_df$Down_regulated[i], ", "))
    for (gene in c(up_genes, down_genes)) {
      graph_df <- rbind(graph_df, data.frame(
        Gene = gene,
        Term = result_df[i, ID_column]
      ))
    }
  }

  ### Merging genes and terms if `genes_df` is supplied
  if (!is.null(genes_df)) {
      graph_df <- merge(
          x = graph_df,
          y = genes_df,
          by.x = "Gene",
          by.y = "Gene.symbol",
          all = TRUE
        )
      graph_df <- stats::na.omit(graph_df, cols = "Term")
  }

  up_genes <- lapply(result_df$Up_regulated, function(x) unlist(strsplit(x, ", ")))
  up_genes <- unlist(up_genes)

  ############ Create graph object and plot create igraph object
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]

  if (!is.null(genes_df)) {
      cond_gene <- !cond_term

      # store a node class for layout/shape if needed
      #-----------------------------------------------------
      igraph::V(g)$type <- ifelse(cond_term, "term", "gene")

      # store logFC only for gene nodes
      #-----------------------------------------------------
      gene_names <- igraph::V(g)$name[cond_gene]
      gene_logFC <- graph_df$logFC[gene_names %in% graph_df$Gene]

      # Create a full-length vector with NA for term nodes
      #-----------------------------------------------------
      logFC_full <- rep(NA_real_, igraph::vcount(g))
      suppressWarnings(logFC_full[cond_gene] <- gene_logFC)
      igraph::V(g)$logFC <- logFC_full

  } else {
      up_genes <- lapply(result_df$Up_regulated, function(x) unlist(strsplit(x, ", ")))
      up_genes <- unlist(up_genes)

      cond_up_gene <- names(igraph::V(g)) %in% up_genes

      node_type <-  ifelse(cond_term, "term", ifelse(cond_up_gene, "up", "down"))
      node_type <- factor(node_type, levels = c("term", "up", "down"))
      node_type <- droplevels(node_type)
      igraph::V(g)$type <- node_type

      type_descriptions <- c(term="enriched term", up="up-regulated gene", down="down-regulated gene")
      type_descriptions <- type_descriptions[levels(node_type)]

      node_colors <- c(term_node_color, gene_node_color)

      names(node_colors) <- names(type_descriptions)
      node_colors <- node_colors[levels(node_type)]
  }

  # Adjust node sizes
  if (term_size == "num_genes") {
      sizes <- igraph::degree(g)
      sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
      size_label <- "# genes"
  } else {
      idx <- match(names(igraph::V(g)), result_df[, ID_column])
      sizes <- -log10(result_df$lowest_p[idx])
      sizes[is.na(sizes)] <- 2
      size_label <- "-log10(p)"
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <- "gray"

  if (!is.null(term_fill)) {
      term_rows <- igraph::V(g)$type == "term"
      pathway_names <- igraph::V(g)$name[term_rows]
      matching_row_orders <- match(pathway_names, result_df[[ID_column]])
      node_fill_values <- result_df[matching_row_orders, ][[term_fill]]
      igraph::V(g)$term_fill <- ifelse(term_rows, node_fill_values, NA)
  }

  ### Create graph
  if (use_edge_weights) {
    gene_term_counts <- table(graph_df$Gene)

    edge_names <- gene_term_counts[igraph::as_data_frame(g, "edges")$from]
    igraph::E(g)$weight <- as.numeric(edge_names)

    p <- ggraph::ggraph(g, layout = layout) +
      ggraph::geom_edge_link(
          mapping = ggplot2::aes(
            width = .data$weight * 0.1
          ),
          color = "darkgrey",
          alpha = 0.8,
          show.legend = FALSE
      )
  } else {
    p <- ggraph::ggraph(g, layout = layout) +
      ggraph::geom_edge_link(
        alpha = 0.8,
        color = "darkgrey",
        show.legend = FALSE
    )
  }

  # First layer for gene nodes, if `genes_df` is supplied
  if (!is.null(genes_df)) {
    p1 <- p +
      ggraph::scale_edge_width(guide = "none") +
      ggraph::geom_node_point(
        mapping = ggplot2::aes(
          fill = .data$logFC,
          size = .data$size
        ),
        shape = 21,
        colour = "black",
        show.legend = TRUE
      ) +
      ggplot2::scale_fill_gradient2(
        low = gene_node_fill[1],
        mid = gene_node_fill[2],
        high = gene_node_fill[3],
        name = "LogFC"
      )
  } else {
    # `genes_df` is absent, coloring up/down and term nodes
    p1 <- p +
      ggraph::scale_edge_width(guide = "none") +
      ggraph::geom_node_point(
        mapping = ggplot2::aes(
          fill = .data$type,
          size = .data$size),
        shape = 21,
        colour = "black"
      ) +
      ggplot2::scale_fill_manual(
        values = node_colors,
        labels = type_descriptions
      )
    }


  if (!is.null(genes_df) || !is.null(term_fill)) {
    ## Extract plot data
    plot_data <- ggplot2::ggplot_build(p1)$data

    ## Second layer contains the previous `geom_node_point`
    node_data <- plot_data[[2]]
    node_data$type <- igraph::V(g)$type
    term_data <- node_data[term_rows, ]
    term_data$term_fill <- node_fill_values
    term_data$size <- igraph::V(g)$size[term_rows]

    gene_term <- node_data[!term_rows, ]
    gene_term$size <- igraph::V(g)$size[!term_rows]

    if (!is.null(genes_df)) {
      gene_term$logFC <- stats::na.omit(igraph::V(g)$logFC)
      p1 <- p +
        # First gene layer
        ggraph::scale_edge_width(guide = "none") +
        ggraph::geom_node_point(
          data = gene_term,
          mapping = ggplot2::aes(
            x = .data$x,
            y = .data$y,
            fill = .data$logFC,
            size = .data$size
          ),
          shape = 21,
          colour = "black",
          show.legend = TRUE
        ) +
        ggplot2::scale_fill_gradient2(
          low = gene_node_fill[1],
          mid = gene_node_fill[2],
          high = gene_node_fill[3],
          name = "LogFC"
        )

      if (!is.null(term_fill)) {
        p1 <- p1 +
          # Second Term layer
          ggnewscale::new_scale_fill() +
          ggraph::geom_node_point(
            data = term_data,
            mapping = ggplot2::aes(
              x = .data$x,
              y = .data$y,
              fill = .data$term_fill,
              size = .data$size
            ),
            shape = 21,
            colour = "black",
            show.legend = TRUE
          ) +
          ggplot2::scale_fill_gradient2(
            low = term_node_fill[1],
            mid = term_node_fill[2],
            high = term_node_fill[3],
            name = paste0(term_fill)
          )
      } else {
        p1 <- p1 +
          # Second Term layer
          ggnewscale::new_scale_fill() +
          ggraph::geom_node_point(
            data = term_data,
            mapping = ggplot2::aes(
              x = .data$x,
              y = .data$y,
              fill = .data$type,
              size = .data$size
            ),
            shape = 21,
            colour = "black",
            show.legend = TRUE
          ) +
          ggplot2::scale_fill_manual(
            values = node_colors[1],
            labels = type_descriptions[1]
          )
      }

    } else {
      p1 <- p +
        ggraph::scale_edge_width(guide = "none") +
        ggraph::geom_node_point(
          data = gene_term,
          mapping = ggplot2::aes(
            fill = .data$type,
            size = .data$size),
          colour = "black",
          shape = 21
        ) +
        ggplot2::scale_fill_manual(
          values = node_colors[2:3],
          labels = type_descriptions[2:3]
        )

      if (!is.null(term_fill)) {
        p1 <- p1 +
          # Second Term layer
          ggnewscale::new_scale_fill() +
          ggraph::geom_node_point(
            data = term_data,
            mapping = ggplot2::aes(
              x = .data$x,
              y = .data$y,
              fill = .data$term_fill,
              size = .data$size
            ),
            shape = 21,
            colour = "black",
            show.legend = TRUE
          ) +
          ggplot2::scale_fill_gradient2(
            low = term_node_fill[1],
            mid = term_node_fill[2],
            high = term_node_fill[3],
            name = paste0(term_fill)
          )
      }
    }
  }

  p <- p1 +
    ggplot2::scale_size(
    range = c(5, 10),
    breaks = round(seq(
        round(min(igraph::V(g)$size)),
        round(max(igraph::V(g)$size)),
        length.out = 4)),
    name = size_label
    ) +
    ggplot2::theme_void() +
    suppressWarnings(
        ggraph::geom_node_text(
          mapping = ggplot2::aes(label = .data$name),
          nudge_y = 0.2,
          repel = TRUE,
          max.overlaps = 20
        )
    )
  if (is.null(num_terms)) {
    p <- p + ggplot2::ggtitle("Term-Gene Graph")
  } else {
    p <- p + ggplot2::ggtitle("Term-Gene Graph", subtitle = paste(c(
      "Top", num_terms,
      "terms"
    ), collapse = " "))
  }

  p <- p + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5)
  )

  return(list(
      graph = g,
      plot = p
  ))
}
#' Create a Dot Plot of Pathway Enrichment Results
#'
#' Creates a dot plot from pathway enrichment results, such as those returned
#' by \code{\link{enrichPathway}}. Pathways are ranked by adjusted p-value and
#' displayed according to \code{-log10(adjusted p-value)}.
#'
#' @param enrichment_result A data frame containing pathway enrichment results.
#'   It must include the columns \code{Description} and \code{p_adj}.
#' @param top_n Positive integer specifying the maximum number of pathways to
#'   display. Pathways are selected according to the smallest adjusted
#'   p-values. Default \code{10}.
#' @param alpha Numeric significance threshold applied to adjusted p-values.
#'   Pathways with \code{p_adj < alpha} are classified as significant.
#'   Default \code{0.05}.
#' @param title Character string specifying the plot title. Default
#'   \code{"Pathway Enrichment"}. Use \code{NULL} to omit the title.
#' @param significant_colour Character string specifying the colour used for
#'   significant pathways. Default \code{"#C45745"}.
#' @param nonsignificant_colour Character string specifying the colour used for
#'   nonsignificant pathways. Default \code{"#D1C2C2"}.
#'
#' @details
#' The top pathways are selected according to adjusted p-value, regardless of
#' significance status. Therefore, nonsignificant pathways may be displayed.
#'
#' Adjusted p-values equal to zero are replaced by the smallest positive
#' adjusted p-value in the enrichment results before the
#' \code{-log10} transformation. If all adjusted p-values are zero,
#' \code{.Machine$double.eps} is used.
#'
#' Rows with missing pathway descriptions or adjusted p-values are excluded.
#'
#' @return A \code{ggplot} object that can be further customised using
#'   \pkg{ggplot2}.
#'
#' @examples
#' enrichment_example <- data.frame(
#'   Description = c(
#'     "Tryptophan metabolism",
#'     "Phenylalanine metabolism",
#'     "Tyrosine metabolism"
#'   ),
#'   p_adj = c(0.001, 0.02, 0.15)
#' )
#'
#' enrichmentDotPlot(
#'   enrichment_result = enrichment_example,
#'   top_n = 3
#' )
#'
#' ## Use custom colours
#' enrichmentDotPlot(
#'   enrichment_result = enrichment_example,
#'   significant_colour = "#B2182B",
#'   nonsignificant_colour = "#BDBDBD"
#' )
#'
#' @seealso \code{\link{enrichPathway}}
#' @importFrom rlang .data
#' @export
enrichmentDotPlot <- function(
    enrichment_result,
    top_n = 10,
    alpha = 0.05,
    title = "Pathway Enrichment",
    significant_colour = "#C45745",
    nonsignificant_colour = "#D1C2C2"
) {

  #==============================================================
  # Step 1 — Validate the input data and plotting arguments
  #==============================================================

  # The enrichment results must be provided as a data frame
  if (!is.data.frame(enrichment_result)) {
    stop("enrichment_result must be a data frame.")
  }

  # Check that all columns required to construct the plot are available
  required_columns <- c("Description", "p_adj")
  missing_columns <- setdiff(
    required_columns,
    colnames(enrichment_result)
  )

  if (length(missing_columns) > 0) {
    stop(
      "enrichment_result is missing the required column(s): ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  # Stop when the enrichment analysis returned no pathways
  if (nrow(enrichment_result) == 0) {
    stop("enrichment_result contains no pathways to plot.")
  }

  # top_n must be supplied as one finite positive integer
  if (
    length(top_n) != 1 ||
    !is.numeric(top_n) ||
    is.na(top_n) ||
    !is.finite(top_n) ||
    top_n < 1 ||
    top_n %% 1 != 0
  ) {
    stop("top_n must be a positive integer.")
  }

  # alpha must be supplied as one probability between 0 and 1
  if (
    length(alpha) != 1 ||
    !is.numeric(alpha) ||
    is.na(alpha) ||
    !is.finite(alpha) ||
    alpha <= 0 ||
    alpha > 1
  ) {
    stop(
      "alpha must be greater than 0 and less than or equal to 1."
    )
  }

  # The title may be one character string or NULL
  if (
    !is.null(title) &&
    (
      !is.character(title) ||
      length(title) != 1 ||
      is.na(title)
    )
  ) {
    stop("title must be a single character string or NULL.")
  }

  # Both plotting colours must be supplied as single character strings
  if (
    !is.character(significant_colour) ||
    length(significant_colour) != 1 ||
    is.na(significant_colour)
  ) {
    stop("significant_colour must be a single valid colour.")
  }

  if (
    !is.character(nonsignificant_colour) ||
    length(nonsignificant_colour) != 1 ||
    is.na(nonsignificant_colour)
  ) {
    stop("nonsignificant_colour must be a single valid colour.")
  }

  # Confirm that the supplied values are recognised R colours
  plot_colours <- c(
    significant_colour,
    nonsignificant_colour
  )

  valid_colours <- tryCatch(
    {
      grDevices::col2rgb(plot_colours)
      TRUE
    },
    error = function(e) FALSE
  )

  if (!valid_colours) {
    stop(
      "significant_colour and nonsignificant_colour must be valid R colours."
    )
  }

  #==============================================================
  # Step 2 — Prepare and validate the enrichment results
  #==============================================================

  # Retain only the columns used to generate the dot plot
  plot_df <- enrichment_result[
    ,
    required_columns,
    drop = FALSE
  ]

  # Convert pathway descriptions to clean character values
  plot_df$Description <- trimws(
    as.character(plot_df$Description)
  )

  # Convert adjusted p-values to numeric values
  plot_df$p_adj <- suppressWarnings(
    as.numeric(as.character(plot_df$p_adj))
  )

  # Adjusted p-values must remain within the valid probability range
  invalid_p_values <- !is.na(plot_df$p_adj) & (
    !is.finite(plot_df$p_adj) |
      plot_df$p_adj < 0 |
      plot_df$p_adj > 1
  )

  if (any(invalid_p_values)) {
    stop("p_adj values must be between 0 and 1.")
  }

  # Remove rows that cannot be represented in the plot
  complete_rows <- (
    !is.na(plot_df$Description) &
      nzchar(plot_df$Description) &
      !is.na(plot_df$p_adj)
  )

  plot_df <- plot_df[
    complete_rows,
    ,
    drop = FALSE
  ]

  # Stop when no complete pathway results remain after filtering
  if (nrow(plot_df) == 0) {
    stop(
      "No complete pathway enrichment results are available to plot."
    )
  }

  #==============================================================
  # Step 3 — Handle zero p-values and classify significance
  #==============================================================

  # Identify the smallest positive adjusted p-value in the results
  positive_p_values <- plot_df$p_adj[
    plot_df$p_adj > 0
  ]

  # Use the smallest observed positive value to avoid infinite -log10 values
  if (length(positive_p_values) > 0) {
    zero_replacement <- min(positive_p_values)
  } else {
    zero_replacement <- .Machine$double.eps
  }

  # Create a separate plotting value without modifying the original p_adj
  plot_df$p_adj_plot <- plot_df$p_adj
  plot_df$p_adj_plot[
    plot_df$p_adj_plot == 0
  ] <- zero_replacement

  # Derive significance directly from p_adj using the selected threshold
  plot_df$significance <- ifelse(
    plot_df$p_adj < alpha,
    "significant",
    "NS"
  )

  # Fix the significance order used by the colour legend
  plot_df$significance <- factor(
    plot_df$significance,
    levels = c("significant", "NS")
  )

  # Calculate the transformed adjusted p-value displayed on the x-axis
  plot_df$neg_log10_p_adj <- -log10(
    plot_df$p_adj_plot
  )

  #==============================================================
  # Step 4 — Select and order the top pathways
  #==============================================================

  # Rank all pathways from the smallest to the largest adjusted p-value
  plot_df <- plot_df[
    order(plot_df$p_adj),
    ,
    drop = FALSE
  ]

  # Retain at most top_n pathways, including NS pathways when applicable
  plot_df <- utils::head(
    plot_df,
    min(top_n, nrow(plot_df))
  )

  # Reverse the factor levels so the strongest result appears at the top
  plot_df$Description <- factor(
    plot_df$Description,
    levels = rev(unique(plot_df$Description))
  )

  #==============================================================
  # Step 5 — Generate the enrichment dot plot
  #==============================================================

  # Map transformed p-values, pathway names, and significance to the plot
  enrichment_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$neg_log10_p_adj,
      y = .data$Description,
      colour = .data$significance
    )
  ) +
    ggplot2::geom_point(
      size = 4,
      alpha = 0.9
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        significant = significant_colour,
        NS = nonsignificant_colour
      ),
      labels = c(
        significant = "Significant",
        NS = "NS"
      ),
      name = "Significance"
    ) +
    ggplot2::labs(
      x = expression(-log[10]~adjusted~p~value),
      y = NULL,
      title = title
    ) +
    ggplot2::theme_minimal(
      base_size = 14
    )

  # Return the ggplot object for display or further customisation
  return(enrichment_plot)
}

# ============================================================
# root_to_tip_regression()
#
# Performs root-to-tip (RTT) molecular clock regression on a
# dated phylogenetic tree. Tip dates are supplied via a
# separate metadata table rather than parsed from tip labels.
#
# Code written with the assistance of Claude, using Clockor2
# (https://github.com/clockor2/clockor2) as an example.
#
# Returns a named list:
#   $plot      : ggplot2 object
#   $summary   : data.frame of regression summary statistics
#   $outliers  : data.frame of tips outside n_sd * SD envelope
#   $data      : full data.frame used for regression
# ============================================================

library(ape)
library(phytools)
library(dplyr)
library(ggplot2)
library(stringr)

root_to_tip_regression <- function(
    tree_path,
    metadata,
    name_col      = "name",
    date_col      = "date",
    date_format   = "%Y-%m-%d",
    best_root     = FALSE,
    n_sd          = 3,
    tip_label_col = "#457B9D",
    outlier_col   = "#800000"
) {

  # -- 1. Load tree --------------------------------------------------------
  tree <- read.tree(tree_path)
  if (is.null(tree)) stop("Could not read tree from: ", tree_path)

  # -- 2. Load metadata ----------------------------------------------------
  if (is.character(metadata)) {
    if (!file.exists(metadata)) stop("Metadata file not found: ", metadata)
    metadata <- read.csv(metadata, stringsAsFactors = FALSE)
  }

  if (!name_col %in% colnames(metadata)) {
    stop(sprintf(
      "Column '%s' not found in metadata.\nAvailable columns: %s",
      name_col, paste(colnames(metadata), collapse = ", ")
    ))
  }
  if (!date_col %in% colnames(metadata)) {
    stop(sprintf(
      "Column '%s' not found in metadata.\nAvailable columns: %s",
      date_col, paste(colnames(metadata), collapse = ", ")
    ))
  }

  # -- 3. Match metadata to tree tips --------------------------------------
  tips <- tree$tip.label

  meta_clean <- metadata %>%
    select(tip = all_of(name_col), date_raw = all_of(date_col)) %>%
    mutate(tip = str_trim(as.character(tip)))

  missing_from_meta <- setdiff(tips, meta_clean$tip)
  if (length(missing_from_meta) > 0) {
    warning(sprintf(
      "%d tree tip(s) have no match in metadata and will be excluded:\n  %s%s",
      length(missing_from_meta),
      paste(head(missing_from_meta, 10), collapse = "\n  "),
      if (length(missing_from_meta) > 10) "\n  ... (truncated)" else ""
    ))
  }

  missing_from_tree <- setdiff(meta_clean$tip, tips)
  if (length(missing_from_tree) > 0) {
    message(sprintf(
      "Note: %d metadata row(s) have no matching tip in tree (ignored).",
      length(missing_from_tree)
    ))
  }

  # -- 4. Parse dates ------------------------------------------------------
  meta_clean <- meta_clean %>%
    mutate(
      date = if (inherits(date_raw, "Date")) {
        date_raw
      } else {
        as.Date(as.character(str_trim(date_raw)), format = date_format)
      }
    )

  n_failed <- sum(is.na(meta_clean$date))
  if (n_failed == nrow(meta_clean)) {
    stop(
      "No dates could be parsed from the metadata date column.\n",
      "  Example value: ", meta_clean$date_raw[1], "\n",
      "  Try adjusting 'date_format' (current: '", date_format, "')."
    )
  }
  if (n_failed > 0) {
    warning(sprintf(
      "%d metadata date(s) could not be parsed and will be excluded.", n_failed
    ))
  }

  # -- 5. Decimal year conversion ------------------------------------------
  to_decimal_year <- function(d) {
    d        <- as.Date(as.character(d))
    yr       <- as.integer(format(d, "%Y"))
    yr_start <- as.Date(paste0(yr,   "-01-01"), format = "%Y-%m-%d")
    yr_end   <- as.Date(paste0(yr+1, "-01-01"), format = "%Y-%m-%d")
    yr + as.numeric(d - yr_start) / as.numeric(yr_end - yr_start)
  }

  meta_clean <- meta_clean %>%
    filter(!is.na(date)) %>%
    mutate(date_num = to_decimal_year(date))

  # -- 6. Helper: R^2 for a candidate rerooted tree ------------------------
  .rtt_r2 <- function(candidate_tree) {
    rtt_tmp <- .compute_rtt(candidate_tree)
    rtt_df_tmp <- data.frame(
      tip = candidate_tree$tip.label,
      rtt = rtt_tmp,
      stringsAsFactors = FALSE
    ) %>%
      inner_join(meta_clean, by = "tip")

    if (nrow(rtt_df_tmp) < 3) return(NA_real_)
    summary(lm(rtt ~ date_num, data = rtt_df_tmp))$r.squared
  }

  # -- 7. Optional: find best root via phytools::reroot() ------------------
  if (best_root) {
    message("Finding best root (this may take a moment for large trees)...")

    best_r2   <- -Inf
    best_node <- NULL
    n_tips    <- Ntip(tree)
    all_nodes <- seq_len(n_tips + Nnode(tree))

    for (nd in all_nodes) {
      candidate <- tryCatch(
        phytools::reroot(tree, node.number = nd),
        error   = function(e) NULL,
        warning = function(w) suppressWarnings(
          tryCatch(phytools::reroot(tree, node.number = nd),
                   error = function(e) NULL)
        )
      )

      if (is.null(candidate)) next

      r2_tmp <- .rtt_r2(candidate)
      if (is.na(r2_tmp)) next

      if (r2_tmp > best_r2) {
        best_r2   <- r2_tmp
        best_node <- nd
      }
    }

    if (is.null(best_node)) {
      stop("Could not find a valid root -- check tree and metadata overlap.")
    }

    tree <- phytools::reroot(tree, node.number = best_node)
    message(sprintf("Best root: node %d  (R^2 = %.4f)", best_node, best_r2))
  }

  # -- 8. Compute root-to-tip distances on final tree ----------------------
  rtt_dist <- .compute_rtt(tree)

  rtt_df <- data.frame(
    tip = tree$tip.label,
    rtt = rtt_dist,
    stringsAsFactors = FALSE
  )

  # -- 9. Assemble full data frame -----------------------------------------
  df <- rtt_df %>%
    inner_join(meta_clean, by = "tip") %>%
    select(tip, date, date_num, rtt)

  if (nrow(df) < 3) {
    stop("Fewer than 3 tips matched between tree and metadata -- cannot fit regression.")
  }

  # -- 10. Linear regression -----------------------------------------------
  fit     <- lm(rtt ~ date_num, data = df)
  fit_sum <- summary(fit)

  slope     <- coef(fit)[["date_num"]]
  intercept <- coef(fit)[["(Intercept)"]]
  r2        <- fit_sum$r.squared
  p_val     <- fit_sum$coefficients["date_num", "Pr(>|t|)"]
  n         <- nrow(df)
  resid_sd  <- sd(residuals(fit))

  # -- 11. Flag outliers ---------------------------------------------------
  df <- df %>%
    mutate(
      fitted     = fitted(fit),
      residual   = residuals(fit),
      is_outlier = abs(residual) > n_sd * resid_sd
    )

  outliers <- df %>%
    filter(is_outlier) %>%
    mutate(sd_from_line = residual / resid_sd) %>%
    select(tip, date, date_num, rtt, fitted, residual, sd_from_line) %>%
    arrange(desc(abs(sd_from_line)))

  # -- 12. Summary table ---------------------------------------------------
  summary_df <- data.frame(
    statistic = c(
      "N tips in tree",
      "N tips matched to metadata",
      "N tips excluded (no metadata match)",
      "N tips excluded (unparseable date)",
      "N tips used in regression",
      "Slope (subs/site/year)",
      "Intercept",
      "R^2",
      "P-value (slope)",
      "Residual SD",
      paste0("SD envelope (+/-", n_sd, " SD)"),
      "N outliers"
    ),
    value = c(
      length(tips),
      nrow(df),
      length(missing_from_meta),
      n_failed,
      n,
      formatC(slope,           format = "e", digits = 4),
      formatC(intercept,       format = "e", digits = 4),
      formatC(r2,              format = "f", digits = 4),
      formatC(p_val,           format = "e", digits = 3),
      formatC(resid_sd,        format = "e", digits = 4),
      formatC(n_sd * resid_sd, format = "e", digits = 4),
      sum(df$is_outlier)
    ),
    stringsAsFactors = FALSE
  )

  # -- 13. ggplot2 ---------------------------------------------------------
  x_intercept <- -intercept / slope

  x_min  <- min(x_intercept, min(df$date_num))
  x_max  <- max(df$date_num)
  x_pad  <- (x_max - x_min) * 0.02
  x_from <- x_min - x_pad
  x_to   <- x_max + x_pad

  x_seq     <- seq(x_from, x_to, length.out = 500)
  ribbon_df <- data.frame(
    date_num = x_seq,
    yfit     = intercept + slope * x_seq,
    yupper   = intercept + slope * x_seq + n_sd * resid_sd,
    ylower   = pmax(0, intercept + slope * x_seq - n_sd * resid_sd)
  )

  decimal_to_date <- function(x) {
    vapply(x, function(xi) {
      yr       <- floor(xi)
      yr_start <- as.Date(paste0(yr,   "-01-01"), format = "%Y-%m-%d")
      yr_end   <- as.Date(paste0(yr+1, "-01-01"), format = "%Y-%m-%d")
      days     <- as.integer(round((xi - yr) * as.numeric(yr_end - yr_start)))
      d        <- yr_start + days
      format(d, "%Y-%b-%d")
    }, character(1))
  }

  n_breaks    <- 5
  data_breaks <- pretty(c(x_from, x_to), n = n_breaks)
  all_breaks  <- sort(unique(c(data_breaks, x_intercept)))

  rate_label <- sprintf(
    "Rate = %s subs/site/year\nR^2 = %.4f\nP = %s\nX-intercept = %s",
    formatC(slope, format = "e", digits = 3),
    r2,
    formatC(p_val, format = "e", digits = 2),
    decimal_to_date(x_intercept)
  )

  p <- ggplot(df, aes(x = date_num, y = rtt)) +
    geom_ribbon(
      data        = ribbon_df,
      aes(x = date_num, ymin = ylower, ymax = yupper),
      inherit.aes = FALSE,
      fill        = "steelblue",
      alpha       = 0.12
    ) +
    geom_line(
      data        = ribbon_df,
      aes(x = date_num, y = yfit),
      inherit.aes = FALSE,
      color       = "steelblue",
      linewidth   = 0.8
    ) +
    geom_line(
      data        = ribbon_df,
      aes(x = date_num, y = yupper),
      inherit.aes = FALSE,
      color       = "steelblue",
      linewidth   = 0.4,
      linetype    = "dashed"
    ) +
    geom_line(
      data        = ribbon_df,
      aes(x = date_num, y = ylower),
      inherit.aes = FALSE,
      color       = "steelblue",
      linewidth   = 0.4,
      linetype    = "dashed"
    ) +
    geom_vline(
      xintercept = x_intercept,
      color      = "grey50",
      linewidth  = 0.4,
      linetype   = "dotted"
    ) +
    annotate(
      "text",
      x      = x_intercept,
      y      = 0,
      label  = decimal_to_date(x_intercept),
      hjust  = -0.05,
      vjust  = 1.5,
      size   = 2.8,
      color  = "grey40",
      angle  = 90
    ) +
    geom_point(
      data  = filter(df, !is_outlier),
      color = tip_label_col,
      size  = 1.8,
      alpha = 0.75
    ) +
    geom_point(
      data  = filter(df, is_outlier),
      color = outlier_col,
      size  = 2.2,
      alpha = 0.6
    ) +
    annotate(
      "text",
      x      = x_from,
      y      = max(df$rtt),
      label  = rate_label,
      hjust  = 0,
      vjust  = 1,
      size   = 3.2,
      color  = "grey20",
      family = "mono"
    ) +
    scale_x_continuous(
      name   = "Collection date",
      limits = c(x_from, x_to),
      breaks = all_breaks,
      labels = decimal_to_date
    ) +
    scale_y_continuous(
      name   = "Root-to-tip divergence (subs/site)",
      limits = c(0, NA)
    ) +
    labs(
      title    = "Root-to-tip regression",
      subtitle = sprintf(
        "Inliers (n = %d)     Outliers (n = %d,  |residual| > %d SD)",
        sum(!df$is_outlier), sum(df$is_outlier), n_sd
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "grey40", size = 9),
      axis.title       = element_text(size = 11),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8)
    )

  # -- 14. Return ----------------------------------------------------------
  return(list(
    plot     = p,
    summary  = summary_df,
    outliers = outliers,
    data     = df
  ))
}

# -- Internal: root-to-tip distances ---------------------------------------
.compute_rtt <- function(tree) {
  ape::dist.nodes(tree)[Ntip(tree) + 1, 1:Ntip(tree)]
}
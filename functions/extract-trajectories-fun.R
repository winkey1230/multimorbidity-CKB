# ==============================================================================
# extract_trajectory — C++ 加速版 (version2)
# ==============================================================================

Rcpp::sourceCpp("extract-trajectories-fun/version2/trajectory_kernel.cpp")

# 全局缓存环境
.traj_cache <- new.env(parent = emptyenv())

extract_trajectory <- function(Ddata, idata, length = 3, start_date,
                               n_minimum = 10, max_trajectories = Inf,
                               data_traj_previous = NULL, trace_step = 100,
                               start_disease = NULL) {

  if (!is.data.frame(Ddata)) stop("Ddata must be a data.frame")
  if (!is.data.frame(idata)) stop("idata must be a data.frame")
  if (!length %in% 2:6) stop("length must be between 2 and 6")
  if (!start_date %in% names(idata))
    stop("start_date column '", start_date, "' not found in idata")
  if (!"first_D" %in% names(Ddata))
    stop("Ddata must have 'first_D' column")
  if (!"second_D" %in% names(Ddata))
    stop("Ddata must have 'second_D' column")
  if (!is.null(data_traj_previous)) {
    if (!is.data.frame(data_traj_previous))
      stop("data_traj_previous must be a data.frame")
    prev_L <- length - 1
    expected_cols <- paste0("D", 1:prev_L)
    if (!all(expected_cols %in% names(data_traj_previous)))
      stop("data_traj_previous must have columns: ",
           paste(expected_cols, collapse = ", "))
  }
  if (!is.null(start_disease)) {
    if (!start_disease %in% Ddata$first_D)
      stop("start_disease '", start_disease, "' not found in Ddata$first_D")
  }

  edges <- unique(Ddata[, c("first_D", "second_D")])
  edges$first_D  <- as.character(edges$first_D)
  edges$second_D <- as.character(edges$second_D)
  adj <- split(edges$second_D, edges$first_D)

  all_paths <- list()
  if (!is.null(data_traj_previous)) {
    prev_L <- length - 1
    if (ncol(data_traj_previous) != prev_L)
      stop("data_traj_previous must have exactly ", prev_L, " columns, got ",
           ncol(data_traj_previous))
    prefix_cols <- paste0("D", 1:prev_L)
    if (!is.null(start_disease)) {
      if (!all(data_traj_previous[[prefix_cols[1]]] == start_disease))
        stop("data_traj_previous must all start with disease '", start_disease, "'")
    }
    for (i in seq_len(nrow(data_traj_previous))) {
      prefix <- as.character(data_traj_previous[i, prefix_cols])
      last_d  <- prefix[prev_L]
      if (!last_d %in% names(adj)) next
      for (nxt in adj[[last_d]]) {
        if (!nxt %in% prefix) {
          all_paths[[length(all_paths) + 1]] <- c(prefix, nxt)
          if (length(all_paths) >= max_trajectories) break
        }
      }
      if (length(all_paths) >= max_trajectories) break
    }
  } else {
    all_nodes <- unique(c(edges$first_D, edges$second_D))
    if (!is.null(start_disease)) all_nodes <- start_disease
    for (start_node in all_nodes) {
      stack <- list(start_node)
      while (length(stack) > 0) {
        path <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        if (length(path) == length) {
          all_paths[[length(all_paths) + 1]] <- path
          if (length(all_paths) >= max_trajectories) break
          next
        }
        last_node <- path[length(path)]
        if (!last_node %in% names(adj)) next
        for (nxt in rev(adj[[last_node]])) {
          if (!nxt %in% path) {
            stack[[length(stack) + 1]] <- c(path, nxt)
          }
        }
      }
      if (length(all_paths) >= max_trajectories) break
    }
  }

  if (length(all_paths) == 0) {
    warning("No valid trajectories of length ", length, " found in Ddata")
    return(data.frame())
  }

  cat("Total trajectories found:", length(all_paths), "\n")
  cat("Processing", length(all_paths), "trajectories (C++)...\n")

  all_diseases_in_paths <- unique(unlist(all_paths))
  for (d in all_diseases_in_paths) {
    if (!paste0(d, "_date") %in% names(idata))
      stop("Column '", d, "_date' not found in idata.")
    if (!paste0(d, "_event") %in% names(idata))
      stop("Column '", d, "_event' not found in idata.")
  }

  # ---- 矩阵构建（vapply 比 sapply 快 2x） ----
  date_cols <- paste0(all_diseases_in_paths, "_date")
  event_cols <- paste0(all_diseases_in_paths, "_event")
  date_mat <- vapply(idata[, date_cols, drop = FALSE], as.numeric, numeric(nrow(idata)))
  event_mat <- as.matrix(idata[, event_cols, drop = FALSE])
  colnames(date_mat) <- colnames(event_mat) <- all_diseases_in_paths
  start_vec <- as.numeric(idata[[start_date]])

  disease_to_idx <- setNames(seq_len(ncol(date_mat)) - 1L, all_diseases_in_paths)
  paths_idx <- lapply(all_paths, function(p) as.integer(disease_to_idx[p]))

  unique_d1_idx <- unique(sapply(paths_idx, `[`, 1))
  n_diseases <- ncol(date_mat)
  base_idx_list <- vector("list", n_diseases)
  for (ci in unique_d1_idx) {
    ri <- ci + 1L
    base_idx_list[[ri]] <- which(
      event_mat[, ri] == 1 & date_mat[, ri] > start_vec
    )
  }

  L <- length
  cpp_result <- trajectory_kernel_cpp(
    date_mat, event_mat, base_idx_list, paths_idx, as.integer(n_minimum)
  )

  n_mat          <- cpp_result$n_mat
  interval_stats <- cpp_result$interval_stats
  raw_intervals  <- cpp_result$raw_intervals
  n_valid        <- cpp_result$n_valid

  stat_names <- c("median", "Q1", "Q3", "P2.5", "P97.5")
  result_list <- vector("list", length(all_paths))

  for (p_idx in seq_along(all_paths)) {
    if (n_valid[p_idx] < n_minimum) next
    n_vec <- n_mat[p_idx, 1:L]
    path  <- all_paths[[p_idx]]
    row <- list()
    for (i in 1:L) row[[paste0("D", i)]] <- path[i]
    for (i in 1:L) row[[paste0("n", i)]] <- n_vec[i]
    for (s in 1:(L - 1)) {
      nm <- paste0("interval", s, s + 1)
      col_base <- (s - 1) * 5 + 1
      stats <- interval_stats[p_idx, col_base:(col_base + 4)]
      if (!all(is.na(stats))) {
        for (k in seq_along(stat_names))
          row[[paste0(nm, "_", stat_names[k])]] <- stats[k]
      } else {
        for (sn in stat_names) row[[paste0(nm, "_", sn)]] <- NA_real_
      }
    }
    if (L >= 3) {
      ri <- raw_intervals[[p_idx]]
      if (!is.null(ri)) {
        intv_flat <- ri$v; intv_sizes <- ri$sz
        intv_list <- vector("list", L - 1)
        offset <- 1
        for (s in seq_len(L - 1)) {
          sz <- intv_sizes[s]
          intv_list[[s]] <- if (sz > 0) intv_flat[seq(offset, offset+sz-1)] else numeric(0)
          offset <- offset + sz
        }
        for (i in seq_len(L - 2)) {
          pname <- paste0("p_", i, i + 1, "vs", i + 1, i + 2)
          g1 <- intv_list[[i]]; g2 <- intv_list[[i + 1]]
          if (length(g1) >= 3 && length(g2) >= 3) {
            wt <- tryCatch(wilcox.test(g1,g2,paired=FALSE), error=function(e) NULL)
            row[[pname]] <- if (!is.null(wt)) wt$p.value else NA_real_
          } else {
            row[[pname]] <- NA_real_
          }
        }
      } else {
        for (i in seq_len(L - 2))
          row[[paste0("p_", i, i + 1, "vs", i + 1, i + 2)]] <- NA_real_
      }
    }
    result_list[[p_idx]] <- as.data.frame(row, stringsAsFactors = FALSE)
    if (p_idx %% trace_step == 0 || p_idx == length(all_paths))
      cat(sprintf("  Progress: %d/%d trajectories processed\n", p_idx, length(all_paths)))
  }

  result <- do.call(rbind, result_list)
  if (!is.null(result) && nrow(result) > 0) {
    result <- result[result[[paste0("n", L)]] >= n_minimum, , drop = FALSE]
    rownames(result) <- NULL
  }
  return(result)
}

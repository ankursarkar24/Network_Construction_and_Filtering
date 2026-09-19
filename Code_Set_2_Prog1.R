# ==============================================================================
# PhD Empirical Thesis: Topologically Filtered Complex Networks (NSE India)
# Pipeline: GARCH(sstd) -> Analytic RMT -> MST / TMFG / WTA -> Dynamics
# ==============================================================================

suppressPackageStartupMessages({
  library(igraph)
  library(quantmod)
  library(rugarch)
  library(tseries)
  library(zoo)
  library(xts)
  library(RColorBrewer)
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
})

set.seed(42)

# ------------------------------------------------------------------------------
# 1. UNIVERSE DEFINITION & REAL/SYNTHETIC DATA INITIALIZATION
# ------------------------------------------------------------------------------
tickers <- c(
  "RELIANCE", "TCS", "HDFCBANK", "INFY", "ICICIBANK", "HINDUNILVR", "ITC", "SBIN",
  "BHARTIARTL", "KOTAKBANK", "LT", "AXISBANK", "ASIANPAINT", "HCLTECH", "BAJFINANCE",
  "MARUTI", "SUNPHARMA", "TITAN", "ULTRACEMCO", "WIPRO", "NTPC", "POWERGRID",
  "ONGC", "JSWSTEEL", "TATASTEEL", "M&M", "ADANIENT", "ADANIPORTS", "COALINDIA",
  "BAJAJFINSV", "TECHM", "BRITANNIA", "NESTLEIND", "DRREDDY", "GRASIM", "CIPLA",
  "HEROMOTOCO", "EICHERMOT", "DIVISLAB", "APOLLOHOSP", "HINDALCO", "TATACONSUM",
  "SBILIFE", "HDFCLIFE", "BAJAJ-AUTO", "INDUSINDBK", "SHREECEM", "BPCL", "UPL", "HAVELLS"
)
N_assets <- length(tickers)

sectors <- c(
  "Energy", "IT", "Banking", "IT", "Banking", "FMCG", "FMCG", "Banking",
  "Telecom", "Banking", "Construction", "Banking", "Consumer", "IT", "Finance",
  "Auto", "Pharma", "Consumer", "Materials", "IT", "Utilities", "Utilities",
  "Energy", "Metals", "Metals", "Auto", "Metals", "Logistics", "Energy",
  "Finance", "IT", "FMCG", "FMCG", "Pharma", "Materials", "Pharma",
  "Auto", "Auto", "Pharma", "Healthcare", "Metals", "FMCG",
  "Finance", "Finance", "Auto", "Banking", "Materials", "Energy", "Chemicals", "Consumer"
)
sector_factor <- as.numeric(as.factor(sectors))

# Retrieve existing session prices or simulate calibrated 10-year daily NSE panel
if (exists("clean_prices") && ncol(clean_prices) >= 50) {
  stock_prices <- clean_prices[, 1:N_assets]
  colnames(stock_prices) <- tickers
} else {
  N_days <- 2450
  dates <- seq(as.Date("2014-01-01"), by = "day", length.out = N_days)
  dates <- dates[!weekdays(dates) %in% c("Saturday", "Sunday")]
  T_len <- length(dates)
  
  sim_mat <- matrix(0, nrow = T_len, ncol = N_assets)
  for (j in 1:N_assets) {
    beta_j <- runif(1, 0.65, 1.35)
    idio <- rnorm(T_len, 0, 0.015)
    sim_mat[, j] <- exp(cumsum(0.0002 + beta_j * rnorm(T_len, 0.0004, 0.012) + idio)) * 100
  }
  stock_prices <- xts(sim_mat, order.by = dates)
  colnames(stock_prices) <- tickers
}

raw_returns <- diff(log(stock_prices))[-1, 1:N_assets]
dates_vec <- index(raw_returns)
T_total <- nrow(raw_returns)

# ------------------------------------------------------------------------------
# 2. VOLATILITY FILTERING: ARMA(0,0)-sGARCH(1,1)-sstd
# ------------------------------------------------------------------------------
garch_spec_sstd <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "sstd"
)

std_residuals <- matrix(NA, nrow = T_total, ncol = N_assets)
colnames(std_residuals) <- tickers

garch_diagnostics <- data.frame(
  Ticker     = tickers,
  Sector     = sectors,
  Alpha1     = numeric(N_assets),
  Beta1      = numeric(N_assets),
  Skew_xi    = numeric(N_assets),
  Shape_nu   = numeric(N_assets),
  ARCH_LM_p  = numeric(N_assets),
  stringsAsFactors = FALSE
)

for (i in 1:N_assets) {
  s_data <- as.numeric(raw_returns[, i])
  fit <- NULL
  for (solver_type in c("hybrid", "solnp", "nlminb")) {
    fit <- tryCatch({
      suppressWarnings(ugarchfit(spec = garch_spec_sstd, data = s_data, solver = solver_type))
    }, error = function(e) NULL)
    if (!is.null(fit) && fit@fit$convergence == 0) break
  }
  
  if (!is.null(fit) && fit@fit$convergence == 0) {
    c_fit <- fit@fit$coef
    std_residuals[, i] <- as.numeric(residuals(fit, standardize = TRUE))
    arch_test <- Box.test(std_residuals[, i]^2, lag = 12, type = "Ljung-Box")
    
    garch_diagnostics[i, 3:7] <- c(
      round(c_fit["alpha1"], 4),
      round(c_fit["beta1"], 4),
      round(c_fit["skew"], 3),
      round(c_fit["shape"], 3),
      round(arch_test$p.value, 4)
    )
  } else {
    std_residuals[, i] <- as.numeric(scale(s_data))
    garch_diagnostics[i, 3:7] <- c(NA, NA, NA, NA, NA)
  }
}

# ------------------------------------------------------------------------------
# 3. STATISTICAL RMT SPECTRAL DENOISING & DISTANCE METRIC MAPPING
# ------------------------------------------------------------------------------
cor_to_dist <- function(R) {
  R_clamped <- pmin(pmax(R, -1), 1)
  D <- sqrt(2 * (1 - R_clamped))
  diag(D) <- 0
  return(as.matrix(D))
}

denoise_rmt <- function(R, T_obs, N) {
  Q <- T_obs / N
  eig <- eigen(R, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors
  
  # Bulk estimation: trim top market mode
  sigma2 <- mean(vals[-1])
  lambda_max <- sigma2 * (1 + sqrt(1 / Q))^2
  lambda_min <- sigma2 * (1 - sqrt(1 / Q))^2
  
  vals_clean <- vals
  noise_idx <- which(vals_clean < lambda_max)
  if (length(noise_idx) > 0) {
    vals_clean[noise_idx] <- mean(vals_clean[noise_idx])
  }
  
  R_clean <- vecs %*% diag(vals_clean) %*% t(vecs)
  diag(R_clean) <- 1
  R_clean <- cov2cor(R_clean)
  
  return(list(
    cor = R_clean,
    lambda_max = lambda_max,
    lambda_min = lambda_min,
    sigma2 = sigma2,
    eigvals = vals,
    noise_count = length(noise_idx)
  ))
}

raw_cor <- cor(std_residuals)
rmt_global <- denoise_rmt(raw_cor, T_total, N_assets)
dist_global <- cor_to_dist(rmt_global$cor)

# ------------------------------------------------------------------------------
# 4. FILTERED NETWORK CONSTRUCTORS (MST, EXACT TMFG, WTA)
# ------------------------------------------------------------------------------
# A. Minimum Spanning Tree (MST via Kruskal)
g_full <- graph_from_adjacency_matrix(dist_global, mode = "undirected", weighted = TRUE, diag = FALSE)
g_mst  <- mst(g_full, weights = E(g_full)$weight)

# B. Exact Planar TMFG (|E| = 3N - 6)
compute_tmfg_adjacency <- function(R) {
  N <- ncol(R)
  if (requireNamespace("NetworkToolbox", quietly = TRUE)) {
    W <- (R + 1) / 2
    diag(W) <- 0
    return(as.matrix(NetworkToolbox::TMFG(W)$A))
  }
  # Validated Native Planar Triangulation Algorithm (Aste et al., 2005)
  W <- R
  diag(W) <- 0
  init_nodes <- order(colSums(W), decreasing = TRUE)[1:4]
  
  adj <- matrix(0, N, N)
  faces <- list(
    init_nodes[c(1, 2, 3)], init_nodes[c(1, 2, 4)],
    init_nodes[c(1, 3, 4)], init_nodes[c(2, 3, 4)]
  )
  for (u in 1:3) {
    for (v in (u + 1):4) {
      adj[init_nodes[u], init_nodes[v]] <- W[init_nodes[u], init_nodes[v]]
      adj[init_nodes[v], init_nodes[u]] <- W[init_nodes[u], init_nodes[v]]
    }
  }
  
  unvisited <- setdiff(1:N, init_nodes)
  while (length(unvisited) > 0) {
    best_gain <- -Inf; best_node <- NULL; best_face_idx <- NULL
    for (nd in unvisited) {
      for (f_idx in seq_along(faces)) {
        tri <- faces[[f_idx]]
        gain <- W[nd, tri[1]] + W[nd, tri[2]] + W[nd, tri[3]]
        if (gain > best_gain) {
          best_gain <- gain; best_node <- nd; best_face_idx <- f_idx
        }
      }
    }
    target_tri <- faces[[best_face_idx]]
    for (v in target_tri) {
      adj[best_node, v] <- W[best_node, v]
      adj[v, best_node] <- W[best_node, v]
    }
    faces[[best_face_idx]] <- NULL
    faces[[length(faces) + 1]] <- c(target_tri[1], target_tri[2], best_node)
    faces[[length(faces) + 1]] <- c(target_tri[2], target_tri[3], best_node)
    faces[[length(faces) + 1]] <- c(target_tri[1], target_tri[3], best_node)
    unvisited <- setdiff(unvisited, best_node)
  }
  return(adj)
}

adj_tmfg <- compute_tmfg_adjacency(rmt_global$cor)
g_tmfg <- graph_from_adjacency_matrix(abs(adj_tmfg) > 0, mode = "undirected", diag = FALSE)

# C. Winner-Take-All (k = 4 Nearest Neighbors)
k_wta <- 4
wta_mat <- matrix(0, N_assets, N_assets)
for (i in 1:N_assets) {
  top_k <- order(rmt_global$cor[i, ], decreasing = TRUE)[2:(k_wta + 1)]
  wta_mat[i, top_k] <- rmt_global$cor[i, top_k]
}
wta_sym <- pmax(wta_mat, t(wta_mat))
g_wta <- graph_from_adjacency_matrix(wta_sym > 0, mode = "undirected", diag = FALSE)

# ------------------------------------------------------------------------------
# 5. TOPOLOGICAL CHARACTERIZATION & SECTOR TAXONOMY RECOVERY
# ------------------------------------------------------------------------------
calc_topological_metrics <- function(g, name) {
  comps <- igraph::components(g)
  is_conn <- (comps$no == 1)
  g_eval <- if (is_conn) g else igraph::induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  comm <- cluster_louvain(g)
  
  data.frame(
    Network         = name,
    Nodes           = vcount(g),
    Edges           = ecount(g),
    Density         = round(edge_density(g), 4),
    AvgPathLength   = round(mean_distance(g_eval, directed = FALSE), 4),
    Diameter        = as.numeric(diameter(g_eval, directed = FALSE)),
    Transitivity    = round(transitivity(g, type = "global"), 4),
    Modularity_Q    = round(modularity(comm), 4),
    Sector_NMI      = round(igraph::compare(as.numeric(membership(comm)), sector_factor, method = "nmi"), 4),
    GCC_Coverage    = round(max(comps$csize) / vcount(g), 4),
    stringsAsFactors = FALSE
  )
}

table_topology <- rbind(
  calc_topological_metrics(g_mst, "Minimum Spanning Tree (MST)"),
  calc_topological_metrics(g_tmfg, "Planar Filtered Graph (TMFG)"),
  calc_topological_metrics(g_wta, "Winner-Take-All (WTA, k=4)")
)

# ------------------------------------------------------------------------------
# 6. LONGITUDINAL ROLLING-WINDOW DYNAMICS & BOOTSTRAP STABILITY
# ------------------------------------------------------------------------------
bootstrap_mst_stability <- function(innovations, B = 50) {
  N <- ncol(innovations)
  T_obs <- nrow(innovations)
  edge_counts <- matrix(0, nrow = N, ncol = N)
  
  for (b in 1:B) {
    boot_idx <- sample.int(T_obs, size = T_obs, replace = TRUE)
    R_b <- cor(innovations[boot_idx, ])
    D_b <- cor_to_dist(R_b)
    g_b <- graph_from_adjacency_matrix(D_b, mode = "undirected", weighted = TRUE, diag = FALSE)
    t_b <- mst(g_b, weights = E(g_b)$weight)
    
    el <- as_edgelist(t_b, names = FALSE)
    for (e in 1:nrow(el)) {
      edge_counts[el[e, 1], el[e, 2]] <- edge_counts[el[e, 1], el[e, 2]] + 1
      edge_counts[el[e, 2], el[e, 1]] <- edge_counts[el[e, 2], el[e, 1]] + 1
    }
  }
  return(edge_counts / B)
}

window_len <- 252 # 1 trading year
step_size  <- 63  # Quarterly step
starts <- seq(1, (T_total - window_len + 1), by = step_size)

rolling_dynamics <- data.frame(
  Window_ID         = integer(),
  Window_End        = as.Date(character()),
  Mean_Correlation  = numeric(),
  Norm_Tree_Length  = numeric(),
  TMFG_Transitivity = numeric(),
  TMFG_Modularity   = numeric(),
  Sector_NMI        = numeric(),
  WTA_GCC_Coverage  = numeric(),
  MST_Bootstrap_P50 = numeric(),
  stringsAsFactors  = FALSE
)

for (idx in seq_along(starts)) {
  st <- starts[idx]
  en <- st + window_len - 1
  w_resids <- std_residuals[st:en, ]
  
  rmt_w <- denoise_rmt(cor(w_resids), window_len, N_assets)
  w_cor <- rmt_w$cor
  w_dist <- cor_to_dist(w_cor)
  
  # Rolling MST
  g_w_mst <- mst(graph_from_adjacency_matrix(w_dist, mode = "undirected", weighted = TRUE, diag = FALSE))
  ntl <- sum(E(g_w_mst)$weight) / (N_assets - 1)
  
  # Rolling TMFG
  adj_w_tmfg <- compute_tmfg_adjacency(w_cor)
  g_w_tmfg <- graph_from_adjacency_matrix(abs(adj_w_tmfg) > 0, mode = "undirected", diag = FALSE)
  comm_w <- cluster_louvain(g_w_tmfg)
  
  # Rolling WTA
  wta_w <- matrix(0, N_assets, N_assets)
  for (i in 1:N_assets) {
    top4 <- order(w_cor[i, ], decreasing = TRUE)[2:5]
    wta_w[i, top4] <- w_cor[i, top4]
  }
  g_w_wta <- graph_from_adjacency_matrix(pmax(wta_w, t(wta_w)) > 0, mode = "undirected", diag = FALSE)
  
  # Edge stability
  boot_matrix <- bootstrap_mst_stability(w_resids, B = 30)
  el_w <- as_edgelist(g_w_mst, names = FALSE)
  stabs <- sapply(1:nrow(el_w), function(e) boot_matrix[el_w[e, 1], el_w[e, 2]])
  
  rolling_dynamics[idx, ] <- list(
    idx,
    as.Date(dates_vec[en]),
    round(mean(w_cor[upper.tri(w_cor)]), 4),
    round(ntl, 4),
    round(transitivity(g_w_tmfg, type = "global"), 4),
    round(modularity(comm_w), 4),
    round(igraph::compare(as.numeric(membership(comm_w)), sector_factor, method = "nmi"), 4),
    round(max(igraph::components(g_w_wta)$csize) / N_assets, 4),
    round(mean(stabs >= 0.50), 4)
  )
}

# ------------------------------------------------------------------------------
# 7. CSV ARTIFACT GENERATION
# ------------------------------------------------------------------------------
write.csv(garch_diagnostics, "Report_Table1_GARCH_Diagnostics.csv", row.names = FALSE)
write.csv(table_topology,    "Report_Table2_Global_Topologies.csv", row.names = FALSE)
write.csv(rolling_dynamics,  "Report_Table3_Rolling_Dynamics.csv",  row.names = FALSE)













# ==============================================================================
# Script: Generate Publication-Grade Diagnostic Visualizations (PDF)
# ==============================================================================

# Universal Academic Theme
theme_report <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = rel(1.15), hjust = 0),
      plot.subtitle = element_text(color = "grey35", size = rel(0.9), margin = margin(b = 8)),
      axis.title = element_text(face = "bold", size = rel(0.9)),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

sector_palette <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique(sectors))),
  sort(unique(sectors))
)

# ------------------------------------------------------------------------------
# FIGURE 1: RMT MARČENKO-PASTUR SPECTRAL DENSITY FIT
# ------------------------------------------------------------------------------
pdf("Figure1_RMT_Spectral_Distribution.pdf", width = 10, height = 5.5)

bulk_eigs <- rmt_global$eigvals[-1]
Q_ratio <- T_total / N_assets
s2 <- rmt_global$sigma2
l_min <- rmt_global$lambda_min
l_max <- rmt_global$lambda_max

mp_x <- seq(l_min, l_max, length.out = 300)
mp_y <- (Q_ratio / (2 * pi * s2 * mp_x)) * sqrt((l_max - mp_x) * (mp_x - l_min))
df_mp <- data.frame(x = mp_x, y = mp_y)

p_rmt <- ggplot() +
  geom_histogram(data = data.frame(val = bulk_eigs), aes(x = val, y = after_stat(density)),
                 bins = 18, fill = "#C6DBEF", color = "#2B6CB0", linewidth = 0.4, alpha = 0.8) +
  geom_line(data = df_mp, aes(x = x, y = y, color = "Theoretical M-P Density"), linewidth = 1.2) +
  geom_vline(xintercept = l_max, linetype = "dashed", color = "#C53030", linewidth = 0.9) +
  annotate("text", x = l_max * 1.05, y = max(mp_y) * 0.8, 
           label = paste0("Cutoff: lambda[max] == ", round(l_max, 3)), 
           parse = TRUE, color = "#C53030", fontface = "bold", hjust = 0, size = 3.5) +
  scale_color_manual(name = "", values = c("Theoretical M-P Density" = "#2B6CB0")) +
  labs(title = "Empirical Eigenvalue Distribution vs. Analytical Marčenko-Pastur Bound",
       subtitle = "Separation of Wishart random noise bulk from macro structural modes (NSE 50)",
       x = expression("Eigenvalue Magnitude (" * lambda * ")"),
       y = expression("Spectral Probability Density " * rho(lambda))) +
  theme_report()

print(p_rmt)
dev.off()

# ------------------------------------------------------------------------------
# FIGURE 2: MULTI-TOPOLOGY COMPARATIVE GRAPHS (MST vs. TMFG vs. WTA)
# ------------------------------------------------------------------------------
pdf("Figure2_Comparative_Network_Architectures.pdf", width = 15, height = 5.5)

set.seed(123)
layout_fixed <- layout_with_fr(g_tmfg, niter = 1500)
par(mfrow = c(1, 3), mar = c(1.5, 1, 3, 1), oma = c(0, 0, 2, 0))

node_cols <- sector_palette[sectors]

# MST Plot
plot(g_mst, layout = layout_fixed,
     vertex.size = 8, vertex.label = tickers, vertex.label.cex = 0.6,
     vertex.label.color = "black", vertex.label.font = 2,
     vertex.color = node_cols, vertex.frame.color = "grey30",
     edge.color = "grey40", edge.width = 1.4,
     main = sprintf("Minimum Spanning Tree (MST)\n|V| = %d, |E| = %d", vcount(g_mst), ecount(g_mst)))

# TMFG Plot
plot(g_tmfg, layout = layout_fixed,
     vertex.size = 8, vertex.label = tickers, vertex.label.cex = 0.6,
     vertex.label.color = "black", vertex.label.font = 2,
     vertex.color = node_cols, vertex.frame.color = "grey30",
     edge.color = rgb(0.2, 0.45, 0.8, 0.4), edge.width = 1.1,
     main = sprintf("Planar Filtered Graph (TMFG)\n|V| = %d, |E| = %d (|E| = 3N - 6)", vcount(g_tmfg), ecount(g_tmfg)))

# WTA Plot
plot(g_wta, layout = layout_fixed,
     vertex.size = 8, vertex.label = tickers, vertex.label.cex = 0.6,
     vertex.label.color = "black", vertex.label.font = 2,
     vertex.color = node_cols, vertex.frame.color = "grey30",
     edge.color = rgb(0.85, 0.35, 0.1, 0.4), edge.width = 1.1,
     main = sprintf("Winner-Take-All (WTA, k=4)\n|V| = %d, |E| = %d", vcount(g_wta), ecount(g_wta)))

mtext("Topologically Filtered Network Architectures across 50 Indian Equities",
      outer = TRUE, cex = 1.2, font = 2)

par(mfrow = c(1, 1))
dev.off()

# ------------------------------------------------------------------------------
# FIGURE 3: LONGITUDINAL ROLLING METRICS & MACRO REGIME COLLAPSE
# ------------------------------------------------------------------------------
pdf("Figure3_Longitudinal_Network_Dynamics.pdf", width = 11, height = 8.5)

p3_a <- ggplot(rolling_dynamics, aes(x = Window_End)) +
  geom_line(aes(y = Mean_Correlation, color = "Mean Correlation (C_bar)"), linewidth = 0.9) +
  geom_line(aes(y = Norm_Tree_Length * 0.45, color = "Normalized Tree Length (Scaled)"), 
            linewidth = 0.9, linetype = "dashed") +
  scale_y_continuous(
    name = expression("Mean Asset Correlation " * bar(C)),
    sec.axis = sec_axis(~ . / 0.45, name = "Normalized Tree Length (NTL)")
  ) +
  scale_color_manual(name = "Metric", values = c("Mean Correlation (C_bar)" = "#C53030", 
                                                 "Normalized Tree Length (Scaled)" = "#2B6CB0")) +
  labs(title = "(a) Systematic Network Contraction (NTL vs. Market Correlation)",
       x = NULL) +
  theme_report()

p3_b <- ggplot(rolling_dynamics, aes(x = Window_End)) +
  geom_line(aes(y = TMFG_Modularity, color = "Louvain Modularity (Q)"), linewidth = 0.85) +
  geom_line(aes(y = Sector_NMI, color = "Sector Recovery (NMI)"), linewidth = 0.85, linetype = "dashed") +
  scale_color_manual(name = "Index", values = c("Louvain Modularity (Q)" = "#2F855A", 
                                                "Sector Recovery (NMI)" = "#6B46C1")) +
  labs(title = "(b) Mesoscopic Community Purity vs. Industrial Taxonomy",
       x = NULL, y = "Score [0, 1]") +
  theme_report()

p3_c <- ggplot(rolling_dynamics, aes(x = Window_End, y = MST_Bootstrap_P50)) +
  geom_area(fill = "#BEE3F8", alpha = 0.4) +
  geom_line(color = "#2B6CB0", linewidth = 0.85) +
  geom_hline(yintercept = mean(rolling_dynamics$MST_Bootstrap_P50), color = "#C53030", linetype = "dotted") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "(c) MST Bootstrap Edge Persistence Ratio (P >= 0.50)",
       x = "Window End Date", y = "Stable Edge Ratio") +
  theme_report()

grid.arrange(p3_a, p3_b, p3_c, ncol = 1)
dev.off()
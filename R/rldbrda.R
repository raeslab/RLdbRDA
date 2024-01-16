
#' Validate Input Data
#'
#' This internal function checks if the input data frame matches the corresponding metadata
#' in terms of row names and also verifies whether the metadata contains any missing
#' values.
#'
#' @param df A data frame containing the data to be validated.
#' @param meta A data frame containing metadata associated with the data. Both `df` and
#'     `meta` must have identical row names.
#'
#' @return None
validate_input <- function(df, meta) {
  # Validate input
  if (!all(rownames(df) == rownames(meta))) {
    stop("Data and metadata don't match, should have the same rownames.")
  }

  if (any(is.na(meta))) {
    stop("Missing values found in metadata. Metadata cannot contain missing values.")
  }
}

#' Calculate R-squared Values for a Single Feature
#'
#' This function computes the R-squared and adjusted R-squared for a single feature
#' from a distance matrix and its corresponding metadata using constrained
#' ordination (capscale) and permutation tests (anova.cca) from the 'vegan' package.
#'
#' @param distmat A distance matrix in which distances between items are stored.
#' @param meta A data frame containing metadata associated with the items in `distmat`.
#' @param feature A character string specifying the column name in `meta` for which the
#'     R-squared values are to be calculated.
#'
#' @return A data frame row with columns for the feature name, the F statistic 'Fa',
#'     the R-squared 'r2', the adjusted R-squared 'r2adj', the number of non-missing
#'     values 'N', and the associated p-value 'pval'.
#'
#' @importFrom vegan capscale
#' @importFrom vegan anova.cca
#' @importFrom vegan RsquareAdj
get_r2_single <- function(distmat, meta, feature) {
  capsc <- capscale(distmat ~ meta[, feature], na.action=na.omit)
  an <- anova.cca(capsc)

  Fa <- an["F"][[1]][[1]]
  r2 <- RsquareAdj(capsc)[[1]]
  r2adj <- RsquareAdj(capsc)[[2]]
  N <- nrow(na.exclude(meta[,feature,drop=FALSE]))
  pval <- an["Pr(>F)"][[1]][[1]]

  output <- cbind(feature, Fa,r2,r2adj,N,pval)

  return(output)
}

#' Calculate R-squared Values for All Features Across Distance Matrix
#'
#' This function applies `get_r2_single` to each feature in the metadata dataframe
#' and calculates the R-squared and adjusted R-squared using constrained ordination
#' and permutation tests for all features. It also adjusts p-values using the
#' Benjamini-Hochberg method.
#'
#' @param distmat A distance matrix containing distances between sampled entities.
#' @param meta A dataframe containing metadata associated with `distmat`.
#'
#' @return A dataframe with rows for each feature and columns for the feature name,
#'     the F statistic 'Fa', the R-squared 'r2', the adjusted R-squared 'r2adj',
#'     the number of non-missing values 'N', the raw p-values 'pval', and the
#'     Benjamini-Hochberg adjusted p-values 'padj'. Row names of the dataframe
#'     are set to the feature names.
#'
#' @importFrom stats p.adjust
get_r2 <- function(distmat, meta) {
  features = colnames(meta)

  all <- c()

  for (feature in features){
    es <- get_r2_single(distmat, meta, feature)
    all <- rbind(all, es)
  }

  all <- data.frame(all)
  all$padj <- p.adjust(all$pval,method="BH")

  rownames(all) <- all$feature

  return(all)
}

#' Perform Cumulative db-RDA Analysis
#'
#' This function performs forward stepwise distance-based redundancy analysis (db-RDA)
#' to determine the cumulative effects of metadata variables on a distance matrix.
#' It utilizes the `capscale` and `ordiR2step` functions from the vegan package.
#'
#' @param distmat A distance matrix where rows and columns correspond to samples.
#' @param meta A dataframe containing metadata for each sample in `distmat`.
#'
#' @return A dataframe with the R-squared values for each variable added sequentially
#'     into the ordination model. The dataframe has columns adjusted for better
#'     naming convention and includes a cumulative R-squared value column and
#'     step count for each variable, along with the number of non-missing samples
#'     'RDAcumul_N'.
#'
#' @importFrom vegan capscale
#' @importFrom vegan ordiR2step
get_cumul <- function(distmat, meta) {

  mod0=capscale(distmat ~ 1) #H0: unconstrained ordination
  mod1=capscale(distmat ~ ., data=meta) #H1: full constrained ordination, all metadata

  attach(meta)

  step.res<-ordiR2step(mod0, scope=formula(mod1), data=meta ,direction="forward", Pin = 1, R2scope = TRUE, pstep = 100, perm.max = 999, permutations=9999, trace = F) #forward stepwise dbRDA
  res=step.res$anova
  row.names(res)=gsub(pattern="\\+ ", "",row.names(res))
  colnames(res)=gsub(pattern="Pr\\(>F\\)", "pval", colnames(res)) # replace column name
  colnames(res)=paste0("RDAcumul_",colnames(res))
  res[,"RDAcumul_N"]=nrow(meta)

  detach(meta)

  return(res)
}

#' Combine Redundant and Non-redundant Effect Sizes
#'
#' This function takes the results from the `get_r2` function (representing redundant
#' effect sizes based on individual metadata variables) and the results from
#' the `get_cumul` function (representing non-redundant, cumulative effect sizes from
#' a stepwise db-RDA) and combines them into a single dataframe. The resulting dataframe
#' is ordered by redundant effect sizes (R-squared) in decreasing order and by
#' cumulative adjusted R-squared values.
#'
#' @param r2 A dataframe that is the output from `get_r2`, representing the
#'    R-squared values for individual features.
#' @param cumul A dataframe that is the output from `get_cumul`, representing
#'    the cumulative R-squared values from stepwise distance-based redundancy analysis.
#'
#' @return A dataframe containing the combined set of effect sizes, with rows
#'    representing variables and ordered by the R-squared values (`r2`) in
#'    descending order, and then by the adjusted cumulative R-squared values
#'    (`RDAcumul_R2.adj`). Column `Row.names` denotes the names of the variables.
combine_data <- function(r2, cumul){

  all=data.frame(merge(r2,cumul,by="row.names",all=T),row.names=1)
  all=all[order(all$r2,decreasing=TRUE),]
  all=all[order(all$RDAcumul_R2.adj),]

  return(all)
}

#' Distance-Based Redundancy Analysis Workflow
#'
#' This function executes a full db-RDA analysis on a given dataset and its corresponding
#' metadata. It first validates the input, calculates pairwise distances using the vegdist
#' function, determines the R-squared values for each metadata variable, filters
#' variables based on a significance cutoff, and then performs cumulative db-RDA on
#' significant variables. It returns a dataframe combining both the individual and
#' cumulative results.
#'
#' @param df A numeric matrix or dataframe with rows as samples and columns as
#'   species/variables for which the distance matrix will be computed. Assumed to be
#'   the community data.
#' @param meta A dataframe where rows correspond to the samples in 'df' and columns
#'   contain metadata variables.
#' @param method Optional; a character string indicating the distance measure to be used
#'   (e.g., "bray" for Bray-Curtis dissimilarity). Default is "bray".
#' @param p_cutoff Optional; a significance level for determining which metadata
#'   variables to include in subsequent analysis. Default is 0.05.
#'
#' @return A dataframe with combined R-squared and cumulative R-squared values for
#'   significant metadata variables. The dataframe is ordered by R-squared values
#'   (decreasing) and adjusted cumulative R-squared values.
#'
#' @importFrom vegan vegdist
#' @export
#' @examples
#' # To perform a db-RDA analysis with example data 'df' and 'meta':
#' res <- rldbrda(df, meta, method = "bray", p_cutoff = 0.05)
#'
rldbrda <- function(df, meta, method="bray", p_cutoff=0.05) {

  validate_input(df, meta)

  distmat=vegdist(df,method=method)

  r2 = get_r2(distmat, meta)

  sign_r2 = rownames(r2[which(r2$padj < p_cutoff),]) # selects only variables significant in step 1

  if (length(sign_r2) < 1) {
    stop("No significant features found, cannot continue!")
  }

  cumul <- get_cumul(distmat, meta[, sign_r2])

  out <- combine_data(r2, cumul)

  return(out)
}

#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @export
prepare_plot_data <- function(dbrda_data) {
  # Define the order of the y-axis
  plot_order = rev(row.names(dbrda_data))

  # Impute missing values in RDAcumul_R2.adj with the largest non-NA value
  insignificant_features <- row.names(dbrda_data)[is.na(dbrda_data$RDAcumul_R2.adj)]
  max_val <- max(dbrda_data$RDAcumul_R2.adj, na.rm = TRUE)
  dbrda_data$RDAcumul_R2.adj[is.na(dbrda_data$RDAcumul_R2.adj)] <- max_val

  # Filter the data and reshape it to long form
  plot_data = dbrda_data %>%
    filter(row.names(dbrda_data) != "<All variables>") %>%
    select(r2adj, RDAcumul_R2.adj) %>%
    mutate(r2adj = as.numeric(r2adj)) %>%
    rownames_to_column(var = "rowname") %>%
    pivot_longer(-rowname, names_to = "variable", values_to = "value") %>%
    mutate(significant = ifelse(rowname %in% insignificant_features, 0, 1))

  plot_data$rowname <- factor(plot_data$rowname, levels = plot_order)

  return(plot_data)
}

#' @import ggplot2
#' @export
plot_dbrda <- function(plot_data) {
  # Plot the data as a horizontal bar plot
  g <- ggplot(data = plot_data, aes(x = value, y = rowname, fill = variable)) +
    geom_bar(aes(alpha = significant), stat = 'identity', position = position_dodge2()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    xlab('Effect size') +
    ylab('Feature') +
    labs(fill = "") + # Hide title of the legend
    guides(alpha="none") + # Hide alpha legend
    scale_fill_manual(values=c("#60A68B", "#1F4068"),
                      labels=c(bquote(R^2), bquote('Cumulative' ~ R^2))) +
    scale_alpha_continuous(range = c(0.5, 1))

  return(g)
}

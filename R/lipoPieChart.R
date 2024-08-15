#' lipoPieChart
#' Pie charts and tables for lipoproteins.
#' 
#' @param data A data frame of the lipoproteins.
#' @param group A character vector of the groups the same length as the data, 
#' for example Female_HighBMI, Female_LowBMI, Male_HighBMI, Male_LowBMI.
#' @param cohort A character vector of the cohorts the same length as the data. 
#' Used to facet rows of the pie charts.
#' @param subfractions Logical for if subfractions of HDL, VLDL and LDL should
#' be produced. Default is TRUE. looks at subfractions considering the lipid 
#' distribution, so what portion do each of V1TG, V2TG, V3TG, V4TG, V5TG 
#' contribute to VLTG
#' @param subcompositions Logical for if subcompositions of HDL, VLDL and LDL 
#' should be produced. Looks at subfractions considering the lipid composition, 
#' so what portion does each of the raw V1TG, V1CH, V1PL contribute to the 
#' calculated V1TL.
#' @import fusion
#' @import stats
#' @import ggplot2
#' @import reshape2
#' @import ggrepel
#' @export

lipoPieChart <- function(data, group, subfractions = T, subcompositions = T, cohort = 1, optns = list()){

  #use fusion, which uses nmr.parser, to extend the lipo data

 df <- extendLipo(data = data) 
 
 #####plots######
 plotCombos <- list(`Lipoprotein Composition` = c("HDTL_calc", "IDTL_calc", "LDTL_calc", "VLTL_calc"),
                    `Particle Numbers` = c("VLPN_pct","IDPN_pct","L1PN_frac","L2PN_frac","L3PN_frac","L4PN_frac","L5PN_frac","L6PN_frac"),
                    `HDL Distribution` = c("HDTG_pct", "HDFC_pct", "HDCE_pct", "HDPL_pct"),
                    `LDL Distribution` = c("LDCE_pct" , "LDFC_pct" , "LDPL_pct", "LDTG_pct" ),
                    `IDL Distribution` = c("IDCE_pct" , "IDFC_pct", "IDPL_pct" , "IDTG_pct"),
                    `VLDL Distribution` = c("VLCE_pct", "VLFC_pct", "VLPL_pct", "VLTG_pct"))
 
 newnames <- list(
   `Lipoprotein Composition` = c("HDTL_calc" = "HDL", "IDTL_calc" = "IDL", "LDTL_calc" = "LDL", "VLTL_calc" = "VLDL"),
   `Particle Numbers` = c("VLPN_pct" = "VLDL", "IDPN_pct" = "IDL", "L1PN_frac" = "LDL1", "L2PN_frac" = "LDL2", "L3PN_frac" = "LDL3", "L4PN_frac" = "LDL4", "L5PN_frac" = "LDL5", "L6PN_frac" = "LDL6"),
   `HDL Distribution` = c("HDTG_pct" = "TG", "HDFC_pct" = "FC", "HDCE_pct" = "CE", "HDPL" = "PL"),
   `LDL Distribution` = c("LDCE_pct" = "CE", "LDFC_pct" = "FC", "LDPL_pct" = "PL", "LDTG_pct" = "TG"),
   `IDL Distribution` = c("IDCE_pct" = "CE", "IDFC_pct" = "FC", "IDPL_pct" ="PL", "IDTG_pct" ="TG"),
   `VLDL Distribution` = c("VLCE_pct" = "CE", "VLFC_pct" = "FC", "VLPL_pct" = "PL", "VLTG_pct" = "TG")
 )
 
 if(subfractions == T){
   plotCombos <- list(`main composition` = plotCombos,
                      `LDL Subfraction` = list(`TG` = c("L1TG_frac", "L2TG_frac", "L3TG_frac", "L4TG_frac", "L5TG_frac", "L6TG_frac"),
                                               `CH` = c("L1CH_frac", "L2CH_frac", "L3CH_frac", "L4CH_frac", "L5CH_frac", "L6CH_frac"),
                                               `FC` = c("L1FC_frac", "L2FC_frac", "L3FC_frac", "L4FC_frac", "L5FC_frac", "L6FC_frac"), 
                                               `CE` = c("L1CE_frac", "L2CE_frac", "L3CE_frac", "L4CE_frac", "L5CE_frac", "L6CE_frac"),
                                               `PL` = c("L1PL_frac", "L2PL_frac", "L3PL_frac", "L4PL_frac", "L5PL_frac", "L6PL_frac")),
                      `HDL Subfraction` = list(`TG` = c("H1TG_frac", "H2TG_frac", "H3TG_frac", "H4TG_frac"), 
                                               `CH` = c("H1CH_frac", "H2CH_frac", "H3CH_frac", "H4CH_frac"), 
                                               `FC` = c("H1FC_frac", "H2FC_frac", "H3FC_frac", "H4FC_frac"), 
                                               `CE` = c("H1CE_frac", "H2CE_frac", "H3CE_frac", "H4CE_frac"), 
                                               `PL` = c("H1PL_frac", "H2PL_frac", "H3PL_frac", "H4PL_frac")), 
                      `VLDL Subfraction` = list(`TG` = c("V1TG_frac", "V2TG_frac", "V3TG_frac", "V4TG_frac", "V5TG_frac"), 
                                                `CH` = c("V1CH_frac", "V2CH_frac", "V3CH_frac", "V4CH_frac", "V5CH_frac"), 
                                                `FC` = c("V1FC_frac", "V2FC_frac", "V3FC_frac", "V4FC_frac", "V5FC_frac"), 
                                                `CE` = c("V1CE_frac", "V2CE_frac", "V3CE_frac", "V4CE_frac", "V5CE_frac"), 
                                                `PL` = c("V1PL_frac", "V2PL_frac", "V3PL_frac", "V4PL_frac", "V5PL_frac")))
 }else{plotCombos <- list(`main composition` = plotCombos)}
 
 if(subcompositions == T){
   `LDL Subcomposition` = list( `L1` = c("L1TG_pct", "L1FC_pct", "L1CE_pct", "L1PL_pct"),
                                `L2` = c("L2TG_pct", "L2FC_pct", "L2CE_pct", "L2PL_pct"),
                                `L3` = c("L3TG_pct", "L3FC_pct", "L3CE_pct", "L3PL_pct"),
                                `L4` = c("L4TG_pct", "L4FC_pct", "L4CE_pct", "L4PL_pct"),
                                `L5` = c("L5TG_pct", "L5FC_pct", "L5CE_pct", "L5PL_pct"),
                                `L6` = c("L6TG_pct", "L6FC_pct", "L6CE_pct", "L6PL_pct"))
   `HDL Subcomposition` = list( `H1` = c("H1TG_pct", "H1FC_pct", "H1CE_pct", "H1PL_pct"),
                                `H2` = c("H2TG_pct", "H2FC_pct", "H2CE_pct", "H2PL_pct"),
                                `H3` = c("H3TG_pct", "H3FC_pct", "H3CE_pct", "H3PL_pct"),
                                `H4` = c("H4TG_pct", "H4FC_pct", "H4CE_pct", "H4PL_pct"))
   `VLDL Subcomposition` = list(`V1` = c("V1TG_pct", "V1FC_pct", "V1CE_pct", "V1PL_pct"),
                                `V2` = c("V2TG_pct", "V2FC_pct", "V2CE_pct", "V2PL_pct"),
                                `V3` = c("V3TG_pct", "V3FC_pct", "V3CE_pct", "V3PL_pct"),
                                `V4` = c("V4TG_pct", "V4FC_pct", "V4CE_pct", "V4PL_pct"),
                                `V5` = c("V5TG_pct", "V5FC_pct", "V5CE_pct", "V5PL_pct"))
   
   plotCombos <- c(plotCombos, 
                   list("LDL Subcomposition" = `LDL Subcomposition`, 
                        "HDL Subcomposition" = `HDL Subcomposition`,
                        "VLDL Subcomposition" = `VLDL Subcomposition`))
   
   
 }else{plotCombos <- plotCombos}
 
 
 df$group <- as.factor(group)
 df$cohort <- as.factor(cohort)
 
 unique_factors <- unique(df$group)
 
 df[is.na(df)] <- 0
 df[sapply(df, is.infinite)] <- 0
 
 #plotting
 piePlot <- list()
 for(j in names(plotCombos)){
   for(i in names(plotCombos[[j]])){
     lipoproteins <- c(plotCombos[[j]][[i]], "group", "cohort")
     
     medians <- aggregate(. ~ group + cohort, data = df[,lipoproteins], function(x) median(x, na.rm = TRUE))
     medians$total <- rowSums(medians[, -which(names(medians) %in% c("group", "cohort"))])
     
     for (col in names(medians)[-which(names(medians) %in% c("group", "cohort", "total"))]) {
       medians[[col]] <- round((medians[[col]] / medians$total) * 100, 2)
     }
     
     medians <- medians[,-which(names(medians) %in% c("total"))]
     
     
     #rename for graph purposes
     # category <- "Lipoprotein Composition"
     if(j == 'main composition'){
       for (category in names(newnames)) {
         for (oldname in names(newnames[[category]])) {
           newname <- newnames[[category]][oldname]
           if (oldname %in% colnames(medians)) {
             colnames(medians)[colnames(medians) == oldname] <- newname
           }
         }
       }
     } 
     if(grepl("Subfraction", x = j)){
       #remove suffix
       end <- c("TG_frac", "CH_frac", "FC_frac", "CE_frac", "PL_frac")
       for(k in end){
         colnames(medians) <- gsub(pattern = k, replacement = "", x = colnames(medians))
         
         #rename prefix
         replacement <- list("H" = "HDL", "L" = "LDL", "V" = "VLDL")
         for (pattern in names(replacement)) {
           colnames(medians) <- gsub(pattern = paste0("(^|[^A-Za-z])", pattern, "(?![A-Za-z])"), replacement = replacement[[pattern]], x = colnames(medians), perl = TRUE)
         }
       }
     }
     
     if(grepl("Subcomposition", x = j)){
       #remove suffix
       end <- c("_pct")
       colnames(medians) <- gsub(pattern = end, replacement = "", x = colnames(medians))
       
       colnames(medians) <- sapply(colnames(medians), function(name) {
         if (name != "group" && name != "cohort") {
           substr(name, 3, nchar(name))  # Remove the first two characters
         } else {
           name  # Leave "group" and "cohort" unchanged
         }
       })
     }
     
     #reshape data
     long_data <- reshape(
       medians,
       varying = list(names(medians)[-which(names(medians) %in% c("group", "cohort"))]),
       v.names = "y",
       idvar = c("group", "cohort"),
       times = names(medians)[-which(names(medians) %in% c("group", "cohort"))],
       timevar = "Subfraction",
       direction = "long"
     )
     
     long_data <- long_data[order(long_data$group, long_data$cohort, long_data$Subfraction, decreasing = TRUE), ]
     
     # Calculate proportions and positions
     long_data$prop <- ave(long_data$y, long_data$group, long_data$cohort, FUN = function(x) x / sum(x))
     long_data$ypos <- ave(long_data$prop, long_data$group, long_data$cohort, FUN = function(x) cumsum(x) - 0.5 * x)
     long_data$label <- scales::percent(long_data$prop, accuracy = 0.1)
     
     if(j == 'main composition'){
       title <- paste0(i)
     }else{
       title <- paste0(j, " (", i, ")")
     }
     
     piePlot[[j]][[i]] <-
       ggplot(data = long_data,
              aes(x = "",
                  y = prop,
                  fill = Subfraction)) +
       geom_bar(width = 1,
                stat = "identity",
                color = "white",
                alpha = 0.8) +
       coord_polar("y",
                   start = 0) +
       theme_void() +
       facet_grid(rows = vars(cohort),
                  cols = vars(group)) +
       labs(title = title) +
       # geom_text(aes(y = ypos, label = label), size = 3, color = "black") +
       geom_text_repel(
         aes(y = ypos,
             label = label),
         nudge_x = 0.5,
         show.legend = FALSE,
         segment.size = 0.2,
         segment.color = 'grey50'
       )
     # scale_fill_brewer(palette = "Set1")
     
   }
 }
 
 #####Tables#########
 df$group <- as.factor(as.numeric(as.factor(group)))
 unique_factors <- unique(df$group)
 
 df[is.na(df)] <- 0
 df[sapply(df, is.infinite)] <- 0
 
 
 ########tables
 tableCombos = list(`Lipoprotein Composition` = c("HDCE_calc", "HDTL_calc", "IDCE_calc", "IDTL_calc", "LDCE_calc", "LDTL_calc", "TBPN_calc", "VLCE_calc", "VLTL_calc"),
                    `Particle Numbers` = c("TBPN_calc", "VLPN_pct","IDPN_pct","L1PN_frac","L2PN_frac","L3PN_frac","L4PN_frac","L5PN_frac","L6PN_frac"),
                    `HDL Distribution` = c("HDTG_pct", "HDCH_pct", "HDFC_pct", "HDCE_pct", "HDPL_pct"),
                    `LDL Distribution` = c("LDCE_pct", "LDCH_pct", "LDFC_pct", "LDPL_pct", "LDTG_pct"),
                    `IDL Distribution` = c("IDCE_pct", "IDCH_pct", "IDFC_pct", "IDPL_pct", "IDPN_pct", "IDTG_pct"),
                    `VLDL Distribution` = c("VLCE_pct", "VLCH_pct", "VLFC_pct", "VLPL_pct", "VLPN_pct", "VLTG_pct")
 )
 
 #subfractions
 if(subfractions == T){
   tableCombos = list(`main composition` = tableCombos,
                      `LDL Subfraction` = list(`TG` = c("L1TG_frac", "L2TG_frac", "L3TG_frac", "L4TG_frac", "L5TG_frac", "L6TG_frac", "LDTG_pct"),
                                               `CH` = c("L1CH_frac", "L2CH_frac", "L3CH_frac", "L4CH_frac", "L5CH_frac", "L6CH_frac", "LDCH_pct"),
                                               `FC` = c("L1FC_frac", "L2FC_frac", "L3FC_frac", "L4FC_frac", "L5FC_frac", "L6FC_frac", "LDFC_pct"), 
                                               `CE` = c("L1CE_frac", "L2CE_frac", "L3CE_frac", "L4CE_frac", "L5CE_frac", "L6CE_frac", "LDCE_pct"),
                                               `PL` = c("L1PL_frac", "L2PL_frac", "L3PL_frac", "L4PL_frac", "L5PL_frac", "L6PL_frac", "LDPL_pct")),
                      `HDL Subfraction` = list(`TG` = c("H1TG_frac", "H2TG_frac", "H3TG_frac", "H4TG_frac", "HDTG_pct"), 
                                               `CH` = c("H1CH_frac", "H2CH_frac", "H3CH_frac", "H4CH_frac", "HDCH_pct"), 
                                               `FC` = c("H1FC_frac", "H2FC_frac", "H3FC_frac", "H4FC_frac", "HDFC_pct"), 
                                               `CE` = c("H1CE_frac", "H2CE_frac", "H3CE_frac", "H4CE_frac", "HDCE_pct"), 
                                               `PL` = c("H1PL_frac", "H2PL_frac", "H3PL_frac", "H4PL_frac", "HDPL_pct")), 
                      `VLDL Subfraction` = list(`TG` = c("V1TG_frac", "V2TG_frac", "V3TG_frac", "V4TG_frac", "V5TG_frac", "VLTG_pct"), 
                                                `CH` = c("V1CH_frac", "V2CH_frac", "V3CH_frac", "V4CH_frac", "V5CH_frac", "VLCH_pct"), 
                                                `FC` = c("V1FC_frac", "V2FC_frac", "V3FC_frac", "V4FC_frac", "V5FC_frac", "VLFC_pct"), 
                                                `CE` = c("V1CE_frac", "V2CE_frac", "V3CE_frac", "V4CE_frac", "V5CE_frac", "VLCE_pct"), 
                                                `PL` = c("V1PL_frac", "V2PL_frac", "V3PL_frac", "V4PL_frac", "V5PL_frac", "VLPL_pct"))
   )
   
 } else {tableCombos = list(`main composition` = tableCombos)}
 
 if(subcompositions == T){
   
   `LDL Subcomposition` = list( `L1` = c("L1TG_pct", "L1FC_pct", "L1CE_pct", "L1PL_pct", "L1TL_calc"),
                                `L2` = c("L2TG_pct", "L2FC_pct", "L2CE_pct", "L2PL_pct", "L2TL_calc"),
                                `L3` = c("L3TG_pct", "L3FC_pct", "L3CE_pct", "L3PL_pct", "L3TL_calc"),
                                `L4` = c("L4TG_pct", "L4FC_pct", "L4CE_pct", "L4PL_pct", "L4TL_calc"),
                                `L5` = c("L5TG_pct", "L5FC_pct", "L5CE_pct", "L5PL_pct", "L4TL_calc"),
                                `L6` = c("L6TG_pct", "L6FC_pct", "L6CE_pct", "L6PL_pct", "L6TL_calc"))
   `HDL Subcomposition` = list( `H1` = c("H1TG_pct", "H1FC_pct", "H1CE_pct", "H1PL_pct", "H1TL_calc"),
                                `H2` = c("H2TG_pct", "H2FC_pct", "H2CE_pct", "H2PL_pct", "H2TL_calc"),
                                `H3` = c("H3TG_pct", "H3FC_pct", "H3CE_pct", "H3PL_pct", "H3TL_calc"),
                                `H4` = c("H4TG_pct", "H4FC_pct", "H4CE_pct", "H4PL_pct", "H4TL_calc"))
   `VLDL Subcomposition` = list(`V1` = c("V1TG_pct", "V1FC_pct", "V1CE_pct", "V1PL_pct", "V1TL_calc"),
                                `V2` = c("V2TG_pct", "V2FC_pct", "V2CE_pct", "V2PL_pct", "V2TL_calc"),
                                `V3` = c("V3TG_pct", "V3FC_pct", "V3CE_pct", "V3PL_pct", "V3TL_calc"),
                                `V4` = c("V4TG_pct", "V4FC_pct", "V4CE_pct", "V4PL_pct", "V4TL_calc"),
                                `V5` = c("V5TG_pct", "V5FC_pct", "V5CE_pct", "V5PL_pct", "V5TL_calc"))
   
   tableCombos <- c(tableCombos, 
                    list("LDL Subcomposition" = `LDL Subcomposition`, 
                         "HDL Subcomposition" = `HDL Subcomposition`,
                         "VLDL Subcomposition" = `VLDL Subcomposition`))
   
 }else{tableCombos <- tableCombos}
 
 #set up 
 if(length(unique_factors) == 2){
   col <- 3
 } else{
   col <- 2 + (sum(1:length(unique_factors))-length(unique_factors))
 }
 
 combinations <- combn(unique_factors, 2)
 
 # Create column names from the combinations with the larger letter first
 column_names <- apply(combinations, 2, function(x) paste(sort(x, decreasing = TRUE), collapse = "-"))
 
 #create a list for the table and fill them with all combinations of groups for significance
 tabs <- list()
 for(j in names(tableCombos)){
   for(i in names(tableCombos[[j]])){
     lipoproteins <- tableCombos[[j]][[i]]
     
     sig <- as.data.frame(matrix(NA, nrow = length(lipoproteins), ncol = col))
     
     colnames(sig) <- c("lipoproteins", "pval", column_names)
     sig$lipoproteins <- lipoproteins
     
     for(lipo in lipoproteins){
       anova_model <- aov(as.formula(paste(lipo, "~ group")), data = df)
       tukey_results <- TukeyHSD(anova_model)
       # Extract p-values
       idx <- which(sig$lipoproteins == lipo)
       sig[idx, "pval"] <- as.numeric(summary(anova_model)[[1]]$`Pr(>F)`[1])
       tukey_pvalues <- tukey_results$group
       
       #put the tukey results in the correct columns
       for (name in rownames(tukey_pvalues)) {
         if (name %in% colnames(sig)) {
           sig[idx, name] <- 
             # tukey_pvalues[name, 4]
             
             if(tukey_pvalues[name, 4] < 0.05 & tukey_pvalues[name, 4] >= 0.01){"*"
             }else if (tukey_pvalues[name, 4] < 0.01 & tukey_pvalues[name, 4] >= 0.001) {
               "**"
             } else if (tukey_pvalues[name, 4] < 0.001) {
               "***"
             } else {
               ""
             }
         }
       }
     }
     
     # Create a mapping between numbers and words, rename columns appropriately 
     mapping <- setNames(unique(group), unique(df$group))
     
     testnames<- as.data.frame(mapping, check.names = F)
     testnames$rowName <- rownames(testnames)
     
     new_colnames <- colnames(sig)
     for (num in names(mapping)) {
       new_colnames <- gsub(num, mapping[num], new_colnames)
     }
     colnames(sig) <- new_colnames
     
     # Calculate mean and sd for each 'group' group and each 'calc' column
     agg_mean <- melt(aggregate(. ~ group, data = df, FUN = mean), id.vars = "group")
     agg_sd <- melt(aggregate(. ~ group, data = df, FUN = sd), id.vars = "group")
     
     new <- merge(x = agg_mean, y = agg_sd, by = c("group", "variable"))
     new$`mean±sd` <- paste0(round(new$value.x, 2), " ± ", round(new$value.y, 2))
     new <- new[,which(!(colnames(new) %in% c("value.x", "value.y")))]
     
     result <- dcast(new, variable ~ group, value.var = "mean±sd")
     #rename using mapping
     new_colnames <- colnames(result)
     for (num in names(mapping)) {
       new_colnames <- gsub(num, mapping[num], new_colnames)
     }
     colnames(result) <- new_colnames
     
     tble <- merge(result, sig, by.x = "variable", by.y = "lipoproteins")
     
     #remove the _pct where necessary
     # tble$variable <- gsub(pattern = "_pct", replacement = "", x = tble$variable)
     tble$variable <- gsub(pattern = "_(pct|calc|frac)", replacement = "", x = tble$variable)
     tabs[[j]][[i]] <- tble
     
   }
   
 }
 
 results<- list("tables" = tabs,
                "pieCharts" = piePlot)
 
 return(results)
 
}

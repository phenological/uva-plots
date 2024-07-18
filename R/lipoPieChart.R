#' lipoPieChart
#' Pie charts and tables for lipoproteins.
#' 
#' @param data A data frame of the lipoproteins.
#' @param group A character vector of the groups the same length as the data, 
#' for example Female_HighBMI, Female_LowBMI, Male_HighBMI, Male_LowBMI.
#' @param cohort A character vector of the cohorts the same length as the data. 
#' Used to facet rows of the pie charts.
#' @param subfractions Logical for if subfractions of HDL, VLDL and LDL should
#' be produced. Default is TRUE. 
#' @import stats
#' @import ggplot2
#' @import reshape2
#' @import ggrepel
#' @export

lipoPieChart <- function(data, group, subfractions = T, cohort = 1, optns = list()){

 #create data frame
  if(is(data)[1] == "dataElement"){
    data <- as.data.frame(apply(data@.Data,2,as.numeric))
  }
  if(is(data)[1] != "data.frame"){
    data <- as.data.frame(data)
  }

  df <- data

  lipo_name<-nmr.parser::getLipoTable()$ID

   #Are the correct lipoproteins present
  if(length(which(lipo_name %in% colnames(df))) != 112){
    stop("some lipoprotein parameters does not match")
  }
  
  #is there multiple cohorts
  # if(length(cohort) == 1){
  #   cohort <- rep_len(x = 1, length.out = nrow(data))
  # }

 #create the required column (all in the names already)
 calc <- data.frame(
   HDTL = df$HDTG + df$HDCH + df$HDPL,
   HDCE = df$HDCH - df$HDFC,
   VLTL = df$VLTG + df$VLCH + df$VLPL,
   VLCE = df$VLCH - df$VLFC,
   IDTL = df$IDTG + df$IDCH + df$IDPL,
   IDCE = df$IDCH - df$IDFC,
   LDTL = df$LDTG + df$LDCH + df$LDPL,
   LDCE = df$LDCH - df$LDFC,
   TBPN = df$VLPN + df$IDPN + df$L1PN + df$L2PN + df$L3PN + df$L4PN + df$L5PN + df$L6PN,
   HDA1 = df$H1A1 + df$H2A1 + df$H3A1 + df$H4A1,
   HDA2 = df$H1A2 + df$H2A2 + df$H3A2 + df$H4A2,
   LDAB = df$L1AB + df$L2AB + df$L3AB + df$L4AB + df$L5AB + df$L6AB
 )

 #create percentage df use calc and df
 #initialisation (6 rows)
 perc <- data.frame(
   HDCE = round(calc$HDCE / calc$HDTL, 4) * 100,
   VLCE = round(calc$VLCE / calc$VLTL, 4) * 100,
   IDCE = round(calc$IDCE / calc$IDTL, 4) * 100,
   LDCE = round(calc$LDCE / calc$LDTL, 4) * 100,
   VLPN = round(df$VLPN / calc$TBPN, 4) * 100,
   IDPN = round(df$IDPN / calc$TBPN, 4) * 100
 )

 letters <- c("H", "V", "L")
 ranges <- list(H = 1:4, V = 1:5, L = 1:6)
 prefixes <- c("HD", "VL", "ID", "LD")
 suffixes <- c("TG", "CH", "FC", "PL")
 suffixes2 <- c("A1", "A2")

 for (letter in letters) {
   for (i in ranges[[letter]]) {
     ch_col <- paste0(letter, i, "CH")
     fc_col <- paste0(letter, i, "FC")
     ce_col <- paste0(letter, i, "CE")
     calc_col <- if (letter == "V") "VLCE" else paste0(letter, "DCE")
     perc[[ce_col]] <- round((df[[ch_col]] - df[[fc_col]]) / calc[[calc_col]], 4) * 100
   }
 }

 for (letter in letters) {
   for (i in ranges[[letter]]) {
     for (suffix in suffixes) {
       for (prefix in prefixes){
         if (prefix == "ID") next  # Skip "ID"
         col <- paste0(letter, i, suffix)
         col2 <- paste0(prefix, suffix)
         perc[[col]] <- round(df[[col]] / df[[col2]], 4) * 100
       }
     }
   }
 }

 for (prefix in prefixes) {
   for (suffix in suffixes) {
     col <- paste0(prefix, suffix)
     calc_col <- paste0(prefix, "TL")
     perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100
   }
 }

   for (i in ranges[["H"]]) {
     for (suffix in suffixes2) {
       col <- paste0("H", i, suffix)
       calc_col <- paste0("HD", suffix)
       perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100
     }
   }

   for (i in ranges[["L"]]) {
   col <- paste0("L", i, "AB")
   calc_col <- paste0("LDAB")
   perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100

   col <- paste0("L", i, "PN")
   calc_col <- paste0("TBPN")
   perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100
 }

 
 colnames(perc)[colnames(perc) %in% c("HDCE", "VLCE", "IDCE", "LDCE")] <- paste0(colnames(perc)[colnames(perc) %in% c("HDCE", "VLCE", "IDCE", "LDCE")], "_perc")
 
 #####plots######
 plotCombos <- list(`Lipoprotein Composition` = c("HDTL", "IDTL", "LDTL", "VLTL"),
                    `Particle Numbers` = c("VLPN","IDPN","L1PN","L2PN","L3PN","L4PN","L5PN","L6PN"),
                    `HDL Distribution` = c("HDTG" , "HDFC", "HDCE_perc" , "HDPL" ),
                    `LDL Distribution` = c("LDCE_perc" , "LDFC" , "LDPL", "LDTG" ),
                    `IDL Distribution` = c("IDCE_perc" , "IDFC", "IDPL" , "IDTG"),
                    `VLDL Distribution` = c("VLCE_perc", "VLFC", "VLPL", "VLTG"))
 
 newnames <- list(
   `Lipoprotein Composition` = c("HDTL" = "HDL", "IDTL" = "IDL", "LDTL" = "LDL", "VLTL" = "VLDL"),
   `Particle Numbers` = c("VLPN" = "VLDL", "IDPN" = "IDL", "L1PN" = "LDL1", "L2PN" = "LDL2", "L3PN" = "LDL3", "L4PN" = "LDL4", "L5PN" = "LDL5", "L6PN" = "LDL6"),
   `HDL Distribution` = c("HDTG" = "TG", "HDFC" = "FC", "HDCE_perc" = "CE", "HDPL" = "PL"),
   `LDL Distribution` = c("LDCE_perc" = "CE", "LDFC" = "FC", "LDPL" = "PL", "LDTG" = "TG"),
   `IDL Distribution` = c("IDCE_perc" = "CE", "IDFC" = "FC", "IDPL" ="PL", "IDTG" ="TG"),
   `VLDL Distribution` = c("VLCE_perc" = "CE", "VLFC" = "FC", "VLPL" = "PL", "VLTG" = "TG")
 )
 
 if(subfractions == T){
   plotCombos <- list(`main` = plotCombos,
                      `LDL Subfraction` = list(`TG` = c("L1TG", "L2TG", "L3TG", "L4TG", "L5TG", "L6TG"),
                                               `CH` = c("L1CH", "L2CH", "L3CH", "L4CH", "L5CH", "L6CH"),
                                               `FC` = c("L1FC", "L2FC", "L3FC", "L4FC", "L5FC", "L6FC"), 
                                               `CE` = c("L1CE", "L2CE", "L3CE", "L4CE", "L5CE", "L6CE"),
                                               `PL` = c("L1PL", "L2PL", "L3PL", "L4PL", "L5PL", "L6PL")),
                      `HDL Subfraction` = list(`TG` = c("H1TG", "H2TG", "H3TG", "H4TG"), 
                                               `CH` = c("H1CH", "H2CH", "H3CH", "H4CH"), 
                                               `FC` = c("H1FC", "H2FC", "H3FC", "H4FC"), 
                                               `CE` = c("H1CE", "H2CE", "H3CE", "H4CE"), 
                                               `PL` = c("H1PL", "H2PL", "H3PL", "H4PL")), 
                      `VLDL Subfraction` = list(`TG` = c("V1TG", "V2TG", "V3TG", "V4TG", "V5TG"), 
                                                `CH` = c("V1CH", "V2CH", "V3CH", "V4CH", "V5CH"), 
                                                `FC` = c("V1FC", "V2FC", "V3FC", "V4FC", "V5FC"), 
                                                `CE` = c("V1CE", "V2CE", "V3CE", "V4CE", "V5CE"), 
                                                `PL` = c("V1PL", "V2PL", "V3PL", "V4PL", "V5PL")))
 }else{plotCombos <- list(`main` = plotCombos)}
 df$group <- (as.factor(group))
 df$cohort <- as.factor(cohort)
 
 unique_factors <- unique(df$group)
 
 all <- cbind(calc, perc)
 all$group <- df$group
 all$cohort <- df$cohort
 
 all[is.na(all)] <- 0
 all[sapply(all, is.infinite)] <- 0
 
 #plotting
 piePlot <- list()
 for(j in names(plotCombos)){
   for(i in names(plotCombos[[j]])){

     lipoproteins <- c(plotCombos[[j]][[i]], "group", "cohort")
     
     medians <- aggregate(. ~ group + cohort, data = all[,lipoproteins], function(x) median(x, na.rm = TRUE))
     medians$total <- rowSums(medians[, -which(names(medians) %in% c("group", "cohort"))])
     
     for (col in names(medians)[-which(names(medians) %in% c("group", "cohort", "total"))]) {
       medians[[col]] <- round((medians[[col]] / medians$total) * 100, 2)
     }
     
     medians <- medians[,-which(names(medians) %in% c("total"))]
     
     
     #rename for graph purposes
     # category <- "Lipoprotein Composition"
     if(j == 'main'){
       for (category in names(newnames)) {
         for (oldname in names(newnames[[category]])) {
           newname <- newnames[[category]][oldname]
           if (oldname %in% colnames(medians)) {
             colnames(medians)[colnames(medians) == oldname] <- newname
           }
         }
       }
     } else{
       #remove suffix
       end <- c("TG", "CH", "FC", "CE", "PL")
       for(k in end){
         colnames(medians) <- gsub(pattern = k, replacement = "", x = colnames(medians))
         
         #rename prefix
         replacement <- list("H" = "HDL", "L" = "LDL", "V" = "VLDL")
         for (pattern in names(replacement)) {
           colnames(medians) <- gsub(pattern = paste0("(^|[^A-Za-z])", pattern, "(?![A-Za-z])"), replacement = replacement[[pattern]], x = colnames(medians), perl = TRUE)
         }
       }
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
     
     if(j == 'main'){
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
 
 all <- cbind(calc, perc)
 all$group <- df$group
 
 all[is.na(all)] <- 0
 all[sapply(all, is.infinite)] <- 0
 
 
 ########tables
 tableCombos = list(`Lipoprotein Composition` = c("HDCE", "HDTL", "IDCE", "IDTL", "LDCE", "LDTL", "TBPN", "VLCE", "VLTL"),
                    `Particle Numbers` = c("TBPN", "VLPN","IDPN","L1PN","L2PN","L3PN","L4PN","L5PN","L6PN"),
                    `HDL Distribution` = c("HDTG", "HDCH", "HDFC", "HDCE_perc", "HDPL"),
                    `LDL Distribution` = c("LDCE_perc", "LDCH", "LDFC", "LDPL", "LDTG"),
                    `IDL Distribution` = c("IDCE_perc", "IDCH", "IDFC", "IDPL", "IDPN", "IDTG"),
                    `VLDL Distribution` = c("VLCE_perc", "VLCH", "VLFC", "VLPL", "VLPN", "VLTG")
 )
 
 #subfractions
 if(subfractions == T){
   tableCombos = list(`main` = tableCombos,
                      `LDL Subfraction` = list(`TG` = c("L1TG", "L2TG", "L3TG", "L4TG", "L5TG", "L6TG", "LDTG"),
                                               `CH` = c("L1CH", "L2CH", "L3CH", "L4CH", "L5CH", "L6CH", "LDCH"),
                                               `FC` = c("L1FC", "L2FC", "L3FC", "L4FC", "L5FC", "L6FC", "LDFC"), 
                                               `CE` = c("L1CE", "L2CE", "L3CE", "L4CE", "L5CE", "L6CE", "LDCE"),
                                               `PL` = c("L1PL", "L2PL", "L3PL", "L4PL", "L5PL", "L6PL", "LDPL")),
                      `HDL Subfraction` = list(`TG` = c("H1TG", "H2TG", "H3TG", "H4TG", "HDTG"), 
                                               `CH` = c("H1CH", "H2CH", "H3CH", "H4CH", "HDCH"), 
                                               `FC` = c("H1FC", "H2FC", "H3FC", "H4FC", "HDFC"), 
                                               `CE` = c("H1CE", "H2CE", "H3CE", "H4CE", "HDCE"), 
                                               `PL` = c("H1PL", "H2PL", "H3PL", "H4PL", "HDPL")), 
                      `VLDL Subfraction` = list(`TG` = c("V1TG", "V2TG", "V3TG", "V4TG", "V5TG", "VLTG"), 
                                                `CH` = c("V1CH", "V2CH", "V3CH", "V4CH", "V5CH", "VLCH"), 
                                                `FC` = c("V1FC", "V2FC", "V3FC", "V4FC", "V5FC", "VLFC"), 
                                                `CE` = c("V1CE", "V2CE", "V3CE", "V4CE", "V5CE", "VLCE"), 
                                                `PL` = c("V1PL", "V2PL", "V3PL", "V4PL", "V5PL", "VLPL")))
 } else {tableCombos = list(`main` = tableCombos)}
 
 
 
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
       anova_model <- aov(as.formula(paste(lipo, "~ group")), data = all)
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
     agg_mean <- melt(aggregate(. ~ group, data = all, FUN = mean), id.vars = "group")
     agg_sd <- melt(aggregate(. ~ group, data = all, FUN = sd), id.vars = "group")
     
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
     
     #remove the _perc where necessary
     tble$variable <- gsub(pattern = "_perc", replacement = "", x = tble$variable)
     
     tabs[[j]][[i]] <- tble
     
   }
   
 }
 
 results<- list("tables" = tabs,
                "pieCharts" = piePlot)
 
 return(results)
 
}


# #create the required column (all in the names already)
# calc <- data.frame(
#   HDTL = df$HDTG + df$HDCH + df$HDPL,
#   HDCE = df$HDCH - df$HDFC,
#   VLTL = df$VLTG + df$VLCH + df$VLPL,
#   VLCE = df$VLCH - df$VLFC,
#   IDTL = df$IDTG + df$IDCH + df$IDPL,
#   IDCE = df$IDCH - df$IDFC,
#   LDTL = df$LDTG + df$LDCH + df$LDPL,
#   LDCE = df$LDCH - df$LDFC,
#   TBPN = df$VLPN + df$IDPN + df$L1PN + df$L2PN + df$L3PN + df$L4PN + df$L5PN + df$L6PN,
#   HDA1 = df$H1A1 + df$H2A1 + df$H3A1 + df$H4A1,
#   HDA2 = df$H1A2 + df$H2A2 + df$H3A2 + df$H4A2,
#   LDAB = df$L1AB + df$L2AB + df$L3AB + df$L4AB + df$L5AB + df$L6AB
# )
#
# #create percentage df use calc and df
# #initialisation (6 rows)
# perc <- data.frame(
#   HDCE = round(calc$HDCE / calc$HDTL, 4) * 100,
#   VLCE = round(calc$VLCE / calc$VLTL, 4) * 100,
#   IDCE = round(calc$IDCE / calc$IDTL, 4) * 100,
#   LDCE = round(calc$LDCE / calc$LDTL, 4) * 100,
#   VLPN = round(df$VLPN / calc$TBPN, 4) * 100,
#   IDPN = round(df$IDPN / calc$TBPN, 4) * 100
# )


#XXCE (15 rows)
# perc$H1CE = round(c(df$H1CH - df$H1FC) / calc$HDCE, 4) * 100
# perc$H2CE = round(c(df$H2CH - df$H2FC) / calc$HDCE, 4) * 100
# perc$H3CE = round(c(df$H3CH - df$H3FC) / calc$HDCE, 4) * 100
# perc$H4CE = round(c(df$H4CH - df$H4FC) / calc$HDCE, 4) * 100
# perc$V1CE = round(c(df$V1CH - df$V1FC) / calc$VLCE, 4) * 100
# perc$V2CE = round(c(df$V2CH - df$V2FC) / calc$VLCE, 4) * 100
# perc$V3CE = round(c(df$V3CH - df$V3FC) / calc$VLCE, 4) * 100
# perc$V4CE = round(c(df$V4CH - df$V4FC) / calc$VLCE, 4) * 100
# perc$V5CE = round(c(df$V5CH - df$V5FC) / calc$VLCE, 4) * 100
# perc$L1CE = round(c(df$L1CH - df$L1FC) / calc$LDCE, 4) * 100
# perc$L2CE = round(c(df$L2CH - df$L2FC) / calc$LDCE, 4) * 100
# perc$L3CE = round(c(df$L3CH - df$L3FC) / calc$LDCE, 4) * 100
# perc$L4CE = round(c(df$L4CH - df$L4FC) / calc$LDCE, 4) * 100
# perc$L5CE = round(c(df$L5CH - df$L5FC) / calc$LDCE, 4) * 100
# perc$L6CE = round(c(df$L6CH - df$L6FC) / calc$LDCE, 4) * 100

# letters <- c("H", "V", "L")
# ranges <- list(H = 1:4, V = 1:5, L = 1:6)
#
# for (letter in letters) {
#   for (i in ranges[[letter]]) {
#     ch_col <- paste0(letter, i, "CH")
#     fc_col <- paste0(letter, i, "FC")
#     ce_col <- paste0(letter, i, "CE")
#     calc_col <- if (letter == "V") "VLCE" else paste0(letter, "DCE")
#     perc[[ce_col]] <- round((df[[ch_col]] - df[[fc_col]]) / calc[[calc_col]], 4) * 100
#   }
# }


#YX Y = HD, VL, ID, LD X = TG, CH, FC, PL (16 rows)
# perc$HDTG = round(df$HDTG / calc$HDTL, 4) * 100
# perc$HDCH = round(df$HDCH / calc$HDTL, 4) * 100
# perc$HDFC = round(df$HDFC / calc$HDTL, 4) * 100
# perc$HDPL = round(df$HDPL / calc$HDTL, 4) * 100
# perc$VLTG = round(df$VLTG / calc$VLTL, 4) * 100
# perc$VLCH = round(df$VLCH / calc$VLTL, 4) * 100
# perc$VLFC = round(df$VLFC / calc$VLTL, 4) * 100
# perc$VLPL = round(df$VLPL / calc$VLTL, 4) * 100
# perc$IDTG = round(df$IDTG / calc$IDTL, 4) * 100
# perc$IDCH = round(df$IDCH / calc$IDTL, 4) * 100
# perc$IDFC = round(df$IDFC / calc$IDTL, 4) * 100
# perc$IDPL = round(df$IDPL / calc$IDTL, 4) * 100
# perc$LDTG = round(df$LDTG / calc$LDTL, 4) * 100
# perc$LDCH = round(df$LDCH / calc$LDTL, 4) * 100
# perc$LDFC = round(df$LDFC / calc$LDTL, 4) * 100
# perc$LDPL = round(df$LDPL / calc$LDTL, 4) * 100

# prefixes <- c("HD", "VL", "ID", "LD")
# suffixes <- c("TG", "CH", "FC", "PL")
#
# for (prefix in prefixes) {
#   for (suffix in suffixes) {
#     col <- paste0(prefix, suffix)
#     calc_col <- paste0(prefix, "TL")
#     perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100
#   }
# }

#Z = A1 or A2 (for H) AB or PN (for L)
#x = 1:4 for H 1:6 for L
#perc$HxZ = round(df$HxZ / calc$HDZ, 4) * 100 (8 rows)
# perc$H1A1 = round(df$H1A1 / calc$HDA1, 4) * 100
# perc$H2A1 = round(df$H2A1 / calc$HDA1, 4) * 100
# perc$H3A1 = round(df$H3A1 / calc$HDA1, 4) * 100
# perc$H4A1 = round(df$H4A1 / calc$HDA1, 4) * 100
# perc$H1A2 = round(df$H1A2 / calc$HDA2, 4) * 100
# perc$H2A2 = round(df$H2A2 / calc$HDA2, 4) * 100
# perc$H3A2 = round(df$H3A2 / calc$HDA2, 4) * 100
# perc$H4A2 = round(df$H4A2 / calc$HDA2, 4) * 100


# ranges <- 1:4
# suffixes <- c("A1", "A2")
#
# for (i in ranges) {
#   for (suffix in suffixes) {
#     col <- paste0("H", i, suffix)
#     calc_col <- paste0("HD", suffix)
#     perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100
#   }
# }




#perc$LxY = round(df$LxY / calc$WAB, 4) * 100
# perc$L1AB = round(df$L1AB / calc$LDAB, 4) * 100
# perc$L2AB = round(df$L2AB / calc$LDAB, 4) * 100
# perc$L3AB = round(df$L3AB / calc$LDAB, 4) * 100
# perc$L4AB = round(df$L4AB / calc$LDAB, 4) * 100
# perc$L5AB = round(df$L5AB / calc$LDAB, 4) * 100
# perc$L6AB = round(df$L6AB / calc$LDAB, 4) * 100
#
# perc$L1PN = round(df$L1PN / calc$TBPN, 4) * 100
# perc$L2PN = round(df$L2PN / calc$TBPN, 4) * 100
# perc$L3PN = round(df$L3PN / calc$TBPN, 4) * 100
# perc$L4PN = round(df$L4PN / calc$TBPN, 4) * 100
# perc$L5PN = round(df$L5PN / calc$TBPN, 4) * 100
# perc$L6PN = round(df$L6PN / calc$TBPN, 4) * 100

# ranges <- 1:6
# for (i in ranges) {
#   col <- paste0("L", i, "AB")
#   calc_col <- paste0("LDAB")
#   perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100
# }
#
# for (i in ranges) {
#   col <- paste0("L", i, "PN")
#   calc_col <- paste0("TBPN")
#   perc[[col]] <- round(df[[col]] / calc[[calc_col]], 4) * 100
# }
#

### new below ###

# YxZW
#Y = H, V, L
#x = 1:4(for H) 1:5(for V) or 1:6(for L)
#Z = TG, CH, FC, PL
#W = D(for Y = H or L), L(for Y = V), (60 rows)
# perc$H1TG = round(df$H1TG / df$HDTG, 4) * 100
# perc$H2TG = round(df$H2TG / df$HDTG, 4) * 100
# perc$H3TG = round(df$H3TG / df$HDTG, 4) * 100
# perc$H4TG = round(df$H4TG / df$HDTG, 4) * 100
#
# perc$H1CH = round(df$H1CH / df$HDCH, 4) * 100
# perc$H2CH = round(df$H2CH / df$HDCH, 4) * 100
# perc$H3CH = round(df$H3CH / df$HDCH, 4) * 100
# perc$H4CH = round(df$H4CH / df$HDCH, 4) * 100
#
# perc$H1FC = round(df$H1FC / df$HDFC, 4) * 100
# perc$H2FC = round(df$H2FC / df$HDFC, 4) * 100
# perc$H3FC = round(df$H3FC / df$HDFC, 4) * 100
# perc$H4FC = round(df$H4FC / df$HDFC, 4) * 100
#
# perc$H1PL = round(df$H1PL / df$HDPL, 4) * 100
# perc$H2PL = round(df$H2PL / df$HDPL, 4) * 100
# perc$H3PL = round(df$H3PL / df$HDPL, 4) * 100
# perc$H4PL = round(df$H4PL / df$HDPL, 4) * 100
#
# perc$V1TG = round(df$V1TG / df$VLTG, 4) * 100
# perc$V2TG = round(df$V2TG / df$VLTG, 4) * 100
# perc$V3TG = round(df$V3TG / df$VLTG, 4) * 100
# perc$V4TG = round(df$V4TG / df$VLTG, 4) * 100
# perc$V5TG = round(df$V5TG / df$VLTG, 4) * 100
#
# perc$V1CH = round(df$V1CH / df$VLCH, 4) * 100
# perc$V2CH = round(df$V2CH / df$VLCH, 4) * 100
# perc$V3CH = round(df$V3CH / df$VLCH, 4) * 100
# perc$V4CH = round(df$V4CH / df$VLCH, 4) * 100
# perc$V5CH = round(df$V5CH / df$VLCH, 4) * 100
#
# perc$V1FC = round(df$V1FC / df$VLFC, 4) * 100
# perc$V2FC = round(df$V2FC / df$VLFC, 4) * 100
# perc$V3FC = round(df$V3FC / df$VLFC, 4) * 100
# perc$V4FC = round(df$V4FC / df$VLFC, 4) * 100
# perc$V5FC = round(df$V5FC / df$VLFC, 4) * 100
#
# perc$V1PL = round(df$V1PL / df$VLPL, 4) * 100
# perc$V2PL = round(df$V2PL / df$VLPL, 4) * 100
# perc$V3PL = round(df$V3PL / df$VLPL, 4) * 100
# perc$V4PL = round(df$V4PL / df$VLPL, 4) * 100
# perc$V5PL = round(df$V5PL / df$VLPL, 4) * 100
#
# perc$L1TG = round(df$L1TG / df$LDTG, 4) * 100
# perc$L2TG = round(df$L2TG / df$LDTG, 4) * 100
# perc$L3TG = round(df$L3TG / df$LDTG, 4) * 100
# perc$L4TG = round(df$L4TG / df$LDTG, 4) * 100
# perc$L5TG = round(df$L5TG / df$LDTG, 4) * 100
# perc$L6TG = round(df$L6TG / df$LDTG, 4) * 100
#
# perc$L1CH = round(df$L1CH / df$LDCH, 4) * 100
# perc$L2CH = round(df$L2CH / df$LDCH, 4) * 100
# perc$L3CH = round(df$L3CH / df$LDCH, 4) * 100
# perc$L4CH = round(df$L4CH / df$LDCH, 4) * 100
# perc$L5CH = round(df$L5CH / df$LDCH, 4) * 100
# perc$L6CH = round(df$L6CH / df$LDCH, 4) * 100
#
# perc$L1FC = round(df$L1FC / df$LDFC, 4) * 100
# perc$L2FC = round(df$L2FC / df$LDFC, 4) * 100
# perc$L3FC = round(df$L3FC / df$LDFC, 4) * 100
# perc$L4FC = round(df$L4FC / df$LDFC, 4) * 100
# perc$L5FC = round(df$L5FC / df$LDFC, 4) * 100
# perc$L6FC = round(df$L6FC / df$LDFC, 4) * 100
#
# perc$L1PL = round(df$L1PL / df$LDPL, 4) * 100
# perc$L2PL = round(df$L2PL / df$LDPL, 4) * 100
# perc$L3PL = round(df$L3PL / df$LDPL, 4) * 100
# perc$L4PL = round(df$L4PL / df$LDPL, 4) * 100
# perc$L5PL = round(df$L5PL / df$LDPL, 4) * 100
# perc$L6PL = round(df$L6PL / df$LDPL, 4) * 100

# Define the parameters for the loop
# letters <- c("H", "V", "L")
# ranges <- list(H = 1:4, V = 1:5, L = 1:6)
# suffixes <- c("TG", "CH", "FC", "PL")
# prefixes <- c("HD", "VL", "LD")
#
# for (letter in letters) {
#   for (i in ranges[[letter]]) {
#     for (suffix in suffixes) {
#       for(prefix in prefixes){
#         col <- paste0(letter, i, suffix)
#         col2 <- paste0(prefix, suffix)
#         perc[[col]] <- round(df[[col]] / df[[col2]], 4) * 100
#       }
#     }
#   }
# }


#packages needed
#reshape2




lipoData$group <- as.factor(as.numeric(as.factor(lipoData$category)))
unique_factors <- unique(lipoData$group)
# perc$group <- lipoData$group
# calc$group <- lipoData$group
# 
# perc[is.na(perc)] <- 0
# perc[sapply(perc, is.infinite)] <- 0
# 
# calc[is.na(calc)] <- 0
# calc[sapply(calc, is.infinite)] <- 0

all <- cbind(calc, perc)
all$group <- lipoData$group

all[is.na(all)] <- 0
all[sapply(all, is.infinite)] <- 0
########tables######
tableCombos = list(`Lipoprotein Composition` = c("HDCE", "HDTL", "IDCE", "IDTL", "LDCE", "LDTL", "TBPN", "VLCE", "VLTL"),
                   `Particle Numbers` = c("TBPN", "VLPN","IDPN","L1PN","L2PN","L3PN","L4PN","L5PN","L6PN"),
                   `HDL Distribution` = c("HDTG", "HDCH", "HDFC", "HDCE_perc", "HDPL"),
                   `LDL Distribution` = c("LDCE_perc", "LDCH", "LDFC", "LDPL", "LDTG"),
                   `IDL Distribution` = c("IDCE_perc", "IDCH", "IDFC", "IDPL", "IDPN", "IDTG"),
                   `VLDL Distribution` = c("VLCE_perc", "VLCH", "VLFC", "VLPL", "VLPN", "VLTG")
                   )

# cat(paste(shQuote((VLDL[["PL"]]), type = "cmd"), collapse=", "))
#all from perc

#subfractions
subfractions <- TRUE
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
trial <- list()
for(j in names(tableCombos)){
  for(i in names(tableCombos[[j]])){
    lipoproteins <- tableCombos[[j]][[i]]
    j = "main"
    i = "Lipoprotein Composition"
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
    mapping <- setNames(unique(lipoData$category), unique(lipoData$group))
    
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
    
    trial[[j]][[i]] <- tble
    
  }
  
}




# for(i in names(tableCombos)){
#   lipoproteins <- tableCombos[[i]]
#   
#   sig <- as.data.frame(matrix(NA, nrow = length(lipoproteins), ncol = col))
#   
#   colnames(sig) <- c("lipoproteins", "pval", column_names)
#   sig$lipoproteins <- lipoproteins
#   
#   for(lipo in lipoproteins){
#     anova_model <- aov(as.formula(paste(lipo, "~ group")), data = all)
#     tukey_results <- TukeyHSD(anova_model)
#     # Extract p-values
#     idx <- which(sig$lipoproteins == lipo)
#     sig[idx, "pval"] <- as.numeric(summary(anova_model)[[1]]$`Pr(>F)`[1])
#     tukey_pvalues <- tukey_results$group
#     
#     #put the tukey results in the correct columns
#     for (name in rownames(tukey_pvalues)) {
#       if (name %in% colnames(sig)) {
#         sig[idx, name] <- 
#           # tukey_pvalues[name, 4]
#           
#           if(tukey_pvalues[name, 4] < 0.05 & tukey_pvalues[name, 4] >= 0.01){"*"
#           }else if (tukey_pvalues[name, 4] < 0.01 & tukey_pvalues[name, 4] >= 0.001) {
#             "**"
#           } else if (tukey_pvalues[name, 4] < 0.001) {
#             "***"
#           } else {
#             ""
#           }
#       }
#     }
#   }
#   
#   # Create a mapping between numbers and words, rename columns appropriately 
#   mapping <- setNames(unique(lipoData$category), unique(lipoData$group))
#   
#   testnames<- as.data.frame(mapping, check.names = F)
#   testnames$rowName <- rownames(testnames)
#   
#   new_colnames <- colnames(sig)
#   for (num in names(mapping)) {
#     new_colnames <- gsub(num, mapping[num], new_colnames)
#   }
#   colnames(sig) <- new_colnames
#   
#   # Calculate mean and sd for each 'group' group and each 'calc' column
#   agg_mean <- melt(aggregate(. ~ group, data = all, FUN = mean), id.vars = "group")
#   agg_sd <- melt(aggregate(. ~ group, data = all, FUN = sd), id.vars = "group")
#   
#   new <- merge(x = agg_mean, y = agg_sd, by = c("group", "variable"))
#   new$`mean±sd` <- paste0(round(new$value.x, 2), " ± ", round(new$value.y, 2))
#   new <- new[,which(!(colnames(new) %in% c("value.x", "value.y")))]
#   
#   result <- dcast(new, variable ~ group, value.var = "mean±sd")
#   #rename using mapping
#   new_colnames <- colnames(result)
#   for (num in names(mapping)) {
#     new_colnames <- gsub(num, mapping[num], new_colnames)
#   }
#   colnames(result) <- new_colnames
#   
#   tble <- merge(result, sig, by.x = "variable", by.y = "lipoproteins")
#   
#   trial[[i]] <- tble
#   
# }


# #if only 2 factors, just need the second column, since first is all NA (control to control)
# if(length(unique_factors) == 2){
#   #cd_df <- as.data.frame(cd_df[,2])
#   cd_df <- cd_df[, 2, drop = FALSE]
# }
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

# lipoData$group <- as.factor((as.factor(lipoData$category)))
# unique_factors <- unique(lipoData$group)
# 
# all <- cbind(calc, perc)
# all$group <- lipoData$group
# 
# all[is.na(all)] <- 0
# all[sapply(all, is.infinite)] <- 0
# 
# #plotting
# piePlot <- list()
# for(j in names(plotCombos)){
#   for(i in names(plotCombos[[j]])){
#     lipoproteins <- c(plotCombos[[j]][[i]], "group")
# 
# medians <- aggregate(. ~ group, data = all[,lipoproteins], function(x) median(x, na.rm = TRUE))
# medians$total <- rowSums(medians[, -which(names(medians) == "group")])
# 
# for (col in names(medians)[-which(names(medians) %in% c("group", "total"))]) {
#   medians[[col]] <- round((medians[[col]] / medians$total) * 100, 2)
# }
# 
# medians <- medians[,-which(names(medians) %in% c("total"))]
# 
# 
#   #rename for graph purposes
#   # category <- "Lipoprotein Composition"
#   if(j == 'main'){
#     for (category in names(newnames)) {
#       for (oldname in names(newnames[[category]])) {
#         newname <- newnames[[category]][oldname]
#         if (oldname %in% colnames(medians)) {
#           colnames(medians)[colnames(medians) == oldname] <- newname
#         }
#       }
#     }
#   } else{
#     #remove suffix
#     end <- c("TG", "CH", "FC", "CE", "PL")
#     for(k in end){
#       colnames(medians) <- gsub(pattern = k, replacement = "", x = colnames(medians))
# 
#       #rename prefix
#       replacement <- list("H" = "HDL", "L" = "LDL", "V" = "VLDL")
#       for (pattern in names(replacement)) {
#         colnames(medians) <- gsub(pattern = paste0("(^|[^A-Za-z])", pattern, "(?![A-Za-z])"), replacement = replacement[[pattern]], x = colnames(medians), perl = TRUE)
#       }
#     }
#   }
#   
# #reshape data
# long_data <- reshape(
#   medians,
#   varying = list(names(medians)[-which(names(medians) %in% c("group"))]),
#   v.names = "y",
#   idvar = "group",
#   times = names(medians)[-which(names(medians) %in% c("group"))],
#   timevar = "Subfraction",
#   direction = "long"
# )
# 
# # Remove the row.names column added by reshape
# long_data <- long_data[order(long_data$group), ]
# long_data$labels <- paste0(long_data$y, " %")
# long_data$Subfraction <- as.factor(long_data$Subfraction)
# long_data$label_pos <- NA
# 
# long_data <-
#   do.call(what = rbind,
#           args = lapply(split(long_data, long_data$group), function(df) {
#   df$label_pos <- cumsum(df$y) - 0.5 * df$y
#   return(df)
# }))
# 
# if(j == 'main'){
#   title <- paste0(i)
# }else{
#   title <- paste0(j, " (", i, ")")
# }
#  
# 
# piePlot[[j]][[i]] <-
# ggplot(data = long_data, aes(x = "", y = y, fill = Subfraction)) +
#   geom_bar(stat = "identity", width = 1, color = "white") +
#   coord_polar("y", start = 0) +
#   theme_void() +
#   facet_wrap(~group)+
#   labs(title = title)
# 
# 
#   }
# }



####trial####

load("~/git/phenological/mva-plots/data/lipoData.rda")

lipoData$cohort<- rep_len(x = "Aus", length.out = nrow(lipoData))
lipoData2 <- lipoData
lipoData2$cohort<- rep_len(x = "USA", length.out = nrow(lipoData))
df <- rbind(lipoData, lipoData2)
group <- df$category
cohort <- df$cohort
df <- df[,1:112]
subfractions = F

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

df$group <- as.factor(group)
df$cohort <- as.factor(cohort)  # Assuming cohort is also a factor

unique_factors <- unique(df$group)

all <- cbind(calc, perc)
all$group <- df$group
all$cohort <- df$cohort  # Add cohort to all

all[is.na(all)] <- 0
all[sapply(all, is.infinite)] <- 0

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
#plotting
piePlot <- list()
for(j in names(plotCombos)){
  for(i in names(plotCombos[[j]])){
    # j <- "main"
    # i <- "Lipoprotein Composition"
    lipoproteins <- c(plotCombos[[j]][[i]], "group", "cohort")
    
    medians <- aggregate(. ~ group + cohort, data = all[,lipoproteins], function(x) median(x, na.rm = TRUE))
    medians$total <- rowSums(medians[, -which(names(medians) %in% c("group", "cohort"))])
    
    for (col in names(medians)[-which(names(medians) %in% c("group", "cohort", "total"))]) {
      medians[[col]] <- round((medians[[col]] / medians$total) * 100, 2)
    }
    
    medians <- medians[,-which(names(medians) %in% c("total"))]
    
    #rename for graph purposes
    if(j == 'main'){
      for (category in names(newnames)) {
        for (oldname in names(newnames[[category]])) {
          newname <- newnames[[category]][oldname]
          if (oldname %in% colnames(medians)) {
            colnames(medians)[colnames(medians) == oldname] <- newname
          }
        }
      }
    } else {
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
    
    # Plotting with ggplot2 and ggrepel
    ggplot(data = long_data, aes(x = "", y = prop, fill = Subfraction)) +
      geom_bar(width = 1, stat = "identity", color = "white", alpha = 0.8) +
      coord_polar("y", start = 0) +
      theme_void() +
      facet_grid(rows = vars(cohort), cols = vars(group)) +
      labs(title = title) +
      geom_text(aes(y = ypos, label = label), size = 3, color = "white") +
      scale_fill_brewer(palette = "Set1")
    
 ##########################   
    
    
    
    
    
    
    long_data <- long_data[order(long_data$group, long_data$cohort), ]
    long_data$labels <- paste0(long_data$y, " %")
    long_data$Subfraction <- as.factor(long_data$Subfraction)
    long_data$label_pos <- NA
    
    long_data <-
      do.call(what = rbind,
              args = lapply(split(long_data, list(long_data$group, long_data$cohort)), function(df) {
                df$label_pos <- cumsum(df$y) - (0.5 * df$y)
                return(df)
              }))
    
    
    if(j == 'main'){
      title <- paste0(i)
    }else{
      title <- paste0(j, " (", i, ")")
    }
    
    # piePlot[[j]][[i]] <-
   pp <-   ggplot(data = long_data, aes(x = "", y = y, fill = Subfraction)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      theme_void() +
      facet_grid(rows = vars(cohort), cols = vars(group)) +
      labs(title = title)  +
     geom_text(aes(y = label_pos, label = labels), size=3, color = "black")
      # geom_text_repel(
      #   aes(y = label_pos, label = labels),
      #   # nudge_x = 0.5,
      #   show.legend = FALSE,
      #   segment.size = 0.2,
      #   segment.color = 'grey50'
      # )
   
   # SubSegment<- c('S1','S2','S3','S4')
   # v <- c(100, 300, 500, 200)
   # df<- cbind.data.frame(SubSegment, v)
   # 
   # #calculations for % labels in chart
   # df <- df %>% 
   #   arrange(desc(SubSegment)) %>%
   #   mutate(prop = v / sum(df$v)) %>%
   #   mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
   #   mutate(label= prop*1)
   # df[5] = sapply(df[5], function(x) scales::percent(x, accuracy = 0.1))
   # 
   # plot.ex <- ggplot(df, aes(x = "", y = prop, fill = SubSegment)) +
   #   geom_bar(width = 1, stat = "identity", color="white", alpha=0.8) +
   #   coord_polar("y", start = 0) +
   #   theme_void() + 
   #   geom_text(aes(y = ypos, label = label), size=3, color = "white") +
   #   scale_fill_brewer(palette="Set1") 
  }
}

#####tables multicohort##########
df$group <- as.factor(as.numeric(as.factor(group)))
df$cohort <- as.factor(cohort)

unique_factors <- unique(df$group)

all <- cbind(calc, perc)
all$group <- df$group
all$cohort <- df$cohort

all[is.na(all)] <- 0
all[sapply(all, is.infinite)] <- 0


tableCombos = list(`Lipoprotein Composition` = c("HDCE", "HDTL", "IDCE", "IDTL", "LDCE", "LDTL", "TBPN", "VLCE", "VLTL"),
                   `Particle Numbers` = c("TBPN", "VLPN","IDPN","L1PN","L2PN","L3PN","L4PN","L5PN","L6PN"),
                   `HDL Distribution` = c("HDTG", "HDCH", "HDFC", "HDCE_perc", "HDPL"),
                   `LDL Distribution` = c("LDCE_perc", "LDCH", "LDFC", "LDPL", "LDTG"),
                   `IDL Distribution` = c("IDCE_perc", "IDCH", "IDFC", "IDPL", "IDPN", "IDTG"),
                   `VLDL Distribution` = c("VLCE_perc", "VLCH", "VLFC", "VLPL", "VLPN", "VLTG")
)

# cat(paste(shQuote((VLDL[["PL"]]), type = "cmd"), collapse=", "))
#all from perc

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
trial <- list()
for(j in names(tableCombos)){
  for(i in names(tableCombos[[j]])){
    lipoproteins <- tableCombos[[j]][[i]]
    j = "main"
    i = "Lipoprotein Composition"
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
    mapping <- setNames(unique(lipoData$category), unique(lipoData$group))
    
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
    
    trial[[j]][[i]] <- tble
    
  }
  
}



############etablished test##############
test_that("lipoprotein data frame made correctly ", {
  
  load("~/git/phenological/mva-plots/data/lipoData.rda")
  
  
  df <- lipoData[,1:112]
  
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
  
  #there should not be an duplicates of names from calc in perc. If there is, change them to have _perc.
  m <- colnames(calc) %in% colnames(perc)
  expect_false(object = TRUE %in% m)
  
  
})

test_that("table is made", {
  
  #set up 3 groups for lipoData instead of 2
  a_indices <- which(lipoData$category == "A")
  a_to_change <- sample(a_indices, ceiling(length(a_indices) / 3))
  lipoData$category[a_to_change] <- "C"
  
  # Change a third of "B" entries to "C"
  b_indices <- which(lipoData$category == "B")
  b_to_change <- sample(b_indices, ceiling(length(b_indices) / 3))
  lipoData$category[b_to_change] <- "C"
})


test_that("works with subfractions",{
  
  load("~/git/phenological/mva-plots/data/lipoData.rda")
  df <- lipoData[,1:112]
  test <- lipoPieChart(data = df, group = lipoData$category, subfractions = T)
  
  expect_equal(object = length(names(test[["tables"]])), expected = 4)
  
  #are the main plots there
  expect_contains(object = names(test[["pieCharts"]][["main"]]), expected = c("Lipoprotein Composition", "Particle Numbers", "HDL Distribution",   "LDL Distribution",  "IDL Distribution" ,"VLDL Distribution"  ))

  #are the sub-fraction plots there
  expect_contains(object = names(test[["pieCharts"]]), expected = c("LDL Subfraction",  "HDL Subfraction",  "VLDL Subfraction"))
  expect_contains(object = names(test[["pieCharts"]][["HDL Subfraction"]]), expected = c("TG", "CH", "FC", "CE", "PL"))
  
  })



test_that("works without subfractions",{
  
  load("~/git/phenological/mva-plots/data/lipoData.rda")
  df <- lipoData[,1:112]
  test <- lipoPieChart(data = df, group = lipoData$category, subfractions = F)
  
  expect_equal(object = length(names(test[["tables"]])), expected = 1)
  expect_false(object = "HDL Subfraction" %in% names(test$pieCharts))
})

test_that("three groups works", {
  load("~/git/phenological/mva-plots/data/lipoData.rda")
  lipoData$category[sample(nrow(lipoData),35)]<-"C"
  df <- lipoData[,1:112]
  test <- lipoPieChart(data = df, group = lipoData$category, subfractions = F)
  groups <- unique(test[["pieCharts"]][["main"]][["Lipoprotein Composition"]][["data"]][["group"]])
  
  expect_contains(object = groups, expected = c("A", "B", "C"))
  
  groups <- colnames(test[["tables"]][["main"]][["Lipoprotein Composition"]])
  expect_contains(object = groups, expected = c("C-A", "C-B", "B-A"))
})

test_that("multiple cohorts works", {
  load("~/git/phenological/mva-plots/data/lipoData.rda")
  lipoData$cohort<- rep_len(x = "Aus", length.out = nrow(lipoData))
  lipoData2 <- lipoData
  lipoData2$cohort<- rep_len(x = "USA", length.out = nrow(lipoData))
  df <- rbind(lipoData, lipoData2)
  # group <- df$category
  # cohort <- df$cohort
  
  test <- lipoPieChart(data = df[,1:112], group = df$category, cohort = df$cohort, subfractions = F)
  
  
  
})

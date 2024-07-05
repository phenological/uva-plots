
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
  
  
  
  
  actual <- colnames(perc)
  expected = c('HDCE', 'VLCE', 'IDCE', 'LDCE', 'VLPN', 'IDPN', 'H1CE', 'H2CE',
               'H3CE', 'H4CE', 'V1CE', 'V2CE', 'V3CE', 'V4CE', 'V5CE', 'L1CE',
               'L2CE', 'L3CE', 'L4CE', 'L5CE', 'L6CE', 'HDTG', 'HDCH', 'HDFC',
               'HDPL', 'VLTG', 'VLCH', 'VLFC', 'VLPL', 'IDTG', 'IDCH', 'IDFC',
               'IDPL', 'LDTG', 'LDCH', 'LDFC', 'LDPL', 'H1A1', 'H2A1', 'H3A1',
               'H4A1', 'H1A2', 'H2A2', 'H3A2', 'H4A2', 'L1AB', 'L2AB', 'L3AB',
               'L4AB', 'L5AB', 'L6AB', 'L1PN', 'L2PN', 'L3PN', 'L4PN', 'L5PN',
               'L6PN', 'H1TG', 'H2TG', 'H3TG', 'H4TG', 'H1CH', 'H2CH', 'H3CH',
               'H4CH', 'H1FC', 'H2FC', 'H3FC', 'H4FC', 'H1PL', 'H2PL', 'H3PL',
               'H4PL', 'V1TG', 'V2TG', 'V3TG', 'V4TG', 'V5TG', 'V1CH', 'V2CH',
               'V3CH', 'V4CH', 'V5CH', 'V1FC', 'V2FC', 'V3FC', 'V4FC', 'V5FC',
               'V1PL', 'V2PL', 'V3PL', 'V4PL', 'V5PL', 'L1TG', 'L2TG', 'L3TG',
               'L4TG', 'L5TG', 'L6TG', 'L1CH', 'L2CH', 'L3CH', 'L4CH', 'L5CH',
               'L6CH', 'L1FC', 'L2FC', 'L3FC', 'L4FC', 'L5FC', 'L6FC', 'L1PL',
               'L2PL', 'L3PL', 'L4PL', 'L5PL', 'L6PL')

  same <- setequal(actual, expected)

  expect_true(object = same)

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


lipoData$cohort <- as.factor(as.numeric(as.factor(lipoData$category)))
unique_factors <- unique(lipoData$cohort)
# perc$cohort <- lipoData$cohort
# calc$cohort <- lipoData$cohort
# 
# perc[is.na(perc)] <- 0
# perc[sapply(perc, is.infinite)] <- 0
# 
# calc[is.na(calc)] <- 0
# calc[sapply(calc, is.infinite)] <- 0

all <- cbind(calc, perc)
all$cohort <- lipoData$cohort

all[is.na(all)] <- 0
all[sapply(all, is.infinite)] <- 0
########tables######
tableCombos = list(`Lipoprotein Composition` = c("HDCE", "HDTL", "IDCE", "IDTL", "LDCE", "LDTL", "TBPN", "VLCE", "VLTL"),
                   `particle numbers` = c("TBPN", "VLPN","IDPN","L1PN","L2PN","L3PN","L4PN","L5PN","L6PN"),
                   `HDL` = c("HDTG", "HDCH", "HDFC", "HDCE_perc", "HDPL"),
                   `LDL` = c("LDCE_perc", "LDCH", "LDFC", "LDPL", "LDTG"),
                   `IDL` = c("IDCE_perc", "IDCH", "IDFC", "IDPL", "IDPN", "IDTG"),
                   `VLDL` = c("VLCE_perc", "VLCH", "VLFC", "VLPL", "VLPN", "VLTG")
                   )

# cat(paste(shQuote((VLDL[["PL"]]), type = "cmd"), collapse=", "))
#all from perc

#subfractions
subfractions <- TRUE
if(subfractions == T){
  tableCombos = list(`main` = tableCombos,
                     `LDL` = list(`TG` = c("L1TG", "L2TG", "L3TG", "L4TG", "L5TG", "L6TG", "LDTG"),
                                  `CH` = c("L1CH", "L2CH", "L3CH", "L4CH", "L5CH", "L6CH", "LDCH"),
                                  `FC` = c("L1FC", "L2FC", "L3FC", "L4FC", "L5FC", "L6FC", "LDFC"), 
                                  `CE` = c("L1CE", "L2CE", "L3CE", "L4CE", "L5CE", "L6CE", "LDCE_perc"),
                                  `PL` = c("L1PL", "L2PL", "L3PL", "L4PL", "L5PL", "L6PL", "LDPL")),
                     `HDL` = list(`TG` = c("H1TG", "H2TG", "H3TG", "H4TG", "HDTG"), 
                                  `CH` = c("H1CH", "H2CH", "H3CH", "H4CH", "HDCH"), 
                                  `FC` = c("H1FC", "H2FC", "H3FC", "H4FC", "HDFC"), 
                                  `CE` = c("H1CE", "H2CE", "H3CE", "H4CE", "HDCE_perc"), 
                                  `PL` = c("H1PL", "H2PL", "H3PL", "H4PL", "HDPL")), 
                     `VLDL` = list(`TG` = c("V1TG", "V2TG", "V3TG", "V4TG", "V5TG", "VLTG"), 
                                   `CH` = c("V1CH", "V2CH", "V3CH", "V4CH", "V5CH", "VLCH"), 
                                   `FC` = c("V1FC", "V2FC", "V3FC", "V4FC", "V5FC", "VLFC"), 
                                   `CE` = c("V1CE", "V2CE", "V3CE", "V4CE", "V5CE", "VLCE_perc"), 
                                   `PL` = c("V1PL", "V2PL", "V3PL", "V4PL", "V5PL", "VLPL")))
}



#set up 
if(length(unique_factors) == 2){
  col <- 3
} else{
  col <- 2 + (sum(1:length(unique_factors))-length(unique_factors))
}

combinations <- combn(unique_factors, 2)

# Create column names from the combinations with the larger letter first
column_names <- apply(combinations, 2, function(x) paste(sort(x, decreasing = TRUE), collapse = "-"))


# Define the lipoproteins vector
# lipoproteins <- c("HDTL", "HDCE", "VLTL", "VLCE", "IDTL", "IDCE", "LDTL", "LDCE", "TBPN") #"HDA1", "HDA2", "LDAB")

trial <- list()
for(j in names(tableCombos)){
  for(i in names(tableCombos[[j]])){
    lipoproteins <- tableCombos[[j]][[i]]
    
    sig <- as.data.frame(matrix(NA, nrow = length(lipoproteins), ncol = col))
    
    colnames(sig) <- c("lipoproteins", "pval", column_names)
    sig$lipoproteins <- lipoproteins
    
    for(lipo in lipoproteins){
      anova_model <- aov(as.formula(paste(lipo, "~ cohort")), data = all)
      tukey_results <- TukeyHSD(anova_model)
      # Extract p-values
      idx <- which(sig$lipoproteins == lipo)
      sig[idx, "pval"] <- as.numeric(summary(anova_model)[[1]]$`Pr(>F)`[1])
      tukey_pvalues <- tukey_results$cohort
      
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
    mapping <- setNames(unique(lipoData$category), unique(lipoData$cohort))
    
    testnames<- as.data.frame(mapping, check.names = F)
    testnames$rowName <- rownames(testnames)
    
    new_colnames <- colnames(sig)
    for (num in names(mapping)) {
      new_colnames <- gsub(num, mapping[num], new_colnames)
    }
    colnames(sig) <- new_colnames
    
    # Calculate mean and sd for each 'cohort' group and each 'calc' column
    agg_mean <- melt(aggregate(. ~ cohort, data = all, FUN = mean), id.vars = "cohort")
    agg_sd <- melt(aggregate(. ~ cohort, data = all, FUN = sd), id.vars = "cohort")
    
    new <- merge(x = agg_mean, y = agg_sd, by = c("cohort", "variable"))
    new$`mean±sd` <- paste0(round(new$value.x, 2), " ± ", round(new$value.y, 2))
    new <- new[,which(!(colnames(new) %in% c("value.x", "value.y")))]
    
    result <- dcast(new, variable ~ cohort, value.var = "mean±sd")
    #rename using mapping
    new_colnames <- colnames(result)
    for (num in names(mapping)) {
      new_colnames <- gsub(num, mapping[num], new_colnames)
    }
    colnames(result) <- new_colnames
    
    tble <- merge(result, sig, by.x = "variable", by.y = "lipoproteins")
    
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
#     anova_model <- aov(as.formula(paste(lipo, "~ cohort")), data = all)
#     tukey_results <- TukeyHSD(anova_model)
#     # Extract p-values
#     idx <- which(sig$lipoproteins == lipo)
#     sig[idx, "pval"] <- as.numeric(summary(anova_model)[[1]]$`Pr(>F)`[1])
#     tukey_pvalues <- tukey_results$cohort
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
#   mapping <- setNames(unique(lipoData$category), unique(lipoData$cohort))
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
#   # Calculate mean and sd for each 'cohort' group and each 'calc' column
#   agg_mean <- melt(aggregate(. ~ cohort, data = all, FUN = mean), id.vars = "cohort")
#   agg_sd <- melt(aggregate(. ~ cohort, data = all, FUN = sd), id.vars = "cohort")
#   
#   new <- merge(x = agg_mean, y = agg_sd, by = c("cohort", "variable"))
#   new$`mean±sd` <- paste0(round(new$value.x, 2), " ± ", round(new$value.y, 2))
#   new <- new[,which(!(colnames(new) %in% c("value.x", "value.y")))]
#   
#   result <- dcast(new, variable ~ cohort, value.var = "mean±sd")
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
                   `particle numbers` = c("VLPN","IDPN","L1PN","L2PN","L3PN","L4PN","L5PN","L6PN"),
                   `HDL` = c("HDTG" , "HDFC", "HDCE_perc" , "HDPL" ),
                   `LDL` = c("LDCE_perc" , "LDFC" , "LDPL", "LDTG" ),
                   `IDL` = c("IDCE_perc" , "IDFC", "IDPL" , "IDTG"),
                   `VLDL` = c("VLCE_perc", "VLFC", "VLPL", "VLTG"))

newnames <- list(
  `Lipoprotein Composition` = c("HDTL" = "HDL", "IDTL" = "IDL", "LDTL" = "LDL", "VLTL" = "VLDL"),
  `particle numbers` = c("VLPN" = "VLDL", "IDPN" = "IDL", "L1PN" = "LDL1", "L2PN" = "LDL2", "L3PN" = "LDL3", "L4PN" = "LDL4", "L5PN" = "LDL5", "L6PN" = "LDL6"),
  `HDL` = c("HDTG" = "TG", "HDFC" = "FC", "HDCE_perc" = "CE", "HDPL" = "PL"),
  `LDL` = c("LDCE_perc" = "CE", "LDFC" = "FC", "LDPL" = "PL", "LDTG" = "TG"),
  `IDL` = c("IDCE_perc" = "CE", "IDFC" = "FC", "IDPL" ="PL", "IDTG" ="TG"),
  `VLDL` = c("VLCE_perc" = "CE", "VLFC" = "FC", "VLPL" = "PL", "VLTG" = "TG")
)

plotComboSubf <- list(`LDL` = list(`TG` = c("L1TG", "L2TG", "L3TG", "L4TG", "L5TG", "L6TG", "LDTG"),
                                   `CH` = c("L1CH", "L2CH", "L3CH", "L4CH", "L5CH", "L6CH", "LDCH"),
                                   `FC` = c("L1FC", "L2FC", "L3FC", "L4FC", "L5FC", "L6FC", "LDFC"), 
                                   `CE` = c("L1CE", "L2CE", "L3CE", "L4CE", "L5CE", "L6CE", "LDCE_perc"),
                                   `PL` = c("L1PL", "L2PL", "L3PL", "L4PL", "L5PL", "L6PL", "LDPL")),
                      `HDL` = list(`TG` = c("H1TG", "H2TG", "H3TG", "H4TG", "HDTG"), 
                                   `CH` = c("H1CH", "H2CH", "H3CH", "H4CH", "HDCH"), 
                                   `FC` = c("H1FC", "H2FC", "H3FC", "H4FC", "HDFC"), 
                                   `CE` = c("H1CE", "H2CE", "H3CE", "H4CE", "HDCE_perc"), 
                                   `PL` = c("H1PL", "H2PL", "H3PL", "H4PL", "HDPL")), 
                      `VLDL` = list(`TG` = c("V1TG", "V2TG", "V3TG", "V4TG", "V5TG", "VLTG"), 
                                    `CH` = c("V1CH", "V2CH", "V3CH", "V4CH", "V5CH", "VLCH"), 
                                    `FC` = c("V1FC", "V2FC", "V3FC", "V4FC", "V5FC", "VLFC"), 
                                    `CE` = c("V1CE", "V2CE", "V3CE", "V4CE", "V5CE", "VLCE_perc"), 
                                    `PL` = c("V1PL", "V2PL", "V3PL", "V4PL", "V5PL", "VLPL")))

#graphCombos
# `Lipoprotein Composition` = c("HDTL","VLTL","IDTL","LDTL")


lipoproteins <- c(plotCombos[["Lipoprotein Composition"]], "cohort")
medians <- aggregate(. ~ cohort, data = all, function(x) median(x, na.rm = TRUE))

medians <- aggregate(. ~ cohort, data = all[,lipoproteins], function(x) median(x, na.rm = TRUE))
medians$total <- rowSums(medians[, -which(names(medians) == "cohort")])

for (col in names(medians)[-which(names(medians) %in% c("cohort", "total"))]) {
  medians[[col]] <- round((medians[[col]] / medians$total) * 100, 2)
}

medians<-
medians[,-which(names(medians) %in% c("total"))]

category <- `Lipoprotein Composition`
  for (category in names(newnames)) {
    for (oldname in names(newnames[[category]])) {
      newname <- newnames[[category]][oldname]
      if (oldname %in% colnames(medians)) {
        colnames(medians)[colnames(medians) == oldname] <- newname
      }
    }
  }

# Reshape the data from wide to long format
long_data <- reshape(
  medians,
  varying = list(names(medians[,-which(names(medians) %in% c("cohort"))])),
  v.names = "y",
  idvar = "cohort",
  times = names(medians)[2:5],
  timevar = "subfraction",
  direction = "long"
)

# Remove the row.names column added by reshape
long_data <- long_data[order(long_data$cohort), ]
long_data$labels <- paste0(long_data$y, " %")
long_data$subfraction <- as.factor(long_data$subfraction)
long_data$label_pos <- NA

#   calculate_label_pos <- function(df) {
#     df$label_pos <- cumsum(df$y) - 0.5 * df$y
#     return(df)
#   }
# 
# # Split the data by cohort, apply the function, and combine the results
# long_data <- do.call(rbind, lapply(split(long_data, long_data$cohort), calculate_label_pos))

long_data <- 
  do.call(what = rbind, 
          args = lapply(split(long_data, long_data$cohort), function(df) {
  df$label_pos <- cumsum(df$y) - 0.5 * df$y
  return(df)
}))

piePlot <-
ggplot(data = long_data, aes(x = "", y = y, fill = subfraction)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  facet_wrap(~cohort)+
  labs(title = "Lipoprotein subfraction distribution")






dat1%>%
  select(cohort,starts_with("perc.HD"))%>%
  rename_with( ~ gsub("perc.HD", "", .x, fixed = TRUE))%>%
  select(!"CH")%>%
  mutate(cohort = factor(cohort,levels = c("Control","COVID","MISC")))%>%
  dplyr::group_by(cohort)%>%
  summarise(
    med_TG = median(TG,na.rm = TRUE),
    med_FC = median(FC,na.rm = TRUE),
    med_CE = median(CE,na.rm = TRUE),
    med_PL = median(PL,na.rm = TRUE),
    total= med_TG+med_FC+med_CE+med_PL
  )%>%
  mutate(
    TG = round((med_TG/total)*100,2),
    FC = round((med_FC/total)*100,2),
    CE = round((med_CE/total)*100,2),
    PL = round((med_PL/total)*100,2),
  )%>%
  select(cohort,TG:PL)
trial <- calc[,which(grepl(pattern = "TL", x = colnames(calc))|colnames(calc) == "cohort")]

colnames(trial)[colnames(trial) == "VLTL"] <- "VLDL"

colnames(trial) <- gsub("T", "", colnames(trial))

agg_med <-  melt(aggregate(. ~ cohort, data = trial, FUN = median), id.vars = "cohort")
agg_med <- dcast(agg_med, variable ~ cohort, value.var = "value")

idx <- nrow(agg_med)+1
agg_med$variable<- as.character(agg_med$variable)
agg_med[5, "variable"] <- "total"

for(i in 1:length(unique_factors)){
  agg_med[idx, i+1] <- sum(na.omit(agg_med[, i+1]))
  
  # agg_med[idx, which(colnames(agg_med) == as.character(i))] <- sum(na.omit(agg_med[, which(colnames(agg_med) == as.character(i))]))
}

total_row <- agg_med[agg_med$variable == "total", ]

# Divide each entry in columns A and B by the corresponding total values and round to 2 decimal places
for(i in 1:length(unique_factors)){
  for(j in 1:nrow(agg_med)){
    agg_med[j,i+1] <- 
      round(agg_med[j,i+1] / total_row[i+1]*100, 2)
  }
}


agg_med <- agg_med[-idx,]
# agg_med$variable <- as.character(agg_med$variable)


###############
agg_med <- 
reshape(agg_med, 
        varying = c(2:ncol(agg_med)), 
        v.names = "value", 
        timevar = "cohort", 
        # times = c(1:length(unique_factors)), 
        idvar = "variable",
        direction = "long")

#rename using mapping
agg_med$cohort <- mapping[as.character(agg_med$cohort)]


# agg_med$label <- paste0(agg_med$value, "%")
plot <- 
ggplot(data = agg_med, 
       aes(x = "", 
           y = value, 
           fill = variable)) +
  geom_bar(stat = "identity", 
           width = 1, 
           color = "white") +
  coord_polar("y", 
              start = 0) +
  theme_void() +
  facet_wrap(~ cohort)+
  labs(title = "Lipoprotein subfraction distribution")




 perc[,which(grepl("PN", colnames(perc))| colnames(perc) == "cohort")]


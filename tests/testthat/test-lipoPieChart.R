
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


############etablished test##############

test_that("works with subfractions and subcompositions",{
  # lipoData <- mva.plots::lipoData
  # load("~/git/phenological/mva-plots/data/lipoData.rda")
  rmlg <- lipoData[,1:112]
  test <- lipoPieChart(data = rmlg, group = lipoData$category, subfractions = T, subcompositions = T)
  
  expect_equal(object = length(names(test[["tables"]])), expected = 7)
  
  #are the main composition plots there
  expect_contains(object = names(test[["pieCharts"]][["main composition"]]), expected = c("Lipoprotein Composition", "Particle Numbers", "HDL Distribution",   "LDL Distribution",  "IDL Distribution" ,"VLDL Distribution"  ))

  #are the sub-fraction and sub-composition plots there
  expect_contains(object = names(test[["pieCharts"]]), expected = c("LDL Subfraction",  "HDL Subfraction",  "VLDL Subfraction", "LDL Subcomposition",  "HDL Subcomposition",  "VLDL Subcomposition"))
  expect_contains(object = names(test[["pieCharts"]][["HDL Subfraction"]]), expected = c("TG", "CH", "FC", "CE", "PL"))
  expect_contains(object = names(test[["pieCharts"]][["HDL Subcomposition"]]), expected = c("H1", "H2", "H3", "H4"))
  
  })



test_that("works without subfractions",{
  # lipoData <- mva.plots::lipoData
  # load("~/git/phenological/mva-plots/data/lipoData.rda")
  df <- lipoData[,1:112]
  test <- lipoPieChart(data = df, group = lipoData$category, subfractions = F, subcompositions = F)
  
  expect_equal(object = length(names(test[["tables"]])), expected = 1)
  expect_false(object = "HDL Subfraction" %in% names(test$pieCharts))
})

test_that("three groups works", {
  # lipoData <- mva.plots::lipoData
  # load("~/git/phenological/mva-plots/data/lipoData.rda")
  lipoData$category[sample(nrow(lipoData),35)]<-"C"
  df <- lipoData[,1:112]
  test <- lipoPieChart(data = df, group = lipoData$category, subfractions = F, subcompositions = F)
  groups <- unique(test[["pieCharts"]][["main composition"]][["Lipoprotein Composition"]][["data"]][["group"]])
  
  expect_contains(object = groups, expected = c("A", "B", "C"))
  
  groups <- colnames(test[["tables"]][["main composition"]][["Lipoprotein Composition"]])
  expect_contains(object = groups, expected = c("C-A", "C-B", "B-A"))
})

test_that("multiple cohorts works", {
  # load("~/git/phenological/mva-plots/data/lipoData.rda")
  # lipoData <- mva.plots::lipoData
  lipoData$cohort<- rep_len(x = "Aus", length.out = nrow(lipoData))
  lipoData2 <- lipoData
  lipoData2$cohort<- rep_len(x = "USA", length.out = nrow(lipoData))
  df <- rbind(lipoData, lipoData2)
  
  test <- lipoPieChart(data = df[,1:112], group = df$category, cohort = df$cohort, subfractions = F)

  cohorts <- levels(test[["pieCharts"]][["main composition"]][["Lipoprotein Composition"]][["data"]][["cohort"]])
  
  expect_equal(cohorts, c("Aus", "USA"))
  
})





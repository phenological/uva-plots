pieChart <- function(data, cohort, multiCohort = F){

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




  data<-data[,which(colnames(data) %in% lipo_name)]
  tdf<-apply(data,1,as.list)
  tdf1<-lapply(tdf, as.data.frame)
  melt_and_rename <- function(df) {
    df_long  = stack(df)
    colnames(df_long) = c("value","id")
    return(df_long)
  }
  tdf2 <- lapply(tdf1, melt_and_rename)
  data<-lapply(tdf2, function(x) extend_lipo(x))
  data<-data.frame(do.call("rbind",data))
  return(data)

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
 

}

#' Title
#'
#' @param interpro
#' @param gene_list
#' @param gene_length
#' @param domain_db
#' @param motif_db
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
interpro2draw <- function(interpro, gene_list,
                          gene_length, domain_db = "Pfam",
                          motif_db = "TMHMM"){
  db_list <- list(domain_db, motif_db)

  drawobj <- interpro %>%
    filter(feature_id %in% gene_list) %>%
    mutate(begin = start, description = sigunature_desc) %>%
    filter(DB %in% db_list) %>%
    mutate(type = ifelse(DB == domain_db, "DOMAIN", "MOTIF")) %>%
    mutate(entryName = feature_id) %>%
    mutate(order = 1) %>%
    select(type, description, begin, end, length, entryName, order)
  return(drawobj)
}

#' Title
#'
#' @param gene_length
#' @param gene_list
#'
#' @return
#' @export
#'
#' @import tibble
#' @import dplyr
#'
#' @examples
lenlist2chain <- function(gene_length, gene_list){
  genelen_hit <- gene_length %>%
    filter(protein_id %in% gene_list)

  genenum <- nrow(genelen_hit)

  chain_data <- tibble(type = rep("CHAIN", genenum),
                       description = rep(NA, genenum),
                       begin = rep(1, genenum),
                       end = genelen_hit$length,
                       length = genelen_hit $length,
                       entryName = genelen_hit$protein_id)
  return(chain_data)
}

#' Title
#'
#' @param domain_data
#' @param chain_data
#' @param splice_data
#'
#' @return
#' @export
#'
#' @import tibble
#' @import dplyr
#'
#' @examples
bind_draws <- function(domain_data, chain_data, splice_data = NA){
  all_draw <- domain_data %>%
    add_row(chain_data)

  if(!is.na(splice_data)){
    all_draw <- all_draw %>%
      add_row(splice_data) %>%
      mutate(entryName_psuedo = str_replace(entryName, pattern=".m1$", replacement="")) %>%
      mutate(.,order = group_indices(.,entryName_psuedo)) %>%
      mutate(entryName = entryName_psuedo)
  } else {
    all_draw <- all_draw %>%
      mutate(entryName_psuedo = str_replace(entryName, pattern=".m1$", replacement="")) %>%
      mutate(.,order = group_indices(.,entryName_psuedo)) %>%
      mutate(entryName = entryName_psuedo)
  }
  return(all_draw)
}

#' Extract the genes corresponding to the domains from the inteproscan results
#'
#' @param interpro
#' @param domain
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
domain_extract <- function(interpro, domain){
  interpro %>%
    filter(ID == domain) %>%
    select(feature_id) %>%
    pull()
}

#' Extract the reference genes corresponding to the query genes from the blast results
#'
#' @param blast_data
#' @param bacteria_gene
#' @param evalues
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
blast_extract <- function(blast_data, bacteria_gene, evalues = 1*10^-10){
  host_hit <- blast_data %>%
    filter(qaccver %in% bacteria_gene & evalue < evalues) %>%
    select(saccver) %>%
    unique() %>%
    pull()
  return(host_hit)
}

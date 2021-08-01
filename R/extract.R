#' Extract the genes corresponding to the domains from the inteproscan results
#'
#'
#'
#' @param interpro interproscan results from tsv output files using read_tsv
#' @param domain Pfam ID you want to specify
#' @return genelist
#'
#' @import dplyr
#'
#' @examples
#'
#' @export
domain_extract <- function(interpro, domain){
  interpro %>%
    filter(ID == domain) %>%
    select(feature_id) %>%
    pull()
}

#' Extract the reference genes corresponding to the query genes from the blast results
#'
#'
#'
#' @param blast_data a blast tibble using read_blast function
#' @param bacteria_gene List of bacteria genes (result of domain_extract)
#' @param evalues threshold for evalues to filter (default = 1*10^-10)
#' @return refenrece gene list coressponding to query genes
#'
#' @import dplyr
#'
#' @examples
#'
#' @export
blast_extract <- function(blast_data, bacteria_gene, evalues = 1*10^-10){
  host_hit <- blast_data %>%
    filter(qaccver %in% bacteria_gene & evalue < evalues) %>%
    select(saccver) %>%
    unique() %>%
    pull()
  return(host_hit)
}

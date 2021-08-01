#' Convert gff to gggenomes data
#'
#' @param gff_data
#' @param gene_list
#' @param aten
#'
#' @return
#' @export
#'
#' @examples
rename_seq <- function(gff_data, gene_list, aten = TRUE){
  renamed_genes <- gff_data %>%
    dplyr::filter(parent_ids %in% gene_list |
                    feat_id %in% gene_list) %>% #gffから特定の遺伝子取得
    arrange(seq_id) %>% #配列順にソート
    {if(aten) mutate(seq_ids_psuedo = ifelse(grepl(".m1$",parent_ids), #mRNAは末端にm1とついているので、
                                          str_replace(parent_ids,  #あるものは除去する
                                                      pattern=".m1$", replacement=""),
                                          parent_ids))
      else .} %>%
    mutate(seq_id = seq_ids_psuedo) #配列名を再定義
  return(renamed_genes)
}

#' Extract cds data from gene datafram
#'
#' @param gene_data
#'
#' @return
#' @export
#'
#' @examples
gene2cds <- function(gene_data){
  cds <- gene_data %>%
    dplyr::filter(type == "CDS") %>%
    mutate(length = abs(start - end)) %>%
    mutate(cds_length = ifelse(length %% 3 > 0, length %/% 3 + 1, length %/% 3))
  return(cds)
}

#' Extract splicing point from cds data
#'
#' @param cds_data
#'
#' @return
#' @export
#'
#' @examples
cds2sp <- function(cds_data){
  splicing <- cds_data %>%
    group_by(seq_id) %>%
    nest() %>%
    mutate(spliced = map(data, splice_point)) %>%
    select(seq_id, spliced) %>%
    unnest(cols = spliced) %>%
    mutate(entryName = seq_id) %>%
    ungroup() %>%
    dplyr::select(-seq_id)
  return(splicing)
}

splice_point <- function(data){
  point <- c()
  num <- nrow(data)

  if(num != 1){
    #    if(data$strand[1] == "-"){
    for (i in 1:(num - 1)){
      point <- c(point, sum(data$cds_length[1:i]))
    }
    #    } else {
    #      for (i in num:2){
    #        point <- c(point, sum(data$cds_length[num:i]))
    #      }
    #    }
  }
  #以下の部分はdrawProteins様に成型している
  out <- tibble(type = rep("MOD_RES", num - 1),
                description = rep("splicing point", num - 1),
                begin = point,
                end = point,
                length = rep(NA, num - 1),
                order = rep(NA, num - 1)
  )
  return(out)
}

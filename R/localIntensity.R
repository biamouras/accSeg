#' Processes the local intensity for each locality
#'
#' @param matrix data frame (origin, destination, weight)
#' @param pop data frame (id, group, population structure)
#' @return data frame with local intensity for each group
#' @import dplyr

localIntensity <- function(matrix, pop){
  matrix %>%
    dplyr::right_join(pop, by=c('destination'='id')) %>%
    dplyr::select(-.data$destination) %>%
    dplyr::mutate(int = .data$weight * .data$population) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(local_int = sum(.data$int, na.rm=T)/sum(.data$weight, na.rm=T)) %>%
    dplyr::ungroup()
}

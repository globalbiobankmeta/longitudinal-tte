# PheWAS::mapPhecodesToExclusions() throws the error:
# Error in View : `tbl_df()` was deprecated in dplyr 1.0.0 and is now defunct.
# Please use `tibble::as_tibble()` instead.

# So this copy of the function replaces tbl_df() with tibble::as_tibble()

mapPhecodesToExclusions <- function (phecodes, ids) 
{
  if (missing(ids)) {
    input = tibble::as_tibble(data.frame(id = 0, exclusion_criteria = phecodes, 
                              stringsAsFactors = FALSE))
  }
  else {
    input = tibble::as_tibble(data.frame(id = ids, exclusion_criteria = phecodes, 
                              stringsAsFactors = FALSE))
  }
  output = inner_join(input, phecode_exclude)
  output = output %>% transmute(id, exclusion = code) %>% distinct()
  output
}

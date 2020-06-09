# Define geometricly centred quantiles function
gm_quantiles <- function (x, q = c(0.25, 0.5, 0.75), na_rm = FALSE, wide=FALSE) {
  #tibble(x = exp(quantile(log(x), q,  na.rm = na_rm)), q=q)
  out <- tibble("{{ x }}" := exp(quantile(log(x), q,  na.rm = na_rm)), "{{ x }}_q" := q)
  if(wide){
    out <- out %>%
      pivot_wider(names_from = 2,
                  values_from=1) %>%
      rename_with(~gsub("0.", "q", .x, fixed = TRUE), starts_with("0."))
  }
  return(out)
}
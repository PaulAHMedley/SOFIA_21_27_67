# Standard Plotting and Table Functions

library("officer")
library("ggplot2")
library("cowplot")


SRA_PrepareCatchData <- function(assess_species) {

  ass_df <- tc_df |>
    filter(Name_En %in% assess_species) |>
    select(stock = Name_En, taxa = Scientific_Name, year=PERIOD, capture=VALUE)
  # add in zeroes for missing years
  df <- expand_grid(stock=assess_species, year=min(tc_df$PERIOD):max(tc_df$PERIOD)) |>
    mutate(capture=0) |>
    left_join(summarise(group_by(ass_df, stock, taxa)), by = "stock")

  ass_df <- rbind(ass_df, df) |>
    group_by(taxa, stock, year) |>
    summarise(capture=sum(capture)) |>
    ungroup()

  return(ass_df %>%
    group_by(stock, taxa) %>%
    nest() %>%
    ungroup())
}

SRA_fit_catchonly <- function(ISSCAAP, stocks) {

  filename = here("output", paste0("area", as.character(area), "_", gsub("([, ])", "", ISSCAAP),  ".rda"))

  if (file.exists(filename)) {
    load(filename)
  } else {
    sfs <- safely(fit_sraplus)
    samps <- nrow(stocks)

    stocks <- stocks |>
      ungroup() |>
      sample_n(samps) |>
      mutate(sraplus_fit = map(
        driors,
        ~ sfs(
          driors = .x,
          engine = "stan",
          model = "sraplus_tmb",
          adapt_delta = 0.95,
          max_treedepth = 10,
          n_keep = 1000,
          chains = 4,
          cores = 4,
          ## q_slope prior=,
          estimate_qslope = FALSE,
          estimate_proc_error = TRUE,
          cleanup = TRUE)
      ))
    stocks <- stocks |>
      mutate(sraplus_worked = map_lgl(map(sraplus_fit,"error"), is.null)) |>
      filter(sraplus_worked) |>
      mutate(sraplus_fit = map(sraplus_fit,"result")) |>
      mutate(sraplus_summary = map(sraplus_fit, summarize_sralpus)) |>
      mutate(sraplus_diagnostics = map2(sraplus_fit, driors, diagnose_sraplus)) |>
      arrange(stock)

    save(stocks, file = here("output", paste0("area", as.character(area), "_", gsub("([, ])", "", ISSCAAP),  ".rda")))
  }
  return(stocks)
}

SRA_results <- function(stocks) {

  StockAssessments <- SA_df |>
    filter(`Scientific Name` %in% stocks$taxa) |>
    mutate(Weight = ifelse(Status=="U", NA, Weight)) |>
    group_by(`Scientific Name`) |>
    summarise(
      `Stocks Scored` = n(),
      Weight = sum(Weight, na.rm=TRUE),
      score_wt = sum(score_wt, na.rm=TRUE)) |>
    ungroup() |>
    mutate(`Average Score` = score_wt/Weight) |>
    mutate(`True Status` = StatusCode(`Average Score`)) |>
    select(taxa = `Scientific Name`, `True Status`) |>
    replace_na(list(`True Status`="U"))

  current_status <- stocks |>
    select(taxa, stock, sraplus_summary) |>
    left_join(StockAssessments, by="taxa") |>
    unnest(cols = sraplus_summary) |>
    filter(variable == "b_div_bmsy") |>
    mutate(`SRA Status` = case_when(mean > 1.9 ~ "U",
                                    mean > 1.2 ~ "N",
                                    mean > 0.8 ~ "F",
                                    TRUE ~ "O")) |>
    select(Stock=stock, Year=year, `B/BMSY` = mean, Lower=lower, Upper=upper, `SRA Status`, `True Status`)
  return(current_status)
}

SRA_table_results <- function(current_status) {
  tab <- current_status |>
    flextable() |>
    colformat_num(j=2, big.mark = "") |>
    colformat_double(j=3:5, digits = 3) |>
    align(j = 6:7, align = "center") |>
    autofit()
}


StatusCode <- function(AvgScore) {
  return(case_when(is.na(AvgScore) ~ "U",
                   AvgScore <= 1.25 ~ "O",
                   AvgScore <= 1.75 ~ "O-F",
                   AvgScore <= 2.25 ~ "F",
                   TRUE ~ "N"))
}



table_stock_list <- function(tISSCAAP) {
  return(SA_df |>
           filter(ISSCAAP == tISSCAAP) |>
           select(Jurisdiction=Country, `Stock Name`, Status) |>
           flextable() |>
           autofit() |>
           width(j = 2, width=11, unit = "cm") )
}


table_stock_area  <- function(tISSCAAP) {
  return(SA_df |>
           filter(ISSCAAP == tISSCAAP) |>
           select(`Stock Name`, Area = `Stock Area`, Status) |>
           flextable() |>
           autofit() |>
           width(j = 1, width=9, unit = "cm") )
}


table_ISSCAAP_scores <- function() {
  SA_summary <-  SA_df |>
    filter(ISSCAAP %in% c(main_ISSCAAP, secondary_ISSCAAP),
           Status != "U",
           ! is.na(Weight))
  return(
    rbind(
      SA_summary |>
        group_by(ISSCAAP) |>
        summarise(
          `Stocks Scored` = n(),
          Weight = sum(Weight),
          score_wt = sum(score_wt)) |>
        ungroup() |>
        mutate(`Average Score` = score_wt/Weight) |>
        mutate(`Average Status` = StatusCode(`Average Score`)) |>
        select(ISSCAAP, `Stocks Scored`, `Average Score`, `Average Status`),
      tibble(ISSCAAP = "Grand Total",
             `Stocks Scored` = nrow(SA_summary),
             `Average Score` = sum(SA_summary$score_wt)/sum(SA_summary$Weight),
             `Average Status` = StatusCode(sum(SA_summary$score_wt)/sum(SA_summary$Weight)))
    ) |>
      flextable() |>
      colformat_double(j=3, digits=2) |>
      align(j = 4, align = "center") |>
      border(part = "body", i = ~ ISSCAAP == "Grand Total", border.top = fp_border(color="grey40", width = 2) ) |>
      bold(part = "body", i = ~ ISSCAAP == "Grand Total") |>
      autofit()
  )
}


table_ISSCAAP_status <- function() {

  SA_summary <- SA_df |>
    group_by(ISSCAAP, Status) |>
    summarise(
      N = n()) |>
    ungroup() |>
    pivot_wider(names_from = Status, values_from=N, values_fill=0)

  if ("O-F" %in% names(SA_summary))
    SA_summary <- dplyr::select(SA_summary, ISSCAAP, O, `O-F`, `F`, N, U)
  else
    SA_summary <- dplyr::select(SA_summary, ISSCAAP, O, `F`, N, U)

  total_df <- summarise(SA_summary, across(-ISSCAAP, sum)) |>
    mutate(ISSCAAP = "Grand Total")

  return(rbind(
    SA_summary,
    total_df) |>
      flextable() |>
      colformat_double(j=3, digits=2) |>
      align(j = 4, align = "center") |>
      #  border(part = "body", i = ~ ISSCAAP == "Grand Total", border.top = fp_border(color="grey40", width = 2) ) |>
      bold(part = "body", i = ~ ISSCAAP == "Grand Total") |>
      autofit()
  )
}




table_ISSCAAP_category <- function() {

  SA_summary <- SA_df |>
    mutate(Category = as.character(Category)) |>
    group_by(ISSCAAP, Category) |>
    summarise(
      N = n()) |>
    ungroup() |>
    pivot_wider(names_from = Category, values_from=N, values_fill=0) |>
    dplyr::select(ISSCAAP, `1`, `2`, `3`, `4`, `5`, `6`)

  total_df <- summarise(SA_summary, across(-ISSCAAP, sum)) |>
    mutate(ISSCAAP = "Grand Total")

  return(rbind(
    SA_summary,
    total_df) |>
      flextable() |>
      colformat_double(j=3, digits=2) |>
      align(j = 4, align = "center") |>
      #  border(part = "body", i = ~ ISSCAAP == "Grand Total", border.top = fp_border(color="grey40", width = 2) ) |>
      bold(part = "body", i = ~ ISSCAAP == "Grand Total") |>
      autofit()
  )
}


plot_total_ISSCAAP <- function() {
  return(tc_df |>
           filter(MEASURE == "Q_tlw") |>
           mutate(Group = factor(ifelse(ISSCAAP %in% main_ISSCAAP, ISSCAAP, "Others"),
                                 levels=main_ISSCAAP)) |>
           dplyr::group_by(Group, PERIOD) |>
           summarise(VALUE=sum(VALUE)/1000000) |>
           ungroup() |>
           ggplot(aes(x=PERIOD, y=VALUE, fill=Group)) +
           geom_area() +
           labs(x="Year", y="Live Weight (million tonnes)") +
           theme_cowplot() +
           scale_fill_viridis_d() +
           theme(text = element_text(size=9),
                 axis.text = element_text(size=10)
           ))

}

plot_ISSCAAP <- function(gISSCAAP, retain_groups) {
  tc_df |>
    filter(ISSCAAP == gISSCAAP) |>
    mutate(Group = factor(ifelse(Name_En %in% retain_groups, Name_En, "Others"),
                          levels=retain_groups)) |>
    dplyr::group_by(Group, PERIOD) |>
    summarise(VALUE=sum(VALUE)/1000000) |>
    ungroup() |>
    ggplot(aes(x=PERIOD, y=VALUE, fill=Group)) +
    geom_area() +
    labs(x="Year", y="Live Weight (million tonnes)") +
    theme_cowplot() +
    scale_fill_viridis_d() +
    theme(text = element_text(size=9),
          axis.text = element_text(size=10)
    )
}

plot_relative_catch <- function(stocks) {
  return(stocks |>
    unnest(cols = c(data)) |>
    group_by(stock) |>
    mutate(capture = capture / max(capture)) |>
    ggplot(aes(year, capture, colour=stock)) +
    geom_point() +
    geom_line() +
    theme_cowplot()
  )
}

plot_others <- function(retain_groups) {
  return(tc_df |>
    filter(MEASURE == "Q_tlw",
           ! (ISSCAAP %in% main_ISSCAAP),
           (Major_Group %in% c("PISCES", "CRUSTACEA", "MOLLUSCA"))) |>
    mutate(Group = factor(ifelse(ISSCAAP %in% retain_groups, ISSCAAP, "Others"),
                          levels=retain_groups)) |>
    dplyr::group_by(Group, PERIOD) |>
    summarise(VALUE=sum(VALUE)/1000000) |>
    ungroup() |>
    ggplot(aes(x=PERIOD, y=VALUE, fill=Group)) +
    geom_area() +
    labs(x="Year", y="Live Weight (million tonnes)") +
    theme_cowplot() +
    scale_fill_viridis_d() +
    theme(text = element_text(size=9),
          axis.text = element_text(size=10)
    )
  )
}

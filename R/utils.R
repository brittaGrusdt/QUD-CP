library(boot)
library(here)

# at most one control-physics trial and none of the control-random trials wrong
# and test-example trial must be correctly replied to
clean_data = function(data){
  controls = data %>% filter(startsWith(type, "control") | type=="test-example") %>% 
    dplyr::select(c(submission_id, starts_with("id"), type, question,
                    selected_pic, expected)) %>% 
    group_by(submission_id, type) %>% 
    mutate(n=n()) %>% 
    mutate(correct = selected_pic == expected, n_correct = sum(correct)) %>% 
    distinct_at(vars(c(type, submission_id)), .keep_all = T)
  
  # failing any criteria is sufficient to be removed
  submission_ids.out = controls %>% 
    filter(!(type == "control-random" & n_correct == n) &
           !(type == "control-physics" & n_correct >= n-1) &
           !(type == "test-example" & n_correct == 1)) %>% 
    pull(submission_id) %>% unique()
  print(paste('remove data from participants:',
              paste(submission_ids.out, collapse = ", ")))
  return(data %>% filter(!submission_id %in% submission_ids.out))
}

data_critical_trials = function(cleaned_data){
  data.critical = cleaned_data %>% 
    filter(type=="critical") %>% 
    select(submission_id, id, type, question, response, expected, id1, id2,
           selected_pic, RT) %>% 
    mutate(expectation = case_when(expected == "none" ~ FALSE,
                                   TRUE ~ TRUE)) %>%
    group_by(id1, id2)
  
  picA = "A. w/oDistractor-internal"
  picB = "B. w/oDistractor-external"
  picC = "C. withDistractor-internal"
  picD = "D. withDistractor-external"
  
  data.critical$picPair = data.critical %>% group_indices()
  data.critical <- data.critical %>% 
    mutate(picPair=paste(id1, id2, sep=" & "),
           picPair.short = case_when(picPair == "if2_unn & if2_unu" ~ picD, 
                                     picPair == "if1_un & if2_unu" ~ picB,
                                     picPair == "if1_un & if1_uu" ~ picA,
                                     picPair == "if2_unn & if1_uu" ~ picC)) %>% 
    group_by(picPair, question) %>% mutate(n=n()) %>% 
    mutate(y=case_when(response == "exhaustive" ~ 1, T ~ 0)) %>%
    group_by(picPair) %>% 
    mutate(group_id = cur_group_id(), 
           pics = chartr('1234', 'ABCD', group_id)) %>% 
    group_by(pics, question) %>%
    mutate(pics = as.factor(pics), 
           question = factor(question, levels = c("neutral", "ifp", "willq"))) %>%
    rename(submissionID = submission_id) %>% 
    unite("picPair.long", "pics", "picPair", sep=". ") %>% 
    dplyr::select(-group_id) %>% 
    rename(picPair = picPair.short) %>% 
    mutate(picPair = as.factor(picPair), picPair.short = picPair) %>%
    separate(picPair.short, into=c("picPair.short", "tmp"), sep="\\.") %>% 
    mutate(picPair.short = as.factor(picPair.short)) %>% 
    dplyr::select(-tmp)
  return(data.critical)
}


# Bootstrapping -----------------------------------------------------------
get_mean_selection_rate <- function(dat, indices) {
  return(dat[indices, ] %>% pull(response) %>% mean())
}
# returns the statistics calculated for each of N bootstrap samples
get_bootstrap_statistics = function(dat, response, N=1000) {
  bootstrap_samples <- group_map(dat, function(df, df.group) {
    bootstrap = boot(df, statistic = get_mean_selection_rate, R=N)$t
    samples = tibble(rate=bootstrap[,1]) 
    samples.ordered = left_join(df.group, samples, by=character()) %>%
      arrange(rate) %>% rowid_to_column("idx")
    return(samples.ordered)
   }) %>% bind_rows() %>% add_column(response = response)
  return(bootstrap_samples)
}
get_bootstrap_cis = function(bootstrap_samples, N=1000) {
  CI_bounds = c(low = ceiling(0.025 * N), up = ceiling(0.975 * N))
  bootstrapped_cis <- bootstrap_samples %>% 
    group_by(across(c(-idx, -rate, -response))) %>% 
    mutate(CI = case_when(idx == CI_bounds[["low"]] ~ "ci.low",
                          idx == CI_bounds[["up"]] ~ "ci.up")) %>% 
    filter(!is.na(CI)) %>% dplyr::select(-idx) %>% 
    pivot_wider(names_from = "CI", values_from = "rate")
  return(bootstrapped_cis)
}

# returns bootstrapped confidence intervals for each stimulus + qud and
# also only for each qud, with empirical mean rate 
bootstrap_rate_exhaustive = function(df.boot, N=1000){
  df = df.boot %>% group_by(picPair, picPair.long, picPair.short, question)
  samples.picPair_qud = get_bootstrap_statistics(df, "exhaustive", N)
  cis.picPair_qud <- left_join(
    get_bootstrap_cis(samples.picPair_qud, N),
    # empriical rate
    df %>% summarize(rate_exhaustive = mean(response))
  )
  samples.qud = samples.picPair_qud %>% group_by(question) %>% 
    arrange(question, rate) %>% mutate(idx = row_number()) %>% 
    dplyr::select(-picPair.long, -picPair.short, -picPair)
  cis.qud = left_join(
    get_bootstrap_cis(samples.qud, N=samples.qud$idx %>% max()),
    df.boot %>% group_by(question) %>% summarize(rate_exhaustive=mean(response))
  )
  return(list(rate_qud_picPair = cis.picPair_qud, rate_qud = cis.qud))
}

library(here)
source(here("R", "utils.R"))
main_data <- read_csv(here("results", "results_67_communicationBlocks_BG.csv"))
pilot_data <- read_csv(here("results", "results_67_communicationBlocks_BG_pilot.csv"))

# Clean Data --------------------------------------------------------------

# join original pilot data and
# remove data from our test-runs if had not been done before
data.original = bind_rows(pilot_data, main_data) %>% 
  filter(!str_detect(prolific_id, "test-")) %>% 
  mutate(question = case_when(question == "if_ant" ~ "ifp", 
                              question == "cons" ~ "willq",
                              T ~ question))
# clean according to predefined exclusion criteria
data.cleaned = clean_data(data.original)

# sanity check whether there are duplicate prolific ids (according to settings
# in prolific this should normally not be the case)
duplicate_prolific_ids = data.cleaned %>% 
  dplyr::group_by(submission_id, prolific_id, startDate) %>% 
  dplyr::count() %>%
  group_by(prolific_id) %>% dplyr::count() %>% filter(n > 1) %>% 
  pull(prolific_id)

message("prolific_ids that were not unique (despite exclusion criteria defined in Prolific):")
message(paste(duplicate_prolific_ids, collapse=", "))

# remove later data from duplicate prolific_ids:
fst_batch = "Sat Oct 30 2021"
snd_batch = "Mon Nov 15 2021"
pilot_batch = "Thu Sep 23 2021"

duplicates.out = data.original %>% 
  filter(prolific_id %in% duplicate_prolific_ids) %>% 
  dplyr::select(submission_id, prolific_id, startDate) %>% 
  distinct() %>%
  filter(str_detect(startDate, snd_batch))

data = anti_join(data.cleaned, duplicates.out)

# check for comments
comments = data %>% select(comments, submission_id) %>% 
  filter(!is.na(comments) & comments!="") %>% distinct()

comments$comments
comments[11,]$comments
comments[24,]$comments
comments[30,]$comments

out.comments = comments[c(11, 24), ] %>% dplyr::select(submission_id)

data = anti_join(data, out.comments)

n_all = data.original$submission_id %>% unique() %>% length
n_kept = data$submission_id %>% unique() %>% length
ratio_kept = round(n_kept/n_all, 4)
msg <- "% of all recorded data after cleaning that is included in analysis:"
message(paste(msg, ratio_kept))
write_csv(data, here("results", "data_cleaned.csv"))


# numbers removed data ----------------------------------------------------
data.controls <- data.original %>%
  group_by(submission_id) %>% 
  filter(startsWith(type, "control-") | type == "test-example") %>% 
  mutate(correct = expected == selected_pic) %>% 
  dplyr::select(correct, type, submission_id, prolific_id, id, timeSpent,
                comments, startDate) %>% 
  filter(!str_detect(prolific_id, "test-"))

# 1. Wrong answer in test example
out.test_ex = data.controls %>% filter(type == "test-example" & !correct) %>%
  dplyr::select(submission_id, startDate, comments)

# 2. more than 1 control-trial wrong
out.control = data.controls %>% filter(type == "control-physics") %>% 
  group_by(submission_id, prolific_id, timeSpent, comments, startDate) %>% 
  summarize(nb_correct = sum(correct), .groups = "drop_last") %>%
  filter(nb_correct < 2)

# 3. failed attention-checks
attention_checks = data.controls %>% filter(type=="control-random")
out.attention = attention_checks %>%
  group_by(submission_id, prolific_id, timeSpent, comments, startDate) %>% 
  summarize(nb_correct = sum(correct), .groups = "drop_last") %>%
  filter(nb_correct != 3)

# 4. total number of participants removed
out.all = c(out.test_ex$submission_id, out.control$submission_id, 
            out.attention$submission_id, out.comments$submission_id,
            duplicates.out$submission_id) %>% unique() 
out.n = out.all %>% length()

data.keep = anti_join(
  data.original %>% dplyr::select(submission_id, prolific_id) %>% distinct(),
  tibble(submission_id=out.all)
)
# Print numbers
message(paste("#submissions with wrong test-example:", 
              out.test_ex$submission_id %>% unique() %>% length()))
message(paste("#submissions with more than 1 wrong control trial:", 
              out.control$submission_id %>% unique() %>% length()))
message(paste("#submissions with at least one wrong attention-check trial:", 
              out.attention$submission_id %>% unique() %>% length()))
message(paste("#submissions removed due to several recordings:", 
              duplicates.out$submission_id %>% unique() %>% length()))
message(paste("#submissions removed due to comments:", 
              out.comments$submission_id %>% unique() %>% length()))

message(paste("#submissions that were excluded for the analysis:", out.n))


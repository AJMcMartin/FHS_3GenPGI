## Load all packages with single call
load_packs <- function(...) {
  invisible(lapply(substitute(list(...))[-1], function(x) 
    library(deparse(x), character.only = TRUE)))
}

## Setup table of desired phenotypes/exams/cohorts
pheno_setup_fun <- function(pheno_info_table, gen_num, exam_names = NULL, data_acc = NULL, var_names = NULL) {
    if (!is.null(exam_names)) {pheno_info_table <- pheno_info_table %>% filter(EXAM %in% exam_names)}
    if (!is.null(data_acc)) {pheno_info_table <- pheno_info_table %>% filter(data_acc %in% data_acc)}
    if (is.null(var_names)) {var_names <- grep("^VAR.", colnames(pheno_info_table), value = T)}
    
    pheno_info_table <- pheno_info_table %>%
        filter(GEN %in% gen_num) %>% select(GEN, EXAM, DATA.ACC, ID.VAR, all_of(var_names)) %>%
        ## Remove sets for which none of the specified variables are available
        filter_at(vars(all_of(var_names)), any_vars(!is.na(.))) %>%
        rowwise() %>%
        ## NEED TO ADJUST TO ACCOMODATE ALL VARIABLES USING BASE R
        mutate(ALL.VARS1 = paste0(na.omit(c(VAR.HRX, VAR.ALCO)), collapse = "; ")) %>%
        ungroup() %>%
        #unite("ALL.VARS1", all_of(var_names), remove = F, na.rm = T, sep = "; ") %>%
        ## If variables from multiple specified exams are found in the same data set (i.e. in a vr_wkthru set),
        ## combine to avoid loading multiple times
        mutate(EXAM = ifelse(grepl("VR", EXAM), "VR", EXAM)) %>%
        group_by(EXAM, DATA.ACC, ID.VAR) %>%
        summarize(ALL.VARS2 = paste(ALL.VARS1, collapse = "; ")) %>%
        ungroup()
    
    # pheno_list <- pheno_info_table %>%
    #   select(-EXAM) %>%
    #   unite("FINAL.VARS", c("DATA.ACC", "ID.VAR", "ALL.VARS2"), na.rm = T, sep = "; ") %>%
    #   sapply(strsplit, "; ")
    # names(pheno_list) <- pheno_info_table$EXAM
    # return(pheno_list)
    
    pheno_list <- pheno_info_table %>%
        rowwise() %>%
        mutate(FINAL.VARS = paste0(na.omit(c(DATA.ACC, ID.VAR, ALL.VARS2)), collapse = "; ")) %>%
        ungroup() %>%
        select(-c(EXAM, DATA.ACC, ID.VAR, ALL.VARS2)) %>%
        sapply(strsplit, "; ")
    names(pheno_list) <- pheno_info_table$EXAM
    return(pheno_list)
}

## To be used with output of pheno_setup_fun: select and format variables from server
pheno_select_fun <- function(setup_output, full_set = F) {
  base_paths <- paste0("~/dbgap/FCoh", c(1,2), 
                       "/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v29.p10",
                       c(".c1.HMB-IRB-MDS", ".c2.HMB-IRB-NPU-MDS"), "/PhenotypeFiles")
  for (i in 1:length(setup_output)) {
    ## Find and store paths to specified phenotype files
    file_pattern <- paste0(".*", setup_output[[i]][1], "\\..*\\.txt[\\.gz]*$")
    path_vec <- c(grep(file_pattern, list.files(base_paths[1], full.names = T), value = T)[1],
                  grep(file_pattern, list.files(base_paths[2], full.names = T), value = T)[1])
    ## Load phenotype file for each consent group
    c1_dat <- read.table(path_vec[1], header = T, skip = 10, sep= "\t", quote = "", stringsAsFactors = F)
    c2_dat <- read.table(path_vec[2], header = T, skip = 10, sep= "\t", quote = "", stringsAsFactors = F)
    colnames(c1_dat) <- toupper(colnames(c1_dat))
    colnames(c2_dat) <- toupper(colnames(c2_dat))
    ## Combine consent groups; subset to selected variables if provided, else keep full data set
    if (full_set == F) {
      var_vec <- setup_output[[i]][-1]
      c1c2_dat <- rbind(c1_dat[,var_vec], c2_dat[,var_vec])
    } else {
      c1c2_dat <- rbind(c1_dat, c2_dat)
    }
    ## Many variables share names across FFQs
    ## To distinguish/avoid conflict, append exam name shorthand (e.g. E3_FF) to each variable name
    if (grepl("_FF$", names(setup_output)[i])) {
      colnames(c1c2_dat)[-c(1,2)] <- paste(colnames(c1c2_dat)[-c(1,2)], names(setup_output[i]), sep = ".")
    }
    ## pht000107 (resp1_6s) has complex data structure: 
    ## 1-2 observations for each SHAREID (corresponding to exam 5 and exam 6)
    ## Convert long to wide
    if (grepl("_RESP$", names(setup_output)[i])) {
      c1c2_dat <- c1c2_dat %>% 
        pivot_wider(id_cols = c(SHAREID, IDTYPE), names_from = EXAM,
                    values_from = "RQ021", names_prefix = "CPD.RESP") 
    }
    if (i == 1) {
      pheno_out <- c1c2_dat
    } else {
      pheno_out <- merge(pheno_out, c1c2_dat, by = c("SHAREID", "IDTYPE"), all = T)
    }
  }
  return(pheno_out)
}

## Returns list of variable name vectors used to construct composite variable(s)
## Useful for variables like CES-D that have many items
varvec_selector <- function (pheno_table, cohort, exam_vec, var) {
  phenotype <- paste0("VAR.", var)
  list_names <- paste0(cohort, exam_vec, ".", var)
  select_tab <- NULL
  for (i in 1:length(exam_vec)) {
    cohort_exam <- pheno_table %>% filter(GEN == cohort & EXAM == exam_vec[i]) %>% pull(phenotype) %>% str_split(., "; ") %>% unlist()
    select_tab[[i]] <- cohort_exam
  }
  names(select_tab) <- list_names
  return(select_tab)
}
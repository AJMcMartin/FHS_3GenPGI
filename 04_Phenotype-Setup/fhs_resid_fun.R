resid_fun <- function(source_data, pheno) {
  exam_pheno <- str_subset(colnames(source_data), paste0("^", pheno)) ## Full pheno names
  
  #<<==========================================================================>>#
  #<<================ CASE 1: NPB COG Measures ================================>>#
  #<<==========================================================================>>#
  if (pheno == "COG") {
    exam_age <- "AGE_NP.COG"
    exam_age.sq <- "AGE_NP.COG_sq"
    resid_ex_sub <- source_data %>% select(SHAREID, SEX, exam_age, exam_age.sq, exam_pheno)
    resid_ex_out <- umx_residualize(as.formula(paste0(exam_pheno, 
                                                      paste("~ SEX", exam_age, exam_age.sq, paste0("SEX*", exam_age), 
                                                            paste0("SEX*", exam_age.sq), sep = "+"))), data = resid_ex_sub)
    
    ## Standardize Residuals: 1 Observation per Respondent
    resid_all <- resid_ex_out %>% mutate(across(all_of(exam_pheno), ~as.vector(scale(.)))) %>%
      rename_with(., ~str_replace(., "raw", "StdRes"), .cols = exam_pheno) %>%
      select(SHAREID, starts_with(pheno))
    
    #<<========================================================================>>#
    #<<================ CASE 2: Max Value + YOB ===============================>>#
    #<<========================================================================>>#
  } else if (pheno %in% c("EA", "EVSMK", "MDD", "NEB", "AFB", "CPD.LT")) {
    resid_ex_sub <- source_data %>% select(SHAREID, SEX, YOB_V2, YOB_V2_sq, exam_pheno)
    
    if (pheno %in% c("NEB", "AFB")) {
      resid_ex_out <- umx_residualize(as.formula(paste0(
        exam_pheno, paste("~ YOB_V2", "YOB_V2_sq", sep = "+"))), data = resid_ex_sub)
    } else {
      resid_ex_out <- umx_residualize(as.formula(paste0(
        exam_pheno, paste("~ SEX", "YOB_V2", "YOB_V2_sq", "SEX*YOB_V2", "SEX*YOB_V2_sq", sep = "+"))), data = resid_ex_sub)
    }
    ## Standardize Residuals: 1 Observation per Respondent
    resid_all <- resid_ex_out %>% mutate(across(all_of(exam_pheno), ~as.vector(scale(.)))) %>%
      rename_with(., ~str_replace(., "raw", "StdRes"), .cols = exam_pheno) %>%
      select(SHAREID, starts_with(pheno))
    
    #<<========================================================================>>#
    #<<================ CASE 3: Multi-Wave (Exams 1 and 2) + Age ==============>>#
    #<<========================================================================>>#
  } else if (pheno %in% c("CPD.CUR", "HGT", "BMI", "DBP", "SBP", "BG", "HDL", "LDL", "NoHDL")) {
    exam_nums <- str_extract(exam_pheno, paste0("(?<=", pheno, "_E)[:digit:]"))
    exam_age <- paste0("AGE", exam_nums)
    exam_age.sq <- paste0(exam_age, "_sq")
    
    resid_list <- NULL
    for (i in 1:length(exam_pheno)) {
      resid_ex_sub <- source_data %>% select(SHAREID, SEX, exam_age[i], exam_age.sq[i], exam_pheno[i])
      resid_ex_out <- umx_residualize(as.formula(paste0(
        exam_pheno[i], paste("~ SEX", exam_age[i], exam_age.sq[i], paste0("SEX*", exam_age[i]), 
                             paste0("SEX*", exam_age.sq[i]), sep = "+"))), data = resid_ex_sub)
      resid_list[[i]] <- resid_ex_out
    }
    
    resid_all <- full_join(resid_list[[1]], resid_list[[2]], by = c("SHAREID", "SEX")) %>% 
      ## Standardize residuals for each exam
      mutate(across(all_of(exam_pheno), ~as.vector(scale(.)))) %>%
      rowwise() %>%
      ## Obtain mean of standardized residuals
      mutate(Avg.StdRes = na_if(mean(c_across(all_of(exam_pheno)), na.rm = T), NaN)) %>%
      ungroup() %>%
      rename_with(., ~str_replace(., "raw", "StdRes"), .cols = exam_pheno) %>%
      rename_with(., ~paste0(pheno, "_", .), .cols = "Avg.StdRes") %>%
      ## Select Avg variable only for now
      select(SHAREID, contains("Avg.StdRes"))
    
    #<<========================================================================>>#
    #<<================ CASE 4: Multi-Wave (Non-Standard Exam Pair) + Age =====>>#
    #<<========================================================================>>#
  } else {
    ## For DEP, FBG, FEV1
    ## Only difference from CASE 3 is phenotype name included in each AGE variable
    exam_nums <- str_extract(exam_pheno, paste0("(?<=", pheno, "_E)[:digit:]"))
    exam_age <- paste0("AGE", exam_nums, ".", pheno)
    exam_age.sq <- paste0(exam_age, "_sq")
    
    resid_list <- NULL
    for (i in 1:length(exam_pheno)) {
      resid_ex_sub <- source_data %>% select(SHAREID, SEX, exam_age[i], exam_age.sq[i], exam_pheno[i])
      resid_ex_out <- umx_residualize(as.formula(paste0(
        exam_pheno[i], paste("~ SEX", exam_age[i], exam_age.sq[i], paste0("SEX*", exam_age[i]), 
                             paste0("SEX*", exam_age.sq[i]), sep = "+"))), data = resid_ex_sub)
      resid_list[[i]] <- resid_ex_out
    }
    
    resid_all <- full_join(resid_list[[1]], resid_list[[2]], by = c("SHAREID", "SEX")) %>% 
      ## Standardize residuals for each exam
      mutate(across(all_of(exam_pheno), ~as.vector(scale(.)))) %>%
      rowwise() %>%
      ## Obtain mean of standardized residuals
      mutate(Avg.StdRes = na_if(mean(c_across(all_of(exam_pheno)), na.rm = T), NaN)) %>%
      ungroup() %>%
      rename_with(., ~str_replace(., "raw", "StdRes"), .cols = exam_pheno) %>%
      rename_with(., ~paste0(pheno, "_", .), .cols = "Avg.StdRes") %>%
      ## Select Avg variable only for now
      select(SHAREID, contains("Avg.StdRes"))
  }
  return(resid_all)
}









#' Simulated data from AntiBiotics for Children with severe Diarrhea (ABCD) trial
#' 
#' Simulated dataset based closely on real data from ABCD trial with 40 covariates 
#' and n = 6692 observations. The outcome variable is length-for-age z-score at 
#' day 90 of the trial (laz90) and the treatment variable is azithromycin (an_grp_01).
#' Details on the distribution of each variable can be found in *cite paper?*
#' 
#' @docType data
#' 
#' @usage data(abcd_data)
#' 
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{pid}{Participant identification number}
#'  \item{an_grp_01}{Binary indicator for azithromycin (antibiotic) treatment, 
#'  1=assigned azithromycin, 0=assigned placebo}
#'  \item{rotavirus_new}{Relative quantity of rotavirus in stool sample based 
#'  on log-10 transformed qPCR cycle threshold values}
#'  \item{rotavirus_bin}{Binary indicator for rotavirus presence in stool sample 
#'  (1=rotavirus detected, 0=rotavirus not detected)}
#'  \item{norovirus_new}{Relative quantity of norovirus in stool sample based on
#'   log-10 transformed qPCR cycle threshold values}
#'  \item{norovirus_bin}{Binary indicator for norovirus presence in stool sample
#'   (1=norovirus detected, 0=norovirus not detected)}
#'  \item{adenovirus_new}{Relative quantity of adenovirus in stool sample based 
#'  on log-10 transformed qPCR cycle threshold values}
#'  \item{adenovirus_bin}{Binary indicator for adenovirus presence in stool sample
#'   (1=adenovirus detected, 0=adenovirus not detected)}
#'  \item{astrovirus_new}{Relative quantity of astrovirus in stool sample based 
#'  on log-10 transformed qPCR cycle threshold values}
#'  \item{astrovirus_bin}{Binary indicator for astrovirus presence in stool sample
#'   (1=astrovirus detected, 0=astrovirus not detected)}
#'  \item{sapovirus_new}{Relative quantity of sapovirus in stool sample based on 
#'  log-10 transformed qPCR cycle threshold values}
#'  \item{sapovirus_bin}{Binary indicator for sapovirus presence in stool sample
#'   (1=sapovirus detected, 0=sapovirus not detected)}
#'  \item{st_etec_new}{Relative quantity of st_etec in stool sample based on 
#'  log-10 transformed qPCR cycle threshold values}
#'  \item{st_etec_bin}{Binary indicator for st_etec presence in stool sample 
#'  (1=st_etec detected, 0=st_etec not detected)}
#'  \item{shigella_new}{Relative quantity of shigella in stool sample based on 
#'  log-10 transformed qPCR cycle threshold values}
#'  \item{shigella_bin}{Binary indicator for shigella presence in stool sample 
#'  (1=shigella detected, 0=shigella not detected)}
#'  \item{campylobacter_new}{Relative quantity of campylobacter in stool sample 
#'  based on log-10 transformed qPCR cycle threshold values}
#'  \item{campylobacter_bin}{Binary indicator for campylobacter presence in stool 
#'  sample (1=campylobacter detected, 0=campylobacter not detected)}
#'  \item{tepec_new}{Relative quantity of tepec in stool sample based on log-10
#'   transformed qPCR cycle threshold values}
#'  \item{tepec_bin}{Binary indicator for tepec presence in stool sample 
#'  (1=tepec detected, 0=tepec not detected)}
#'  \item{v_cholerae_new}{Relative quantity of v_cholerae in stool sample based 
#'  on log-10 transformed qPCR cycle threshold values}
#'  \item{v_cholerae_bin}{Binary indicator for v_cholerae presence in stool sample
#'   (1=v_cholerae detected, 0=v_cholerae not detected)}
#'  \item{salmonella_new}{Relative quantity of salmonella in stool sample based 
#'  on log-10 transformed qPCR cycle threshold values}
#'  \item{salmonella_bin}{Binary indicator for salmonella presence in stool sample 
#'  (1=salmonella detected, 0=salmonella not detected)}
#'  \item{cryptosporidium_new}{Relative quantity of cryptosporidium in stool sample
#'   based on log-10 transformed qPCR cycle threshold values}
#'  \item{cryptosporidium_bin}{Binary indicator for cryptosporidium presence in 
#'  stool sample (1=cryptosporidium detected, 0=cryptosporidium not detected)}
#'  \item{dy1_scrn_vomitall}{Factor with two levels for vomiting at screening 
#'  ("No" = not vomit, "Yes" = vomit)}
#'  \item{dy1_scrn_lstools}{Number of loose stools in 24 hours prior to enrollment
#'   (continuous)}
#'  \item{dy1_scrn_sstools}{Number of solid stools in 24 hours prior to enrollment 
#'  (continuous)}
#'  \item{dy1_scrn_diardays}{Duration (days) of diarrhea illness prior to enrollment
#'   (continuous)}
#'  \item{dy1_scrn_dehydr}{Dehydration status at screening 
#'  (ordinal, 1 = "No dehydration", 2 = "Some dehydration", 3 = "Severe dehydration")}
#'  \item{avemuac}{Middle upper arm circumference (continuous)}
#'  \item{wfazscore}{Weight for age z-score (continuous)}
#'  \item{lfazscore}{Length for age z-score (continuous)}
#'  \item{wflzscore}{Weight for length z-score (continuous)}
#'  \item{site}{Study site (nominal, 2 = Bangladesh, 3 = Kenya, 4 = Malawi,
#'   5 = Mali, 6 = India, 7 = Tanzania, 8 = Pakistan)}
#'  \item{dy1_ant_sex}{Sex (1 = Male, 2 = Female)}
#'  \item{agemchild}{Age (in months) (continuous)}
#'  \item{an_ses_quintile}{SES quintile (ordinal, 1 = 1st quintile, 
#'  2 = 2nd quintile, 3 = 3rd quintile, 4 = 4th quintile, 5 = 5th quintile)}
#'  \item{an_tothhlt5}{Number of children under age 5 years in household (continuous)}
#'  \item{rotaseason}{Binary indicator for enrollment during rotavirus season 
#'  (1 = enrolled during rotavirus season, 0 = not enrolled during rotavirus season)}
#'  \item{month_en}{Month enrolled (1-12 correspond to January-December)}
#'  \item{lazd90}{Length for age z-score at day 90 post-enrollment}
#'  }
#' 
#' @keywords datasets
#' 
#' @references cite the paper?
#' 
#' @examples
#' data(abcd_data)
#' head(abcd_data)
#' hist(abcd_data$lazd90)
"abcd_data"
if(!require(stringr)){
  install.packages("stringr")
}
if(!require(lubridate)){
  install.packages("lubridate")
}
library(stringr)
library(lubridate)

########################################
## Setting Area
# Project informations

project_id <- "TBD171053"
reporting_date <- as.Date("2018-11-01")
data_path <- paste("/Volumes/TBI_siyoo/TBI_NonHumanTeam/Report-repository", "TBD171053-KRIBB-Monkey-MethylSeq-TestReport-sabaeus-20181031", sep = '/')

# Analysis informations
species_id <- "C_sabaeus_ENS93"
species_name <- "Chlorocebus sabaeus"
species_alias <- "ChlSab1.1"
assembly_acc <- "GCA_000409795.2"
geneset_ver <- "Ensembl 93"

# mine_report informations
mine_report_path <- "/Volumes/TBI_siyoo/TBI_Development/MINE/mine_report"


########################################

xmonth <- month(as.Date(reporting_date))
xday <- day(as.Date(reporting_date))
dmonth <- stringr::str_pad(xmonth,2,pad="0")
dday <- stringr::str_pad(xday,2,pad="0")

try(
  rmarkdown::render(paste(mine_report_path, "report_template/mine_report_typeA.Rmd", sep = '/'),
                    params = list(
                      dmonth = dmonth,
                      dday = dday,
                      project_id = project_id,
                      data_path = data_path,
                      species_id = species_id,
                      species_name = species_name,
                      species_alias = species_alias,
                      assembly_acc = assembly_acc,
                      geneset_ver = geneset_ver,
                      mine_report_path = mine_report_path
                      ),
                    output_format = "pdf_document",
                    output_file = paste(data_path, paste0(project_id, "_bisulfite-seq_analysis_report", ".pdf"), sep = '/'),
                    encoding = "UTF-8"),
  silent = TRUE
)

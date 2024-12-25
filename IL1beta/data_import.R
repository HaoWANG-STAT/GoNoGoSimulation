library(tidyverse)

d <- read.csv("interimSampleSize/IL1beta/UBIOPREDjohn.csv") %>% 
  select(-X, -TRIALNAME, -PROBE, -Time.Point, -Sample.Code, -Assay.ID, -Platform.ID, -Sample.Type, -Tissue.Type, -GENE.ID, -GENE.SYMBOL) %>% 
  rename_with(
    function(x) {
      stringr::str_replace_all(
        x,
        c(
          "Demographic\\.Data\\." = "",
          "Subject.ID" = "SubjectID",
          "Study\\.Groups\\." = "",
          "Biomarker\\.Data\\.Baseline\\.Visit\\." = "",
          "\\.Master\\." = "",
          "BL\\." = "",
          "Serum\\." = "",
          "SputumPCT." = "",
          "ACQ5\\.Total\\.\\.Imputed\\.\\." = "ACQ5Imputed",
          "ACQ5\\.Total\\.\\.Raw\\.\\." = "ACQ5Raw",
          "Average\\.\\.ACQ1\\.ACQ5\\.\\." = "ACQ5Average",
          "FEV1.PCT" = "FEV1",
          "hsCRP.mg.L" = "hsCRP",
          "\\." = "",
          "VALUE" = "IL1beta",
          "LOG2E" = "LogIL1beta",
          "ZSCORE" = "ZIL1beta"
        )
      )
    }
  ) %>% 
  as_tibble()

saveRDS(d, "interimSampleSize/IL1beta/UBIOPRED.Rds")

# import CRP data provided by Richard

library(here)

d2 <- readxl::read_xlsx(here("IL1beta", "Lute_verse_CRP.xlsx"))

saveRDS(d2, here("IL1beta", "LuteVerse.Rds"))

d3 <- readxl::read_xlsx(here("IL1beta", "Lute_Verse_CRP_EOS_NEUT.xlsx"))

saveRDS(d3, here("IL1beta", "LuteVerse_new.Rds"))

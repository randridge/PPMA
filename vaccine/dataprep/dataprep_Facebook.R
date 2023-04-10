library(tidyverse)
library(lubridate)

fb1 <- read_csv("./data_Facebook/2021-01-race-ethnicity.csv", guess_max=500000)
fb2 <- read_csv("./data_Facebook/2021-02-race-ethnicity.csv")
fb3 <- read_csv("./data_Facebook/2021-03-race-ethnicity.csv")
fb4 <- read_csv("./data_Facebook/2021-04-race-ethnicity.csv")
fb5 <- read_csv("./data_Facebook/2021-05-race-ethnicity.csv")

fb <- bind_rows(fb1, fb2, fb3, fb4, fb5) %>% 
  dplyr::select(StartDatetime,EndDatetime,A3,V1,V3,D1,D2,D8,raceethnicity,wave,fips,weight) %>%
  mutate(date=date(StartDatetime), 
         week=epiweek(StartDatetime), 
         state=substr(fips,1,2),
         vaccinated=(V1==1),
         vacstatus=ifelse(V3 %in% c(1,2), 2, 
                          ifelse(V3 %in% c(3,4), 3, 
                                  ifelse(V1==1, 1, NA))),
         male=(D1==1),
         male_nomiss=ifelse(D1==5, NA, male),
         educ_lthsgrad=(D8==1),
         educ_hs=(D8==2),
         educ_somecoll=(D8 ==3|D8==4),
         educ_bach=(D8==5),
         educ_graddeg=(D8==8|D8==6|D8==7),
         raceth_whitenh=(raceethnicity=="NonHispanicWhite"),
         raceth_blacknh=(raceethnicity=="NonHispanicBlackAfricanAmerican"),
         raceth_aiannh=(raceethnicity=="NonHispanicAmericanIndianAlaskaNative"),
         raceth_asianaapinh=(raceethnicity=="NonHispanicAsian"|raceethnicity=="NonHispanicNativeHawaiianPacificIslander"),
         raceth_othernh=(raceethnicity=="NonHispanicMultipleOther"),
         raceth_hisp=(raceethnicity=="Hispanic"),
         age_18_24=(D2==1),
         age_25_34=(D2==2),
         age_35_44=(D2==3),
         age_45_54=(D2==4),
         age_55_64=(D2==5),
         age_65_74=(D2==6),
         age_75plus=(D2==7)) %>%
  filter(week >= 2 & week <= 18) %>%
  select(wave,weight,date:age_75plus)

write_csv(fb, file="./data/microdata_facebook_week2to18.csv")

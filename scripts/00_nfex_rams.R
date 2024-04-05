#This script reads the data from my google drive (RAW DATA I got from the QE), creates RT bounds to be able to cut
#the peaks we want then you look at the peaks and get areas which can then be used for the Targeted Pipeline


#Load libraries
library(stringr)
library(RaMS)
library(plotly)
library(data.table)
library(tidyverse)
library(janitor)


nfexdata <- grabMSdata(list.files("/Users/sgarcia/Library/CloudStorage/GoogleDrive-sggarcia@uw.edu/My Drive/NFEX", full.names = TRUE))

mass_calcs <- read_csv("intermediates/RAMS_mass_calcs.csv")

#Setting the bounds for each chromotograph
rtbounds <- tribble(
  ~compound_parent_name, ~rtstart, ~rtend,
  "Choline", 14.6, 15.6,
  "Carnitine", 9.7, 10.5,
  "Proline", 9, 10,
  "Proline Betaine", 6.7, 7.5,
  "Lysine", 17, 19.5,
  "Trimethylamine N-Oxide", 6.9, 7.7,
  "Glycine", 11.6, 12.4,
  "Dimethylglycine", 8, 9,
  "Ectoine", 9.5, 10.5)

#Data wrangling and setting up the df I want to plot
peak_bounds_complete <- read_csv("intermediates/RAMS_mass_calcs.csv") %>%
  clean_names() %>%
  mutate( mz = value) %>%
  filter(!is.na(mz)) %>%
  select ("compound_name","mz") %>%
  #dplyr::rename(compound_name=`compound name`) %>%
  mutate(N=str_extract(compound_name, "(?<= )\\d+(?=\\(15N\\))")) %>%
  mutate(C=str_extract(compound_name, "(?<= )\\d+(?=\\(13C\\))")) %>%
  mutate(N = as.numeric(N),
         C = as.numeric(C)) %>%
  mutate(compound_parent_name=str_remove(compound_name, " \\d.*")) %>%
  left_join(rtbounds)

#Use PMAP function to go through row by row by chopping by the ppm, then the intensity and then the retention time too
all_unlabeled_eics <- peak_bounds_complete %>%
  filter(N == 0, C == 0) %>%
  pmap(function(...){
    row_data <- list(...)
    nfexdata$MS1[mz%between%pmppm(row_data$mz, 10)] %>%
      filter(rt%between%c(row_data$rtstart, row_data$rtend)) %>%
      mutate(compound_name=row_data$compound_name,
             N=row_data$N, C=row_data$C) #this part makes sure we have all the data and now just the first row
  }) %>%
#This part focuses on data wrangling again
  bind_rows() %>%
  mutate(treatment=str_extract(filename, "NH4|Urea|NO3")) %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(timepoint_num=as.numeric(str_remove(timepoint, "T"))) %>%
  arrange(timepoint_num) %>%
  mutate(timepoint=fct_inorder(timepoint))

#Semilabled for just 1 N and 0 C
all_N1_C0_labeled_eics <- peak_bounds_complete %>%
  filter(N == 1, C == 0) %>%
  pmap(function(...){
    row_data <- list(...)
    nfexdata$MS1[mz%between%pmppm(row_data$mz, 10)] %>%
      filter(rt%between%c(row_data$rtstart, row_data$rtend)) %>%
      mutate(compound_name=row_data$compound_name,
             N=row_data$N, C=row_data$C) #this part makes sure we have all the data and now just the first row
  }) %>%
  #This part focuses on data wrangling again
  bind_rows() %>%
  mutate(treatment=str_extract(filename, "NH4|Urea|NO3")) %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(timepoint_num=as.numeric(str_remove(timepoint, "T"))) %>%
  arrange(timepoint_num) %>%
  mutate(timepoint=fct_inorder(timepoint))

#Semilabled for just 1 N and 0 C
all_N1_C1_labeled_eics <- peak_bounds_complete %>%
  filter(N == 1, C == 1) %>%
  pmap(function(...){
    row_data <- list(...)
    nfexdata$MS1[mz%between%pmppm(row_data$mz, 10)] %>%
      filter(rt%between%c(row_data$rtstart, row_data$rtend)) %>%
      mutate(compound_name=row_data$compound_name,
             N=row_data$N, C=row_data$C) #this part makes sure we have all the data and now just the first row
  }) %>%
  #This part focuses on data wrangling again
  bind_rows() %>%
  mutate(treatment=str_extract(filename, "NH4|Urea|NO3")) %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(timepoint_num=as.numeric(str_remove(timepoint, "T"))) %>%
  arrange(timepoint_num) %>%
  mutate(timepoint=fct_inorder(timepoint))

#Semilabled for just 2 N and 0 C
all_N2_C0_labeled_eics <- peak_bounds_complete %>%
  filter(N == 2, C == 0) %>%
  pmap(function(...){
    row_data <- list(...)
    nfexdata$MS1[mz%between%pmppm(row_data$mz, 10)] %>%
      filter(rt%between%c(row_data$rtstart, row_data$rtend)) %>%
      mutate(compound_name=row_data$compound_name,
             N=row_data$N, C=row_data$C) #this part makes sure we have all the data and now just the first row
  }) %>%
  #This part focuses on data wrangling again
  bind_rows() %>%
  mutate(treatment=str_extract(filename, "NH4|Urea|NO3")) %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(timepoint_num=as.numeric(str_remove(timepoint, "T"))) %>%
  arrange(timepoint_num) %>%
  mutate(timepoint=fct_inorder(timepoint))

#Semilabled for just 2 N and 1 C *this doesnt work idk why
all_N2_C1_labeled_eics <- peak_bounds_complete %>%
  filter(N == 2, C == 1) %>%
  pmap(function(...){
    row_data <- list(...)
    nfexdata$MS1[mz%between%pmppm(row_data$mz, 10)] %>%
      filter(rt%between%c(row_data$rtstart, row_data$rtend)) %>%
      mutate(compound_name=row_data$compound_name,
             N=row_data$N, C=row_data$C) #this part makes sure we have all the data and now just the first row
  }) %>%
  #This part focuses on data wrangling again
  bind_rows() %>%
  mutate(treatment=str_extract(filename, "NH4|Urea|NO3")) %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(timepoint_num=as.numeric(str_remove(timepoint, "T"))) %>%
  arrange(timepoint_num) %>%
  mutate(timepoint=fct_inorder(timepoint))



#Semilabled for just 0 N and 1 C
all_N0_C1_labeled_eics <- peak_bounds_complete %>%
  filter(N == 0, C == 1) %>%
  pmap(function(...){
    row_data <- list(...)
    nfexdata$MS1[mz%between%pmppm(row_data$mz, 10)] %>%
      filter(rt%between%c(row_data$rtstart, row_data$rtend)) %>%
      mutate(compound_name=row_data$compound_name,
             N=row_data$N, C=row_data$C) #this part makes sure we have all the data and now just the first row
  }) %>%
  #This part focuses on data wrangling again
  bind_rows() %>%
  mutate(treatment=str_extract(filename, "NH4|Urea|NO3")) %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(timepoint_num=as.numeric(str_remove(timepoint, "T"))) %>%
  arrange(timepoint_num) %>%
  mutate(timepoint=fct_inorder(timepoint))


#Binding all the previous sheets together!

all_eics <- rbind(all_unlabeled_eics, all_N0_C1_labeled_eics, all_N1_C0_labeled_eics, all_N1_C1_labeled_eics,
                 all_N2_C0_labeled_eics, all_N2_C1_labeled_eics)


# Plotting them all together
ggplot(all_eics) +
  geom_line(aes(x=rt, y=int, group=filename)) +
  facet_wrap(~compound_name, ncol=1, scales = "free_y")

#Looking at each compound to truly find their retention time bounds
single_eics <- all_eics %>%
  mutate(samp_type=str_extract(filename, "Poo|Blk|Smp|Std")) %>%
  filter (compound_name == "Dimethylglycine 2(15N), 1(13C)",
          samp_type == "Std")
 #Plots them and makes them interactive to be able to see which Mix
gp <- ggplot(single_eics) +
  geom_line(aes(x=rt, y=int, group=filename)) +
  facet_wrap(~compound_name, ncol=1, scales = "free_y")
ggplotly(gp)

wip <- all_eics %>%
  group_by(filename, compound_name) %>%
  summarise(area=trapz(rt, int)) #%>%
  #pivot_wider(names_from = compound_name, values_from = area)


write_csv(wip, "data_raw/NFEX2_parti_rawdata.csv")


#Rams example code by W. kumler ----

msdata_dir <- system.file("extdata", package = "RaMS")
msdata_files <- list.files(msdata_dir, pattern = "mzML", full.names=TRUE)

msdata <- grabMSdata(files = msdata_files[2:4], grab_what = c("BPC", "MS1"))

knitr::kable(head(msdata$BPC, 3))

plot(msdata$BPC$rt, msdata$BPC$int, type = "l", ylab="Intensity")

ggplot(msdata$BPC) + geom_line(aes(x = rt, y=int, color=filename)) +
  facet_wrap(~filename, scales = "free_y", ncol = 1) +
  labs(x="Retention time (min)", y="Intensity", color="File name: ") +
  theme(legend.position="top")

knitr::kable(head(msdata$MS1, 3))

M <- 118.0865
M_13C <- M + 1.003355
M_15N <- M + 0.997035

iso_data <- imap_dfr(lst(M, M_13C, M_15N), function(mass, isotope){
  peak_data <- msdata$MS1[mz%between%pmppm(mass) & rt%between%c(7.6, 8.2)]
  cbind(peak_data, isotope)
})

iso_data %>%
  group_by(filename, isotope) %>%
  summarise(area=trapz(rt, int)) %>%
  pivot_wider(names_from = isotope, values_from = area) %>%
  mutate(ratio_13C_12C = M_13C/M) %>%
  mutate(ratio_15N_14N = M_15N/M) %>%
  select(filename, contains("ratio")) %>%
  pivot_longer(cols = contains("ratio"), names_to = "isotope") %>%
  group_by(isotope) %>%
  summarize(avg_ratio = mean(value), sd_ratio = sd(value), .groups="drop") %>%
  mutate(isotope=str_extract(isotope, "(?<=_).*(?=_)")) %>%
  knitr::kable()














# Will Code ----
all_eics %>%
  filter(treatment=="NO3") %>%
  #filter(compound_name!="Homarine") %>%
  group_by(filename, compound_name, N, C, timepoint, timepoint_num) %>%
  summarize(area=trapz(rt, int)) %>%
  ungroup() %>%
  filter(C==0) %>%
  ggplot() +
  geom_col(aes(x=timepoint, y=area, fill=N), position="fill", color="black") +
  facet_wrap(~compound_name)


all_eics %>%
  filter(treatment=="NH4") %>%
  filter(C<=2) %>%
  filter(compound_name=="Guanine") %>%
  mutate(N=paste0(N, "N")) %>%
  mutate(C=paste0(C, "C")) %>%
  filter(rt%between%c(7, 10)) %>%
  qplotMS1data() +
  ggh4x::facet_grid2(N~C, scales="free_y", independent = "y")


all_eics %>%
  filter(treatment=="NH4") %>%
  filter(C<=2) %>%
  filter(compound_name=="Guanine") %>%
  mutate(N=paste0(N, "N")) %>%
  mutate(C=paste0(C, "C")) %>%
  filter(rt%between%c(7, 10)) %>%
  ggplot() +
  geom_point(aes(x=rt, y=mz, color=log10(int))) +
  ggh4x::facet_grid2(N~C, scales="free_y", independent = "y")

all_eics %>%
  filter(C==1) %>%
  filter(N==2) %>%
  filter(compound_name=="Guanine") %>%
  filter(treatment=="NH4") %>%
  qplotMS1data() +
  facet_wrap(~filename)






nfexdata$MS1[rt%between%c(7.6, 8.1)] %>%
  .[filename=="230810_Smp_G5_NFEX_T24_NH4_Metab_C.mzML"] %>%
  .[mz%between%pmppm(155.0547, 100000)] %>%
  ggplot() +
  geom_point(aes(x=rt, y=mz, color=log10(int))) +
  scale_color_viridis_c()


nfexdata$MS1[rt%between%c(7.6, 8.1)] %>%
  .[filename=="230810_Smp_G5_NFEX_T24_NH4_Metab_C.mzML"] %>%
  .[mz%between%pmppm(155.0547, 10000)] %>%
  arrange(desc(int)) %>%
  mutate(mz_group=mz_group(mz, 10))  %>%
  ggplot() +
  geom_point(aes(x=rt, y=mz, color=factor(mz_group)))


grp_msdata <- nfexdata$MS1[rt%between%c(7.6, 8.1)] %>%
  .[filename=="230810_Smp_G5_NFEX_T24_NH4_Metab_C.mzML"] %>%
  arrange(desc(int)) %>%
  mutate(mz_group=mz_group(mz, 10)) %>%
  group_by(mz_group) %>%
  filter(n()>10)
grp5_cors <- grp_msdata %>%
  select(-filename, -mz) %>%
  arrange(rt) %>%
  pivot_wider(names_from=rt, values_from=int, values_fn=max) %>%
  mutate(mz_group=paste("Grp", mz_group)) %>%
  column_to_rownames("mz_group") %>%
  data.matrix() %>%
  t() %>%
  cor(use = "pairwise") %>%
  .[,"Grp 5"]
sus_groups <- as.numeric(str_remove(names(which(grp5_cors>0.9)), "Grp "))

grp_msdata %>%
  filter(mz_group%in%sus_groups) %>%
  # qplotMS1data() + facet_wrap(~mz_group, ncol=2, scales="free_y")
  group_by(mz_group) %>%
  summarize(med_mz=median(mz)) %>%
  arrange(med_mz)

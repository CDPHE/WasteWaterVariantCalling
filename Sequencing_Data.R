#Sequencing Data compiler----

SequencingDataPath <- list.files(path = "C:/Users/casalata/Documents/R Codes/imported data/CDPHE Lab/Sequencing Data/") %>%
  as.data.frame() %>%
  `colnames<-`("filename") %>%
  filter(grepl("WWT_mutation_frequencies ",filename,fixed = TRUE)) %>%
  dplyr::mutate(date = str_remove(filename, "WWT_mutation_frequencies ") %>%
                  str_remove(" .xlsx") %>%
                  mdy()) %>%
  filter(date == max(date)) %>%
  dplyr::select(filename) %>%
  unlist() %>%
  paste0("C:/Users/casalata/Documents/R Codes/imported data/CDPHE Lab/Sequencing Data/",.)

#see https://readxl.tidyverse.org/articles/articles/multiple-header-rows.html

#Get the names of the mutations and the associated change
mutationChange <- read_excel(path = SequencingDataPath, n_max = 1, skip =0)%>%
  dplyr::select(-AA_change) %>%
  pivot_longer(everything()) %>%
  `colnames<-`(c("Mutation","change"))

#Extract the variant names
LineageNames <- read_excel(path = SequencingDataPath, n_max = 1, skip =2) %>%
  dplyr::select(-(1:2)) %>%
  colnames()

UniquemutationLineages <- LineageNames %>%
  as.data.frame() %>%
  `colnames<-`("AllVariants") %>%
  filter(AllVariants != "Lineages") %>%
  rowwise() %>%
  mutate(AllVariants = case_when(
    grepl(pattern = "...",x=AllVariants,fixed=TRUE) ~ gsub("\\.\\.\\..*","",AllVariants),
    TRUE ~ AllVariants)) %>%
  mutate(listVariants = str_split(AllVariants, pattern = ";",simplify=FALSE)) %>%
  bind_cols(mutationChange)

AllVariants <- unique(UniquemutationLineages$listVariants) %>%
  unlist() %>%
  unique()

#Get the sequencing data
SequencingData <- read_excel(path = SequencingDataPath, skip = 2)

#Make the sequencing data long and add in the lineage names. May need to update line 48 with new column lengths added to C:\Users\casalata\Documents\R Codes\imported data\CDPHE Lab\sequencing data
SequencingDataLong <- SequencingData %>%
  dplyr::select(-3,-(starts_with("..."))) %>% 
  `colnames<-`(c("Utility","date",mutationChange$Mutation)) %>%
  dplyr::mutate(date = ymd(date)) %>%
  dplyr::mutate(across(3:ncol(.),.fns = as.numeric)) %>%
  pivot_longer(-(Utility:date),
               names_to = "Mutation",
               values_to = "Frequency") %>%
  left_join(UniquemutationLineages,by="Mutation")

ShortSystemData <- read_excel("C:/Users/casalata/Documents/R Codes/reference data/System Data.xlsx")%>%
  mutate(CSULong = str_replace_all(CSULong, "[[:space:]]", "")) %>%
  dplyr::select(wwtp_name,`CSU Internal Identifier`)



# Sequencing Calls 2 ------------------------------------------------------
#Create a new blank data frame
SequencingCalls2 <- data.frame()

#Set the thresholds
detectionThreshold <- 0.05

callThreshold <- 0.3 

VariantsToDisplay <- c("B.1.1.7","P.1",
                       "B.1.351","B.1.617.2",
                       #"A.27","AY.1",
                       "B.1.621",
                       "B.1.1.529", "BA.1", "BA.2")

#Loop through each variant
for(i in 1:length(AllVariants)){
  #Choose the variant you care about for this loop
  VariantOfInterest <- AllVariants[i]
  
  #Identify the mutations associated with that variant
  MutationsOfInterest <- UniquemutationLineages %>%
    rowwise() %>%
    filter(VariantOfInterest %in% listVariants) %>%
    dplyr::select(Mutation) %>%
    unlist()
  

  #This generates the output data for each variant
  Output <- SequencingDataLong %>%
    
    #Look only at the mutations you care about
    filter(Mutation %in% MutationsOfInterest) %>%
    
    #Consider each utility / date combo separately
    group_by(Utility,date) %>%
    
    #Check if each mutaton is over the detection threshold (Defined above)
    dplyr::mutate(Detection = (Frequency > detectionThreshold)) %>%
    
    #Then create your output. First, write down the variant
    dplyr::summarize(Variant = VariantOfInterest,
                     #Record how many associated mutations were detected
                     Detections = sum(Detection,na.rm=TRUE),
                     #Record how many mutatons you analyzed
                     MutationsAnalyzed = n()) %>%
    
    #Calculate the proportion of mutations detected
    dplyr::mutate(PercentDetections = Detections/MutationsAnalyzed)
  
  #Bind these results to the growing dataframe
  SequencingCalls2 <- bind_rows(SequencingCalls2,Output)
}
#Make a graph
SequencingCalls2 %>%
  #Look at only variants where enough mutatons were detected (i.e. the percent
  #of mutations detected exceeds the "call threshold")
  filter(PercentDetections > callThreshold) %>%
  
  #Change the name so this will bind with the system names
  dplyr::rename(`CSU Internal Identifier` = Utility) %>% 
  
  #Bind in the plain-language system names
  left_join(ShortSystemData,by="CSU Internal Identifier") %>%
  
  #Choose only the variants you care about today
  filter(Variant %in% VariantsToDisplay) %>%
  
  #Remove utilities without a valid name
  filter(!is.na(wwtp_name)) %>%
  
  #Sort the utilities and make them a factor so they display nice
  arrange(desc(wwtp_name)) %>%
  dplyr::mutate(wwtp_name = factor(wwtp_name,levels=unique(wwtp_name))) %>%
  
  #Initialize the graph
  ggplot(aes(x=date,y=wwtp_name)) +
  
  #Make a dot plot, with some horizontal jitter to get rid of overlaps
  geom_jitter(aes(color=Variant),width=1,height=0,size=3,alpha=0.5) +
  
  #Make the theme classy
  theme_fivethirtyeight()+
  
  #Add a cool title
  labs(title="Variant Mutation Detections by Utility")


#Generate the output file
SequencingOutput2 <- SequencingCalls2 %>%
  #Choose only the variable syou want  
  transmute(`CSU Internal Identifier` = Utility,
            date = date,
            Variant = Variant,
            #Change this to true/false based on the call threshold
            Present = (PercentDetections > callThreshold)) %>%
  
  #Add in system data
  left_join(ShortSystemData,by="CSU Internal Identifier") %>%
  
  #Choose only the variables
  filter(Variant %in% VariantsToDisplay)
# # Replace values in cells
SequencingOutput2$Variant[SequencingOutput2$Variant== "B.1.1.7"] <- "Alpha (B.1.1.7)"
SequencingOutput2$Variant[SequencingOutput2$Variant=="P.1"] <- "Gamma (P.1)"
SequencingOutput2$Variant[SequencingOutput2$Variant=="B.1.351"] <- "Beta (B.1.351)"
SequencingOutput2$Variant[SequencingOutput2$Variant=="B.1.617.2"] <- "Delta (B.1.617.2)"
SequencingOutput2$Variant[SequencingOutput2$Variant=="B.1.1.529"] <- "Omicron (B.1.1.529)"


#Alright, if that looks good it's time to export it
OutputName = paste0("output/SequencingData2",Sys.Date(),".csv")
write_csv(SequencingOutput2,OutputName,na = "")

OutputName = ("InternalDashboard/WWInternalDashboard/data/SequencingData2.csv")
write_csv(SequencingOutput2,OutputName,na = "")



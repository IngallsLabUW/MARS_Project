## 


raw_stans <- read.csv(paste0("https://raw.githubusercontent.com/",
                             "IngallsLabUW/Ingalls_Standards/",
                             "b098927ea0089b6e7a31e1758e7c7eaad5408535/",
                             "Ingalls_Lab_Standards_NEW.csv")) %>%
  filter(Column=="HILIC") %>%
  mutate(polarity=tolower(gsub(pattern = "HILIC", "", .$Fraction1))) %>%
  mutate(m.z=as.numeric(m.z)) %>%
  select(compound_type=Compound.Type, compound_name=Compound.Name,
         compound_name_old=Compound.Name_old,
         formula=Emperical.Formula, rt=RT..min., mz=m.z, ionization_form,
         charge=z, kegg_id=C0, polarity, date_added=Date.added, mix=HILICMix) %>%
  add_row(compound_type="Custom",
          compound_name="Sulfobetaine", compound_name_old="Sulfobetaine",
          formula="C4H8O2S", rt=6.9, mz=120.024501+1.007276,
          ionization_form="[M+H]", charge=1, kegg_id=NA) %>%
  mutate(date_added=strptime(x = date_added, format = "%y%m%d"))


# And here are the tables! The anno_table has all of my annotation information:
#"compound_name" is my level 1 confidence annotation, using our internal standards and m/z + rt matching; 
# "name_deduced" is basically my level 2 confidence annotation, using m/z and MSMS data; 
# "confidence" expresses how certain I am I've got the correct name for it, with more negative numbers being less likely; 
# "classes" denotes the ClassyFire classes the compound is suspected to belong to, which can be quite extensive; 
# "formula_deduced" is kinda like a confidence level 2.5 where I'm pretty darn sure that's the MF but I'm only using m/z data to get it; 
# and "feature" is just an index column that I use to refer back to the actual peak areas in the other table.


# Then there's the actual peak table with area information: 
# "feature" is the index column same as above; 
# peak_id is another index column to a third, even bigger data frame; 
# mz, mzmin, mzmax are all hopefully intuitive - they're specific to each peak; 
# rt, rtmin rtmax also hopefully intuitive, in seconds; into = integrated area; 
# intb is mostly useless but hypothetically is "integration to baseline"; 
# maxo is max peak height; sn is signal-to-noise, with "good" peaks having sn values > 20;
# ignore sample; and file_name is the file the data came from. Whew.


# Typically I'll only be working with a subset of the data, 
# so I can select only the columns I'll need before left_joining the two together and plotting areas. 
# I have a fourth table containing cruise metadata that uses file_name as a key if you're interested in that one too.
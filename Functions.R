MakeScantable <- function(scan) {
  # Create a filtered MS2 scantable from a concatenated scanlist of MS2s.
  #
  # Args
  #   scan: Single observation from a dataframe containing MS2 mz and intensity 
  #         spectra, separated by semicolons (;).
  #
  # Returns
  #   scantable: Tiny dataframe, containing columns of mz and intensity. 
  #              Intensity is scaled to 100 and filtered to drop all intensity
  #              values below 0.5.
  #
  scantable <- read.table(text = as.character(scan),
                          col.names = c("mz", "intensity"), fill = TRUE) %>%
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", "")),
           intensity = round(intensity/max(intensity)*100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity))
  
  return(scantable)
}


MS2CosineSimilarity <- function(scan1, scan2) {
  # Finds the weighted cosine similarity between two sets of MS2 spectra.
  #
  # Args
  #   scan1 & scan2: Tiny dataframes of MS2. First column is mz, second column is intensity.
  #                  These are outputs of the MakeScantable() function.
  #
  # Returns
  #   cosine.similarity: A weighted similarity score between 0 and 1, indicating the cosine
  #                      relationship of the two vectors.
  #   
  mz.tolerance <- 0.02
  
  weight1 <- (scan1[, 1] ^ 2) * sqrt(scan1[, 2])
  weight2 <- (scan2[, 1] ^ 2) * sqrt(scan2[, 2])
  
  diff.matrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x) 
  same.index <- which(abs(diff.matrix) < mz.tolerance, arr.ind = TRUE)
  cosine.similarity <- sum(weight1[same.index[, 2]] * weight2[same.index[, 1]]) / 
    (sqrt(sum(weight2 ^ 2)) * sqrt(sum(weight1 ^ 2)))
  
  return(cosine.similarity)
}


# I'd look at MSMScosine_1_df to see how you can feed in a df of mass1, mass2, scan1 (in my format), scan2 and get an output of df$cosine.
MSMScosine1_df <- function(df) {
  scan1 <- scantable(df["scan1"])
  scan2 <- scantable(df["scan2"])
  mass1 <- df["mass1"]
  mass2 <- df["mass2"]
  
  mztolerance<-0.02
  
  w1<-(scan1[,1]^2)*sqrt(scan1[,2])
  w2<-(scan2[,1]^2)*sqrt(scan2[,2])
  
  diffmatrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  sameindex <- which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity)
}

MakeCandidates <- function(MoNA.Mass) {
  Candidates <- MoNA.Spectra.MHMass %>% 
    filter(MH_mass > MoNA.Mass - 0.020,  
           MH_mass < MoNA.Mass + 0.020) %>% 
    difference_inner_join(Experimental.Spectra.ForJoin, by = "MH_mass", max_dist = 0.02) %>% 
    mutate(scan1 = spectrum_KRHform_filtered, # scan1 is MS2 from MoNA
           scan2 = MS2, # scan2 is MS2 from the experimental data
           mass1 = MH_mass.x, # mass1 is the primary mass from MoNA
           mass2 = MH_mass.y) # mass2 is the primary mass from experimental data
  
  if (length(Candidates$ID) == 0) {
    print("There are no potential candidates.")
    No.Match.Return <- MF.Fraction %>%
      mutate(MassBankMatch = NA,
             MassBankppm = NA,
             MassBankCosine1 = NA)
    
    return(No.Match.Return)
  }
  
  # Add cosine similarity scores
  print("Making potential candidates")
  
  Candidates$Cosine1 <- apply(Candidates, 1, FUN=function(x) MSMScosine1_df(x)) 
  
  Candidates.Filtered.Cosine <- Candidates %>%
    filter(Cosine1 > Cosine.Score.Cutoff) %>%
    arrange(desc(Cosine1))
  
  Final.Candidates <- Candidates.Filtered.Cosine %>%
    mutate(MassBankMatch = paste(Names, ID, sep = " ID:"),
           MassBankppm = abs(mass1 - mass2) / mass2 * 10^6,
           MassBankCosine1 = Cosine1) %>%
    unique() %>%
    filter(MassBankppm < MassBank.ppm.Cutoff) #%>% 
  #select(MF.Fraction, MassBankMatch:MassBankCosine1)
  
  return(Final.Candidates)
}

MSMScosine1_df <- function(df) {
  scan1 <- scantable(df["scan1"])
  scan2 <- scantable(df["scan2"])
  mass1 <- df["mass1"]
  mass2 <- df["mass2"]
  print("scan1")
  print(scan1)
  print("scan2")
  print(scan2)
  
  mztolerance<-0.02
  
  w1<-(scan1[,1]^2)*sqrt(scan1[,2])
  w2<-(scan2[,1]^2)*sqrt(scan2[,2])
  
  diffmatrix <- sapply(scan1[, 1], function(x) scan2[, 1] - x)
  print("diffmatrix")
  print(diffmatrix)
  sameindex <- which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity)
}

scantable <- function(Scan) {
  Try <- read.table(text=as.character(Scan),
                    col.names = c("mz","intensity"), fill = TRUE) %>% 
    mutate(mz = as.numeric(mz %>% str_replace(",", "")),
           intensity = as.numeric(intensity %>% str_replace(";", ""))) %>%
    mutate(intensity = round(intensity/max(intensity)*100, digits = 1)) %>%
    filter(intensity > 0.5) %>%
    arrange(desc(intensity))
  
  return(Try)
}  


# oldfunctions ------------------------------------------------------------

## OLD Function definitions: TBD if they will all be useful ##
ChangeClasses <- function(df, start.column, end.column) {
  # Change specified columns from factors to numeric.
  #
  # Args
  #   df: MSDial dataframe containing sample columns.
  #
  # Returns
  #   df: MSDial dataframe, with specified columns changed to a numeric class. 
  for (i in c(start.column:end.column)) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

ChangeXClasses <- function(df) {
  # Identifies columns starting with X and changes their class to numeric.
  #
  # Args
  #   df: MSDial dataframe reorganized to drop all empty rows at the top.
  #
  # Returns
  #   df: MSDial dataframed with modified sample column classes.
  #
  col.test <- grepl("^X", names(df))
  for (i in which(col.test == TRUE)) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

CheckStandards <- function (df) {
  # Mutates a new column identifying standard run types, then prints number of unique run types.
  #
  # Args
  #   df: Dataset of containing a Replicate.Name column, pre-filtered to include only standard runs.
  #
  # Returns
  #   df.checked: Dataset with a new column describing run types, and a printed message stating how many 
  #               unique types there are.
  #
  df.checked <- df %>%
    mutate(Type = paste(Env = ifelse(str_detect(Replicate.Name, "StdsMix|InH2O"), "Standards", "Water"),
                        Matrix = ifelse(str_detect(Replicate.Name, "InMatrix"), "Matrix", "Water"), sep = "_"))
  
  print(paste("Number of standard run types:", length(unique(df.checked$Type))))
  print(unique(df.checked$Type))
  
  return(df.checked)
}

CheckStandards2 <- function (df) {
  # Mutates a new column identifying standard run types, then prints number of unique run types.
  #
  # Args
  #   df: Dataset of containing a Replicate.Name column, pre-filtered to include only standard runs.
  #
  # Returns
  #   df.checked: Dataset with a new column describing run types, and a printed message stating how many 
  #               unique types there are.
  #
  df.checked <- df %>%
    mutate(Type = paste(Env = ifelse(str_detect(Replicate.Name, "Stds"), "Standards", "Water"),
                        Matrix = ifelse(str_detect(Replicate.Name, "Matrix"), "Matrix", "Water"), sep = "_"))
  
  print(paste("Number of standard run types:", length(unique(df.checked$Type))))
  print(unique(df.checked$Type))
  
  return(df.checked)
}

IdentifyDuplicates <- function(df) {
  # If data is HILIC, Determine which compounds are detected in both positive and negative HILIC runs.
  # Otherwise, the function will return a printed message.
  # 
  # Args
  #   df: MSDial dataframe, containing all required parameters (MZ, SN, Area, etc),
  #       and modified to long form instead of wide.
  # 
  # Returns
  #   duplicates: Simple dataframe of listed compounds that have been identified as duplicates.
  #
  if ("Column" %in% colnames(df)) {
    duplicates <- df %>%
      group_by(Metabolite.Name, Replicate.Name) %>%
      mutate(number = 1) %>%
      mutate(ticker = cumsum(number)) %>%
      filter(ticker == 2) %>%
      ungroup() %>%
      select(Metabolite.Name) %>%
      unique()
    print("HILIC duplicates table created.")
    
    return(duplicates)
  } else {
    print("No instrument column data found.")
  }
}

IdentifyRunTypes <- function(df) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   df: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(df$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  if (length(skyline.runtypes != 4)) {
    stop("This run does not contain the four standard run types!")
  } else {
    print(paste("Your runtypes are:", toString(unique(run.type))))
  }
  return(run.type)
}

RearrangeDatasets <- function(df, parameter) {
  # Shortcut for altering multiple datasets using the tidyr::gather() function.
  #
  # Args
  #   df: MSDial dataframe with first n empty rows removed.
  #   parameter: Table value. This parameter will become the column name when 
  #              changed to long format.
  #
  # Returns
  #   df: MSDial dataframe, changed to long format and with a custom-named value column.
  df <- df %>%
    tidyr::gather(
      key = "Replicate.Name",
      value = "parameter",
      starts_with("X")) %>%
    select(Replicate.Name, parameter, everything())
  
  names(df)[2] <- parameter
  
  return(df)
}

RemoveCsv <- function(full.filepaths) {
  # Gathers all files in given directory and drops the csv extension.
  #
  # Args
  #   full.filepaths: list of files in a directory matching given patterns.
  #
  # Returns
  #   no.path: list of files, with filepath and csv extension removed.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

SetHeader <- function(df) {
  # Test for blank rows in MSDial output, and filter them out. Replace column names with syntactically correct headers.
  #
  # Args
  #   df: MSDial dataframe. 
  #
  # Returns
  #   df: MSDial dataframe, first n blank rows removed and headers set.
  #
  df <- df[!(is.na(df[1]) | df[1] == ""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}

StandardizeMetabolites <- function(df) {
  # Remove any "Ingalls_" prefixes that may be present in a dataframe.
  # Remove "X" prefixes in syntactically correct Replicate Names.
  #
  # Args
  #   df: MSDial dataframe.
  #
  # Returns
  #   df.standardized: Dataframe with above modifications.
  #
  df.standardized <- df %>%
    mutate(Metabolite.Name = ifelse(str_detect(Metabolite.Name, "Ingalls_"), 
                                    sapply(strsplit(Metabolite.Name, "_"), `[`, 2), Metabolite.Name)) 
  
  df.standardized$Replicate.Name <- gsub("^.{0,1}", "", df.standardized$Replicate.Name)
  
  return(df.standardized)
}

TrimWhitespace <- function (x) gsub("^\\s+|\\s+$", "", x)

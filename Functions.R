## Function definitions ##
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



CheckFragments <- function(skyline.file, runtype) { 
  # Modifies transformed skyline output to prepare for standard ion ratio detection by running several tests.
  # 1. Isolate standards.
  # 2. Isolate unique Product.Mz fragments per compound.
  # 3. Check if each compound has two unique fragments.
  # 4. If it does, identify which fragment is quantitative and which is secondary by comparing them to the master compound list.
  # 5. Find 5% of the quantitative fragment.
  # 6. Determine if the secondary trace is > the 5% value from step 5.
  #
  # Args:
  #   skyline.file: Output from skyline that has had its variables modified to numeric values.
  #
  # Returns:
  #   fragments.checked: Modified data frame with added columns reflecting the above tests.
  #
  fragment.check <- skyline.file %>%
    filter(str_detect(Replicate.Name, runtype)) %>%
    select(Replicate.Name, Precursor.Ion.Name, Area, Precursor.Mz, Product.Mz)
  
  fragment.unique <-unique(fragment.check %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))
  
  fragment.multi.unique <- fragment.unique %>%
    dplyr::count(Precursor.Ion.Name) %>%
    mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
    select(-n)
  
  fragments.checked <- fragment.check %>%
    left_join(fragment.multi.unique, by = "Precursor.Ion.Name") %>%
    merge(y = master.file,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    select(Replicate.Name, Precursor.Ion.Name, Area, Precursor.Mz, Product.Mz, Two.Fragments, Quan.Trace, Second.Trace) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    mutate(QT.Five.Percent = ifelse((Two.Fragments == TRUE & Quan.Trace == TRUE), 0.05 * Product.Mz, NA)) %>%
    mutate(Significant.Size = QT.Five.Percent < Product.Mz) %>%
    group_by(Replicate.Name, Precursor.Ion.Name, Product.Mz) %>%
    summarise_all(first) %>%
    arrange(Precursor.Ion.Name) %>%
    select(Replicate.Name, Precursor.Ion.Name, Precursor.Mz, Product.Mz, Area, Two.Fragments:Significant.Size)
  
  return(fragments.checked)
}

IdentifyRunTypes <- function(skyline.file) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   skyline.file: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(skyline.file$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
}



CheckSmpFragments <- function(skyline.file) {
  # Modifies transformed skyline output to prepare for standard ion ratio detection by running several tests.
  # 1. Isolate samples and pooled runs.
  # 2. Isolate unique Product.Mz fragments per compound.
  # 3. Check if each compound has two unique fragments.
  # 4. If it does, identify which fragment is quantitative and which is secondary by comparing them to the master compound list.
  # 5. Find 5% of the quantitative fragment.
  # 6. Determine if the secondary trace is > the 5% value from step 5.
  # 7. Find ion ratio by dividing the area of the quantitative trace by the area of the secondary trace.
  #
  # Args:
  #   areas.transformed: Output from skyline that has had its variables modified to numeric values.
  #
  # Returns:
  #   all.samples.IR: Modified data frame with added columns reflecting the above tests.
  unique.smp.frags <- unique(all.standards %>% select(Precursor.Ion.Name, Precursor.Mz, Product.Mz))
  
  unique.smp.frags2 <- unique.smp.frags %>%
    count(Precursor.Ion.Name) %>%
    mutate(Two.Fragments = ifelse((n==1), FALSE, TRUE)) %>%
    select(-n)
  
  all.samples.IR <- all.standards %>%
    left_join(unique.smp.frags2, by = "Precursor.Ion.Name" ) %>%
    merge(y = master.file,
          by.x = c("Precursor.Ion.Name", "Product.Mz"),
          by.y = c("Compound.Name", "Daughter"),
          all.x = TRUE) %>%
    select(Replicate.Name, Precursor.Ion.Name, Protein.Name, Precursor.Mz, Product.Mz, Area, Two.Fragments, Quan.Trace, Second.Trace) %>%
    mutate(Second.Trace = ifelse(Second.Trace == "", FALSE, TRUE)) %>%
    mutate(Quan.Trace = ifelse(Quan.Trace == "no", FALSE, TRUE)) %>%
    mutate(QT.Five.Percent = ifelse((Two.Fragments == TRUE & Quan.Trace == TRUE), 0.05 * Product.Mz, NA)) %>%
    mutate(Significant.Size = QT.Five.Percent < Product.Mz) %>%
    group_by(Precursor.Ion.Name) %>%
    mutate(IR.Ratio = ifelse(TRUE %in% Significant.Size, (Area[Quan.Trace == TRUE] / Area[Second.Trace == TRUE]), NA))
  
  return(all.samples.IR)
}

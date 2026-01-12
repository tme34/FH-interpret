library(shiny)
library(plotly)
library(data.table)
library(tidyverse)
library(bslib)


#Datasets
LDLR<-data.table::fread("data/LDLR_all_df_share.csv")
PCSK9<-data.table::fread("data/PCSK9_all_df_share.csv")
APOB<-data.table::fread("data/APOB_all_df_share.csv")

data_list = list(LDLR, PCSK9, APOB)


ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      /* 1. Unselected State */
      .btn-default {
        background-color: #f8f9fa !important;
        border-color: #ddd !important;
        color: #444 !important;
        
        /* INCREASED SPACING HERE */
        margin-right: 10px !important;  /* Horizontal gap */
        margin-bottom: 10px !important; /* Vertical gap if they wrap */
        
        border-radius: 8px !important;
        padding: 12px 24px !important;
        transition: all 0.2s ease-in-out !important;
        border-width: 2px !important;
        
        font-style: italic !important;
      }

      /* 2. Selected State */
      .btn-default.active, .btn-default.active:hover {
        background-color: #4A90E2 !important;
        border-color: #357ABD !important;
        color: white !important;
      }
      
      /* 3. Hover effect */
      .btn-default:hover {
        background-color: #e2e6ea !important;
      }

      /* Fix for spacing inside the button container */
      .btn-group-container-sw {
        display: flex;
        flex-wrap: wrap; /* Allows buttons to flow to next line if needed */
      }
    "))
  ),
  
  # 1. Add a title
  titlePanel("Familial hypercholesterolemia web-app"),
  
  # 2. Define the main layout
  sidebarLayout(
    # A. Sidebar panel for inputs
    sidebarPanel(
      width = 3,
      h3("FH gene browser"),
      shinyWidgets::radioGroupButtons
      (
        inputId = "dataset_choice",
        label = "Choose gene",
        choices = c("LDLR" = 1,"PCSK9" = 2, "APOB"= 3),
        selected = 1,
        individual = T,
        status = "default"
      ),
      
      h3("\nGraph appearance"),
      
      numericInput(inputId = "red_line",
                  label = "Elevation of red line above wild-type mean",
                  value = 0.50,
                  min = 0,
                  max = 3),
      
      checkboxInput(inputId = "show_red_line",
                    label = "Remove red line",
                    value = F), # Default to F
      
      checkboxInput(inputId = "show_percentile",
                    label = "Add statin-adjusted LDL percentiles for the population",
                    value = F), # Default to F
      
    ),
    
    # B. Main panel for outputs
    mainPanel(
      
      navset_tab(
      
      nav_panel("app", 
      h3("Automatic interpretation of the chosen variant"),
      shiny::htmlOutput("interpretOut"),
      shiny::plotOutput("graphOut"),
      DT::DTOutput("tableOut"),
      h6("The app is developed and maintained by Helene Gellert-Kristensen, Tim Møller Eyrich og Stefan Stender from the Department of Clinical Biochemistry at Rigshospitalet, Denmark. Contact: stefan.stender@regionh.dk")
      ),
      
      nav_panel("documentation", 
                h5("The app contains data for use in classifying variants in the LDLR, APOB, and PCSK9 genes. For the selected variant, the app generates a draft interpretation text. The app describes the available evidence but does not make suggestions for the final classification. The following data is included:"),
                h4(HTML("The Genome Aggregation Database, <a href='https://gnomad.broadinstitute.org/' target='_blank'>gnomAD</a>:")),
                h5("A database of publicly available aggregate data from large exome and genome sequencing projects from around the world. The interpretation text presents the total allele frequency for the selected variant in the approximately 1.6 million sequenced alleles currently available in gnomAD (v4.1.0). There is also a hyperlink to allele frequencies stratified by ethnicity for the selected variant (click on 'gnomAD')."),
                h4("UK Biobank:"),
                h5("UK Biobank is a large cohort of individuals from the British population. The app contains UK Biobank data for variants in the coding regions of the genes LDLR, APOB, and PCSK9 and their association with statin-adjusted plasma LDL-C in approximately 470,000 UK Biobank participants. The number of alleles and heterozygotes (but not homozygotes) as well as the mean LDL-C concentration for each genotype are described in the interpretation text. The effects on LDL-C are also presented graphically (mean +/-95% confidence interval) in the app itself.  Participants on statins have had their measured LDL-C multiplied by 1.42 to correct for the statin-induced reduction in LDL-C of an average of 30%. Statistical testing for differences in LDL-C between genotypes is via linear regression adjusted for sex, age, and ethnicity (principal genetic components 1-10). These analyses were conducted under UK Biobank Application Number: 104807."),
                h4(HTML("Global Lipids Genetics Consortium (GLGC) exome chip GWAS meta-analyse from 2017 (PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/29083408/' target='_blank'>29083408</a>):")),
                h5("Per-allele effects on LDL-C in approximately 300,000 individuals from 73 different cohorts. The UK Biobank was not included in the meta-analysis, so these data can be considered independent. Only a small group of variants in LDLR , APOB and PCSK9 were measured in this GWAS. For the variants that were not measured, the statement regarding these data is omitted in the app's interpretation text."),
                h4("AlphaMissense:"),
                h5(HTML("AlphaMissense is an AI tool developed by Google Deepmind that predicts the harmfulness of missense variants (PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/37733863/' target='_blank'>37733863</a>). AlphaMissense is based on AlphaFold's ability to predict protein structures, as well as evidence for how rarely variation in a specific amino acid occurs across humans and primates; changes in amino acids that show little or no variation are more likely to be harmful. AlphaMissense harmfulness is expressed as a value from 0 (benign) to 1 (pathogenic), with recommended cutoff values: <0.34: benign, 0.34-0.56: uncertain; >0.56: likely pathogenic. The closer to 0 or 1, respectively, the more certain the prediction.")),
                h4(HTML("Human Gene Mutation Database (<a href='https://www.hgmd.cf.ac.uk/ac/index.php' target='_blank'>HGMD</a>):")),
                h5("Is a curated database of genetic variants in disease-associated genes. For each variant, there is information about which articles the variant has appeared in, which phenotype it associates with, and the overall classification of the variant. The app includes data on all missense, nonsense and splice variants in LDLR, APOB, and PCSK9 that have been strongly associated with familial hypercholesterolemia (Variant Class 'DM', red box in HGMD). Latest extract: December 2025. For the selected variant, the PMIDs of articles referencing the variant are included as evidence in the interpretation text. For variants with more than 4 articles, 4 are mentioned, followed by 'and x other references'."),
                h4(HTML("<a href='https://www.ncbi.nlm.nih.gov/clinvar/' target='_blank' >ClinVar</a>:")),
                h5("A publicly available database of genetic variants and their association with diseases. All missense, nonsense, and splice variants in LDLR and PCSK9 that have been associated with familial hypercholesterolemia phenotypes have been extracted (latest extract: December 2025) and included as evidence in the app. Variant classification is divided into benign, likely benign, VUS, likely pathogenic, pathogenic, or 'unclear'. The app links to the ClinVar page for the selected variant in the interpretation response (click on ClinVar). ClinVar data for APOB is not included due to inconsistent naming of phenotype associations for variants in this gene (it is not possible to create a rule that reliably captures FH)."),
                h4("Other missense variants at the same amino acid position:"),
                h5("If another missense variant at the same amino acid position as the patient's variant is known to cause FH, this counts as moderate evidence of pathogenicity (ACMG criterion PM5). The app reports if other missense variants at the same position have been strongly associated with FH in HGMD (variant class 'DM') or in ClinVar (likely pathogenic or pathogenic)."),
                h4("Effect on LDL uptake and LDLR abundance on the cell surface:"),
                h5(HTML("In a study from 2025, the functional consequence of 17,000 different missense variants (nearly all theoretically possible) in LDLR were investigated (PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/41166440/' target='_blank'>41166440</a>). The variants were generated one at a time using CRISPR/CAS9 on modified HeLa cells. Each cell line with a specific mutation was then examined for cellular uptake of LDL-C (two different assays, one without and one with the presence of VLDL) and for abundance of the LDLR receptor on the cell surface. The function was quantified for each assay on a scale from 0 (no function) to 1 (normal, wild-type function). For each of the three assays, a score <0.5 is interpreted as impaired function. The app's interpretation text describes the overall result of these assays for the selected variant. Please note that the interpretation text only relates to data from this new study. Thus, older functional data that is not described in the app may exist for a given variant (such data may be identified through manual lookup in HGMD or ClinVar, which can be done using the links in the interpretation text).")),
                h5(HTML("The app is developed and maintained by MD, PhD <a href='https://scholar.google.com/citations?user=0zcd41YAAAAJ&hl=en&oi=ao/' target='_blank'>Helene Gellert-Kristensen</a>, cand. scient. <a href='https://scholar.google.com/citations?user=4H5xhzgAAAAJ&hl=en&oi=ao/' target='_blank'>Tim Møller Eyrich</a>, and MD, associate professor, PhD <a href='https://scholar.google.com/citations?user=mtgbiKoAAAAJ&hl=en&oi=ao/' target='_blank'>Stefan Stender</a>, all from the Department of Clinical Biochemistry at Rigshospitalet, Denmark. Contact: <a href='mailto:stefan.stender@regionh.dk'>stefan.stender@regionh.dk</a>"))
                )
      
      )
      
    )
  )
)


server <- function(input, output, session) {
  
  # --- 1. Reactive dataset selected by the user ---
  selected_dataset <- reactive({
    
    selected_data<-data_list[[as.numeric(input$dataset_choice)]]
    
    selected_data
  })
  
  selected_row_data <- reactive({
  
    selected_index <- input$tableOut_rows_selected
  
    # Check if a row is selected; otherwise, return NULL or the first row as a default
    if (length(selected_index) == 0) {
      # Return the first row as a default selection if none is selected
      return(selected_dataset()[1, ])
    } else {
      selected_dataset()[as.numeric(selected_index), ]
    }
  })
  
  
  output$tableOut <- DT::renderDT({
    
    data_df <- selected_dataset()
    req(data_df) 
    
    # Columns to EXCLUDE from rounding
    cols_to_exclude <- c("chr", "pos", "n0_carrier", "n1_carrier", "n2_carrier")
    
    data_df <- data_df %>%
      mutate(across(
        # 1. Select all numeric columns
        where(is.numeric) & 
          # 2. Exclude the specific columns
          !c(all_of(cols_to_exclude)), 
        ~ if_else(
          abs(.) < 0.01,         # Condition: Is the absolute value less than 0.01?
          signif(., digits = 3), # If TRUE: Apply signif with 3 digits
          round(., digits = 2) # If FALSE: Return the original value rounded(.)
        )
      ))
    
    # --- Column Definitions ---
    # Define the order of original column names to keep
    cols_to_keep_final_order <- c(
      # 1. Variant Identification & Location
      "chr",                    # Chromosome
      "pos",                    # Position Base Pair
      "ref",                    # Reference Allele
      "alt",                    # Alternative Allele
      "HGVS",                   # HGVS Notation (Start)
      "transcript",             # nucleotid position
      
      # 2. Sample Sizes
      "n0_carrier",                     # N Genotype 0
      "n1_carrier",                     # N Genotype 1
      "n2_carrier",                     # N Genotype 2
      
      # 3. LDL Means (Unadjusted & Adjusted)
      "mean_ldl_0",             # Mean LDL Genotype 0
      "se_ldl_0",               # Standard Error LDL Genotype 0
      "mean_ldl_1",             # Mean LDL Genotype 1
      "se_ldl_1",               # Standard Error LDL Genotype 1
      "mean_ldl_2",             # Mean LDL Genotype 2
      "se_ldl_2",               # Standard Error LDL Genotype 2
      "diff_1",                 # Difference Genotype 1 vs 0
      "diff_2",                 # Difference Genotype 2 vs 0
      "mean_ldl_0_adj",         # Adjusted Mean LDL Genotype 0
      "se_ldl_0_adj",           # Adjusted SE LDL Genotype 0
      "mean_ldl_1_adj",         # Adjusted Mean LDL Genotype 1
      "se_ldl_1_adj",           # Adjusted SE LDL Genotype 1
      "mean_ldl_2_adj",         # Adjusted Mean LDL Genotype 2
      "se_ldl_2_adj",           # Adjusted SE LDL Genotype 2
      "diff_1_adj",             # Adjusted Difference Genotype 1 vs 0
      "diff_2_adj",             # Adjusted Difference Genotype 2 vs 0
      
      # 4. Core Association Results (Beta and P-Value)
      "beta_genotype",          # Beta Effect Genotype
      "se_beta",                # Standard Error Beta Effect
      "p_value",                # P Value Association
      "beta_genotype_adj",      # Adjusted Beta Effect Genotype
      "se_beta_adj",            # Adjusted Standard Error Beta Effect
      "p_value_adj",            # Adjusted P Value Association
      
      # 5. External GWAS Results
      "rsid",                   # RS ID (rs number)
      "beta",                   # External Beta Effect
      "se",                     # External Standard Error
      "p",                      # External P Value
      
      # 6. Supplementary Scores
      "score_s1",               # Score Tool 1 Value
      "se_s1",                  # Score Tool 1 Standard Error
      "score_s2",               # Score Tool 2 Value
      "se_s2",                  # Score Tool 2 Standard Error
      "score_s3",               # Score Tool 3 Value
      "se_s3",                  # Score Tool 3 Standard Error
      
      # 7. Alpha missense pathogenicity
      "am_pathogenicity"
      
    )
    
    # Define the corresponding New Headers
    new_headers_final_order <- c(
      # 1. Variant Identification & Location
      "Chromosome",
      "Base Pair Position",
      "Reference Allele",
      "Alternative Allele",
      "HGVS Notation",
      "Nucleotide position",
      
      # 2. Sample Sizes
      "N Wild-type",
      "N Heterozygotes",
      "N Homozygotes",
      
      # 3. LDL Measurements (Unadjusted & Adjusted)
      "Mean LDL Wild-type",
      "SE LDL Wild-type",
      "Mean LDL Heterozygotes",
      "SE LDL Heterozygotes",
      "Mean LDL Homozygotes",
      "SE LDL Homozygotes",
      "Difference Heterozygotes vs Wild-type",
      "Difference Homozygotes vs Wild-type",
      "Adjusted Mean LDL Wild-type",
      "Adjusted SE LDL Wild-type",
      "Adjusted Mean LDL Heterozygotes",
      "Adjusted SE LDL Heterozygotes",
      "Adjusted Mean LDL Homozygotes",
      "Adjusted SE LDL Homozygotes",
      "Adjusted Difference Heterozygotes vs Wild-type",
      "Adjusted Difference Homozygotes vs Wild-type",
      
      # 4. Core Association Results (Beta and P-Value)
      "Beta Effect Genotype",
      "SE Beta Effect",
      "P Value Association",
      "Adjusted Beta Effect Genotype",
      "Adjusted SE Beta Effect",
      "Adjusted P Value Association",
      
      # 5. External GWAS Results
      "RS ID",
      "External Beta Effect",
      "External SE",
      "External P Value",
      
      # 6. Supplementary Scores
      "LDL Uptake Score",
      "LDL Uptake Score SE",
      "LDLR Cell-Surface Abundance Score",
      "LDLR Cell-Surface Abundance Score SE",
      "LDL Uptake in the Presence of VLDL Score",
      "LDL Uptake in the Presence of VLDL Score SE",
      
      # 7. alpha missense pathogenecity score
      "Alphamissense Pathogenicity"
    )
    
    # --- 2. APPLY CLEANUP LOGIC ---
    
    # Filter the list of columns to only include those present in the data
    cols_to_keep_present <- cols_to_keep_final_order[cols_to_keep_final_order %in% names(data_df)]
    
    # Filter the headers to match the present columns
    # We find the index of present columns in the original order list, 
    # and use those indexes to select the headers.
    header_indices <- match(cols_to_keep_present, cols_to_keep_final_order)
    headers_present <- new_headers_final_order[header_indices]
    
    cleaned_data <- data_df %>%
      
      # Use select(all_of()) to filter for the present columns in the desired order
      # This step implicitly removes ALL other columns, including the URLs, AC, TA, etc.
      select(all_of(cols_to_keep_present)) %>%
      
      # Use rename_with() to apply the new headers to the selected columns
      rename_with(~headers_present, all_of(cols_to_keep_present)) %>%
      mutate(value_num_n1 = ifelse(`N Heterozygotes` == "<5", 1, as.numeric(`N Heterozygotes`))) %>%
      mutate(value_num_n2 = ifelse(`N Homozygotes` == "<5", 1, as.numeric(`N Homozygotes`)))
    
    
    DT::datatable(cleaned_data, 
                                                 selection = "single",
                                                 options = list(
                                                   pageLength = 10, 
                                                   scrollX = TRUE,
                                                   # Center the table title
                                                   dom = 'Bfrtip',
                                                   columnDefs = list(
                                                     list(
                                                       targets = which(names(cleaned_data) == "N Heterozygotes"),
                                                       orderData = which(names(cleaned_data) == "value_num_n1")
                                                     ),
                                                     list(
                                                       targets = which(names(cleaned_data) == "N Homozygotes"),
                                                       orderData = which(names(cleaned_data) == "value_num_n2")
                                                     ),
                                                     list(
                                                       targets = 0,
                                                       visible = FALSE
                                                     ),
                                                     list(
                                                       targets = c("value_num_n1","value_num_n2"),
                                                       visible = FALSE
                                                     )
                                                    )
                                                   )
                                                
    )
  })
  
  
  output$graphOut <- renderPlot({
    
    # Get the filtered data
    plot_data <- selected_row_data()
    
    plot_data <- dplyr::select(
      plot_data, 
      mean_ldl_0_adj, se_ldl_0_adj, 
      mean_ldl_1_adj, se_ldl_1_adj, 
      mean_ldl_2_adj, se_ldl_2_adj
    ) %>%
      # Convert from wide to long format
      tidyr::pivot_longer(
        # Select all columns starting with 'mean_ldl' or 'se_ldl'
        cols = starts_with("mean_ldl") | starts_with("se_ldl"),
        # '.value' creates the new columns ('mean_ldl' and 'se_ldl'). 
        # 'genotype' holds the 0, 1, 2.
        names_to = c(".value", "genotype"),
        # Defines the column name structure: (content)_(genotype)
        names_pattern = "(.*)_(\\d)_adj" 
      ) %>%
      # Convert genotype to a factor for proper ordering
      dplyr::mutate(
        genotype = factor(genotype, levels = c("0", "1", "2"))
      )
    
    gg <- plot_data %>%
      ggplot() 
    
    gg<-
      if(input$show_red_line) {
        gg
      } else {
        gg+ geom_hline(yintercept = plot_data$mean_ldl[1]+input$red_line, color = "red", linetype = 2) 
      } 
     
    
    gg<-
    if(input$show_percentile) {
      
      gg + 
        
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 4.869, ymax = 5.264), fill = "red", alpha = 0.05, inherit.aes = FALSE) + 
        
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 4.438, ymax = 4.869), fill = "yellow", alpha = 0.05, inherit.aes = FALSE) + 
        
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 3.048, ymax = 4.438), fill = "green", alpha = 0.05, inherit.aes = FALSE) + 
        
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 2.729, ymax = 3.048), fill = "yellow", alpha = 0.05, inherit.aes = FALSE) + 
        
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 2.476, ymax = 2.729), fill = "red", alpha = 0.05, inherit.aes = FALSE) +
        
        
        
        # Adding Text Labels using annotate()
        
        # 1. 90 - 95 percentile (Top Red)
                annotate("text", x = 3.3, y = (4.869 + 5.264) / 2, 
                 
                 label = "90 - 95", color = "red", alpha = 0.8, fontface = "bold") +
                # 2. 80 - 90 percentile (Upper Yellow)
        
        annotate("text", x = 3.3, y = (4.438 + 4.869) / 2, 
                 label = "80 - 90", color = "darkgoldenrod3", alpha = 0.7, fontface = "bold") + # Use a darker yellow for visibility
        
        # 3. 20 - 80 percentile (Green)
        annotate("text", x = 3.3, y = (3.048 + 4.438) / 2, 
                  label = "20 - 80", color = "darkgreen", alpha = 0.7, fontface = "bold") + # Use a darker green for visibility
        
        # 4. 10 - 20 percentile (Lower Yellow)
        annotate("text", x = 3.3, y = (2.729 + 3.048) / 2, 
                 label = "10 - 20", color = "darkgoldenrod3", alpha = 0.7, fontface = "bold") +
        
        # 5. 5 - 10 percentile (Bottom Red)
        annotate("text", x = 3.3, y = (2.476 + 2.729) / 2, 
                 label = "5 - 10", color = "red", alpha = 0.8, fontface = "bold") +
        
        xlab("genotype")
      
      } else {
      gg
    } 
      
    #  geom_text(label = "0.5 mmol/L\nLDL-elevation", x = 3.3, y = plot_data$mean_ldl[1]+0.4, color = "red", size = 4) +
      gg +
      geom_errorbar(aes(x=genotype, ymin=mean_ldl-1.96*se_ldl, ymax=mean_ldl+1.96*se_ldl), width = 0.5) +
      geom_point(aes(x=genotype, y=mean_ldl), size = 8) +
      ylab(label = "Mean\nstatin-adjusted LDL (mmol/L)") +
      scale_x_discrete(labels = c("wild-type", "heterozygotes", "homozygotes")) +
      theme_minimal() +
      theme(text = element_text(size = 20))
    
  })


  
  interpretText<- reactive({
    

    variant  <- selected_row_data()$alph_variant
    
    #Definere chromosome ud fra variant navn
    prefix <- sub(":.*", "", variant)
    #Definere gen ud fra chromosome nummer
    gen <- ifelse(prefix == "19", "LDLR",
                  ifelse(prefix == "2",  "APOB",
                         ifelse(prefix == "1",  "PCSK9", NA)))
    #Assign df variabel ud fra gen navn
    if (gen == "LDLR") {
      df <- LDLR
    } else if (gen == "APOB") {
      df <- APOB
    } else if (gen == "PCSK9") {
      df <- PCSK9
      
    }
    #Assign score variabler (score 3 er ikke inkluderet endnu)
    s1<-df$score_s1[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    s2<-df$score_s2[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    s3<-df$score_s3[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    #Scores sættes til NA, hvis varianten ikke er i LDLR
    if (gen =="APOB" | gen == "PCSK9") {
      s1 <- NA
      s2 <- NA
      s3 <- NA
    }
    #Definere alpha missense score variabel
    am <- ifelse(is.na(df$am_pathogenicity[df$alph_variant==variant & !is.na(df$alph_variant==variant)]), "NA", df$am_pathogenicity[df$alph_variant==variant & !is.na(df$alph_variant==variant)])
    #Sætter am til empty string, hvis varianten er i APOB genet. (AM Scores findes ikke for APOB)
    if (gen == "APOB") {
      am <- ""
    }
    #am variablen skiftes til en string alt efter om scoren er over 0.564 (high pathogenicity).
    if (am>0.564 & am!="NA") {
      am <- paste("AlphaMissense predicts that the variant is deleterious (score: ",am,").", sep = "")
    } else {
      am <- ""
    }
    #hgvs variablen defineres og sættes som string hvis den eksisterer (ellers empty string).
    hgvs<-ifelse(is.na(df$HGVS[df$alph_variant==variant & !is.na(df$alph_variant==variant)]),"NA",df$HGVS[df$alph_variant==variant & !is.na(df$alph_variant==variant)])
    hgvs<- ifelse(hgvs=="NA","",paste(", ",hgvs, sep=""))
    #Definere antal heterozygote.
    het<-df$n1_carrier[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    het <- ifelse(is.na(het),"0",het)
    het_prox <- het
    #Definere total antal i UKbiobank med data på varianten.
    n <- df$ntotal[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    n <- ifelse(is.na(n),0,n)
    n_prox <- n
    #Definere Statin-justerede ldl værdier for vildtype og bærere for varianten.
    ldl_adj_wt<-round(df$mean_ldl_0_adj[df$alph_variant==variant & !is.na(df$alph_variant==variant)],2)
    ldl_adj_carrier<-round(df$mean_ldl_1_adj[df$alph_variant==variant & !is.na(df$alph_variant==variant)],2)
    #Definere forskellen mellem de 2 ldl variable.
    diff <- round(df$diff_1_adj[df$alph_variant==variant & !is.na(df$alph_variant==variant)],2)
    #Definere p-værdi for forskellen mellem bærere og vildtype.
    p <- df$p_value_adj[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    
    if (!is.na(p)){
      if (p > 0.01) {
        # Case 1: p is large (e.g., P = 0.02)
        # Returns: "P = 0.02"
        p<-round(p, digits = 2)
        
      } else {
        # Case 2: p is small, use scientific notation with HTML superscript
        
        val_s <- format(p, scientific = TRUE, digits = 1)
        
        superscript_map <- c(
          '0' = '\u2070', '1' = '\u00B9', '2' = '\u00B2', 
          '3' = '\u00B3', '4' = '\u2074', '5' = '\u2075', 
          '6' = '\u2076', '7' = '\u2077', '8' = '\u2078', 
          '9' = '\u2079'
        )
        
        if (grepl("e-", val_s)) {
          parts <- strsplit(val_s, "e-")[[1]]
          base <- parts[1]
          exponent<-sub("^0", "", parts[2])
          exponent <- strsplit(exponent, split = "")[[1]]
          superscript_exponenet <- superscript_map[exponent]
          superscript_exponenet <- paste(superscript_exponenet, collapse = "")
          
          p_text<-paste0(base, " x ", "10","\u207B",superscript_exponenet)}
        
        # Use HTML tags for superscript (<sup>) and the times symbol (&times;)
        # Returns: "P-trend = 4.5 &times; 10<sup>-6</sup>"
        p<-p_text
      }}
    else {
      # Fallback for unexpected format (e.g., positive exponent, which shouldn't happen with p-values)
      p <- p
    }
    
    #Definere transcript navn for varianten.
    t_name <- df$transcript[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    if (is.na(t_name)) {
      t_name <- df$variant[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    }
    #Definere variant allele count variabel.
    ac <- df$ac[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    ac <- ifelse(is.na(ac),0,ac)
    ac_prox<-ac
    ac <- format(ac, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
    #Definere variant total allele count variabel.
    ta <- df$ta[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    ta <- ifelse(is.na(ta),0,ta)
    ta_prox<-ta
    ta <- format(ta, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
    
    plus_minus <- ifelse(diff > 0 , "higher", "lower")
    diff <- abs(diff)
    ### Har jeg sat ind for at fikse APOB fejl.
    germ_class <- df$Germline.classification[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    if(gen=="APOB"){
      germ_class <- NA
    }
    ###
    if (het=="0" & n==0 | is.na(het) & is.na(n)) {
      ukb <- ""
    } else { 
      
      if (is.na(ldl_adj_carrier)) {
        het <- format(het, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
        n <- format(n, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
        
        ukb <- glue::glue ("In UK Biobank there are {het} heterozygotes out of ~470.000 participants.
      There are no LDL measurements for this variant in UK Biobank.")   
        
      } else {
        het <- format(het, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
        n <- format(n, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
        ukb <- glue::glue(
          "In UK Biobank, there are {het} heterozygotes out of ~470.000 participants. These {het} heterozygotes have on average {diff} mmol/L {plus_minus} statin-adjusted LDL-C than non-carriers ({ldl_adj_carrier} vs {ldl_adj_wt} mmol/L, p = {p})."
        )
      }
    }
    
    
    #If else conditional struktur til at definere korrekt string ud fra s1 og s2 scores.
    if (is.na(s1) & is.na(s3)) {
      if (is.na(s2)) {
        functional_statement <- ""
      } else {
        if (s2 > 0.5) {
          functional_statement <- "Functional studies have demonstrated an effect of the variant on LDLR cell-surface abundance. (PMID: 41166440)."
        } else {
          functional_statement <- "Functional studies have not demonstrated an effect of the variant on LDLR cell-surface abundance. (PMID: 41166440)."
        }
      }
    } else {
      if (s1 < 0.5 & is.na(s2)| isTRUE(s3 < 0.5) & is.na(s2)) {
        functional_statement <- "Functional studies have demonstrated an effect of the variant on LDL uptake (PMID: 41166440)."
      } else if (s1 > 0.5 & is.na(s2) | isTRUE(s3 > 0.5) & is.na(s2))  {
        functional_statement <- "Functional studies have not demonstrated an effect of the variant on LDL uptake (PMID: 41166440)."
      } else if (!is.na(s1) & !is.na(s2) | !is.na(s3) & !is.na(s2) ) {
        functional_statement <- if (s1 < 0.5 & s2 < 0.5 & !is.na(s1) & !is.na(s2) | s3 < 0.5 & s2 < 0.5 & !is.na(s3) & !is.na(s2)) { 
          "Functional studies have demonstrated an effect of the variant on LDL uptake and on LDLR cell-surface abundance (PMID: 41166440)."
        } else if (s1 > 0.5 & s2 > 0.5 & !is.na(s1) & !is.na(s2) | s3 > 0.5 & s2 > 0.5 & !is.na(s3) & !is.na(s2)) {
          "Functional studies have not demonstrated an effect of the variant on LDL uptake or LDLR cell-surface abundance (PMID: 41166440)."
        } else if (s1 < 0.5 & s2 > 0.5 & !is.na(s1) & !is.na(s2) | s3 < 0.5 & s2 > 0.5 & !is.na(s3) & !is.na(s2)) {
          "Functional studies have demonstrated an effect of the variant on LDL uptake (PMID: 41166440)."
        } else if (s1 > 0.5 & s2 < 0.5 & !is.na(s1) & !is.na(s2) | s3 > 0.5 & s2 < 0.5 & !is.na(s3) & !is.na(s2)) {
          "Functional studies have demonstrated an effect of the variant on LDLR cell-surface abundance (PMID: 41166440)."
        } else {
          ""
        }
      }
    }
    #PMID logic
    if (is.na(df$Ref_URL_5[df$alph_variant==variant & !is.na(df$alph_variant==variant)])) {
      if (is.na(df$Ref_URL_4[df$alph_variant==variant & !is.na(df$alph_variant==variant)])) {
        if (is.na(df$Ref_URL_3[df$alph_variant==variant & !is.na(df$alph_variant==variant)])) {
          if (is.na(df$Ref_URL_2[df$alph_variant==variant & !is.na(df$alph_variant==variant)])) {
            if (is.na(df$Ref_URL_1[df$alph_variant==variant & !is.na(df$alph_variant==variant)])) {
              pmid <- ""
            } else {
              pmid <- paste("The variant has been associated with FH in the literature (PMID: ",
                            df$Ref_URL_1[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                            ")",
                            ".",
                            sep="")
            }
          } else {
            pmid <- paste("The variant has been associated with FH in the literature (PMID: ",
                          df$Ref_URL_1[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                          ", ",
                          df$Ref_URL_2[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                          ")",
                          ".",
                          sep="")
          }
        } else {
          pmid <- paste("The variant has been associated with FH in the literature (PMID: ",
                        df$Ref_URL_1[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                        ", ",
                        df$Ref_URL_2[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                        ", ",
                        df$Ref_URL_3[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                        ")",
                        ".",
                        sep="")
        }
      } else {
        pmid <- paste("The variant has been associated with FH in the literature (PMID: ",
                      df$Ref_URL_1[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                      ", ",
                      df$Ref_URL_2[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                      ", ",
                      df$Ref_URL_3[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                      ", ",
                      df$Ref_URL_4[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                      ")",
                      ".",
                      sep="")
      }
    } else {
      pmid <- paste("The variant has been associated with FH in the literature (PMID: ",
                    df$Ref_URL_1[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                    ", ",
                    df$Ref_URL_2[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                    ", ",
                    df$Ref_URL_3[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                    ", ",
                    df$Ref_URL_4[df$alph_variant==variant & !is.na(df$alph_variant==variant)],
                    " and ",
                    gsub("andre reference\\(r\\)","more reference\\(s\\)",gsub("\\...$", "", df$Ref_URL_5[df$alph_variant==variant & !is.na(df$alph_variant==variant)])),
                    ").",
                    sep="")
    }
    
    
    if ((ac_prox==0 | ta_prox==0) & (het_prox=="0" | n_prox==0)) {
      gnomad <- glue::glue("{hgvs} is not found in 1.6 million alleles in gnomAD or among ~470.000 participants in UK Biobank.")
    } else if ((ac_prox!=0 | ta_prox!=0) & (het_prox=="0" | n_prox==0)){
      gnomad_var<- df$variant[df$alph_variant==variant & !is.na(df$alph_variant)]
      gnomad <- glue::glue(paste("{hgvs} is found in {ac} out of {ta} alleles in <a href='https://gnomad.broadinstitute.org/variant/",gnomad_var,"?dataset=gnomad_r4' target='_blank'>","gnomAD</a>, but not among ~470.000 participants in UK Biobank.", sep = ""))
      ukb <- ""
    } else if ((ac_prox==0 | ta_prox==0) & (het_prox!="0" & n_prox!=0)){
      gnomad <- glue::glue("{hgvs} is not found in 1.6 million alleles in gnomAD.")
    } else {
      gnomad_var<- df$variant[df$alph_variant==variant & !is.na(df$alph_variant)]
      gnomad <- glue::glue(paste("{hgvs} is found in {ac} out of {ta} alleles in <a href='https://gnomad.broadinstitute.org/variant/",gnomad_var,"?dataset=gnomad_r4' target='_blank'>","gnomAD</a>.", sep = ""))
    }
    
    clin_desc <- ""
    if (!is.na(germ_class)) {
      var_id <- df$VariationID[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
      if (germ_class=="Benign") {
        clin_desc <- paste("In <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>, the variant is classified as benign.", sep="")
      } else if (germ_class=="Benign/Likely benign") {
        clin_desc <- paste("In <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>, the variant is classified as benign / likely benign.", sep="")
      } else if (germ_class=="Conflicting classifications of pathogenicity") {
        clin_desc <- paste("The variant is reported with conflicting evidence of pathogenicity in <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>.", sep="")
      } else if (germ_class=="Likely benign") {
        clin_desc <- paste("In <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>, the variant is classified as likely benign.", sep="")
      } else if (germ_class=="Likely pathogenic") {
        clin_desc <- paste("In <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>, the variant is classified as likely pathogenic.", sep="")
      } else if (germ_class=="not provided") {
        clin_desc <- ""
      } else if (germ_class=="Pathogenic") {
        clin_desc <- paste("In <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>, the variant is classified as pathogenic.", sep="")
      } else if (germ_class=="Pathogenic/Likely pathogenic") {
        clin_desc <- paste("In <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>,the variant is classified as pathogenic/likely pathogenic.", sep="")
      } else if (germ_class=="Uncertain significance") {
        clin_desc <- paste("In <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",var_id,"/' target='_blank'>","ClinVar</a>, the variant is classified as uncertain significance.", sep="")
      } else {
        clin_desc <- ""
      }
    }
    
    others <- df$similar[df$alph_variant==variant & !is.na(df$alph_variant==variant)]
    
    if (others == 1) {
      others_desc <- "Another missense variant at the same amino acid position has been associated with FH."
    } else if (others > 1) {
      others_desc <- "Other missense variants at the same amino acid position have been associated with FH."
    } else {
      others_desc <- ""
    }
    
    glgc_rsid <- ifelse(!is.na(df$rsid[df$alph_variant==variant & !is.na(df$alph_variant==variant)]),df$rsid[df$alph_variant==variant & !is.na(df$alph_variant==variant)], NA)
    glgc_beta <- ifelse(!is.na(df$beta[df$alph_variant==variant & !is.na(df$alph_variant==variant)]),df$beta[df$alph_variant==variant & !is.na(df$alph_variant==variant)], NA)
    glgc_p <- ifelse(!is.na(df$p[df$alph_variant==variant & !is.na(df$alph_variant==variant)]),df$p[df$alph_variant==variant & !is.na(df$alph_variant==variant)], NA)
    if (!is.na(glgc_p)){
      if (glgc_p > 0.01) {
        # Case 1: p is large (e.g., P = 0.02)
        # Returns: "P = 0.02"
        glgc_p_text<-round(glgc_p, digits = 2)
        
      } else {
        # Case 2: p is small, use scientific notation with HTML superscript
        
        val_s <- format(glgc_p, scientific = TRUE, digits = 1)
        
        superscript_map <- c(
          '0' = '\u2070', '1' = '\u00B9', '2' = '\u00B2', 
          '3' = '\u00B3', '4' = '\u2074', '5' = '\u2075', 
          '6' = '\u2076', '7' = '\u2077', '8' = '\u2078', 
          '9' = '\u2079'
        )
        
        if (grepl("e-", val_s)) {
          parts <- strsplit(val_s, "e-")[[1]]
          base <- parts[1]
          exponent<-sub("^0", "", parts[2])
          exponent <- strsplit(exponent, split = "")[[1]]
          superscript_exponenet <- superscript_map[exponent]
          superscript_exponenet <- paste(superscript_exponenet, collapse = "")
          
          p_text<-paste0(base, " x ", "10","\u207B",superscript_exponenet)}
        
        # Use HTML tags for superscript (<sup>) and the times symbol (&times;)
        # Returns: "P-trend = 4.5 &times; 10<sup>-6</sup>"
        glgc_p_text<-p_text
      }}
    else {
      # Fallback for unexpected format (e.g., positive exponent, which shouldn't happen with p-values)
      glgc_p_text <- glgc_p
    }
    
    glgc<-""
    if (!is.na(glgc_p)) {
      if(glgc_p>=0.05) {
        glgc_beta <- round(glgc_beta, 2)
        glgc<-glue::glue(
          "The variant did not associate with LDL-C (beta: {glgc_beta}, p = {glgc_p_text}) in a meta-analysis of ~300,000 individuals that did not include UK Biobank (PMID: 29083408)."
        )}
      else if(glgc_p<0.05 & glgc_beta>0) {
        glgc_beta <- round(glgc_beta, 2)
        glgc<-glue::glue(
          "The variant was associated with increased LDL-C (beta: {glgc_beta}, p = {glgc_p_text}) in a meta-analysis of ~300,000 individuals, which did not include UK Biobank (PMID: 29083408)."
        )}
      else if(glgc_p<0.05 & glgc_beta<0) {
        glgc_beta <- round(glgc_beta, 2)
        glgc<-glue::glue(
          "The variant was associated with reduced LDL-C (beta: {glgc_beta}, p = {glgc_p_text}) in a meta-analysis of ~300,000 individuals, which did not include UK Biobank (PMID: 29083408)."
        )}
      else {glgc<-""}
    }
    
    
    
    #Definere endelig string ved brug af glue pakken, hvor de forskellige definerede variable indsættes.
    text <- glue::glue(
      "{em(gen)}: {t_name} {gnomad} {ukb} {glgc} {am} {others_desc} {functional_statement} {pmid} {clin_desc}")
    
    text
  })
  
  
output$interpretOut <- renderText({
  interpretText()
})



}

shinyApp(ui = ui, server = server)
# Load libraries
library(shiny)
library(UniprotR)
library(Biostrings)
library(forcats)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(shinydashboard)
library(thematic)
library(bio3d)
library(RColorBrewer)
library(stringdist)
library(shinyjs)
library(NGLVieweR)
library(pwalign)

# MUSCLE path
Sys.setenv(MUSCLE_BIN = "/usr/local/bin/muscle")

# URL existence checker
url.exists <- function(url) {
  tryCatch({ httr::HEAD(url)$status_code == 200 }, error = function(e) FALSE)
}

# Multi-line sequence formatter
format_alignment_multiline <- function(ref_seq, target_seq, width = 60) {
  lines <- ceiling(nchar(ref_seq) / width)
  output <- ""
  for (i in seq_len(lines)) {
    start <- (i - 1) * width + 1
    end <- min(i * width, nchar(ref_seq))
    ref_sub <- substr(ref_seq, start, end)
    tgt_sub <- substr(target_seq, start, end)
    color_line <- mapply(function(r, t) {
      if (r == t) paste0("<span style='color:green'><b>", r, "</b></span>")
      else if (r == "-" || t == "-") paste0("<span style='color:grey'><b>", t, "</b></span>")
      else paste0("<span style='color:red'><b>", t, "</b></span>")
    }, strsplit(ref_sub, "")[[1]], strsplit(tgt_sub, "")[[1]])
    output <- paste0(output, "<br><b>Ref:</b> ", ref_sub,
                     "<br><b>Tgt:</b> ", paste0(color_line, collapse = ""), "<br><br>")
  }
  return(output)
}

# UI
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "AbBridge – Bridging Antibodies Across Species", titleWidth = 450),
  dashboardSidebar(width = 450,
                   sidebarMenu(
                     menuItem("Input Reference sequence",
                              menuSubItem(tabName = "Ref",
                                          selectInput("ref_inputType", "Select Reference Input Type",
                                                      choices = c("", "Uniprot Accession ID", "File")
                                          )
                              ),
                              conditionalPanel(
                                condition = "input.ref_inputType == 'Uniprot Accession ID'",
                                textInput("uniprot", "Enter Reference Uniprot ID", placeholder = "e.g. P01308")
                              ),
                              conditionalPanel(
                                condition = "input.ref_inputType == 'File'",
                                fileInput("fileInput", "Upload Reference Fasta File:", accept = c("text/plain"))
                              ),
                              textInput("epitope_input", "Epitope sequence (optional)", "")
                     ),
                     menuItem("Input Target sequences",
                              menuSubItem(tabName = "Targets",
                                          selectInput("Target_inputType", "Select Target Input Type:",
                                                      choices = c("", "Uniprot Accession IDs", "Fasta files", "Both")
                                          )
                              ),
                              conditionalPanel(
                                condition = "input.Target_inputType == 'Fasta files' || input.Target_inputType == 'Both'",
                                fileInput("filesinput", "Upload Target Fasta File(s)", multiple = TRUE)
                              ),
                              conditionalPanel(
                                condition = "input.Target_inputType == 'Uniprot Accession IDs' || input.Target_inputType == 'Both'",
                                textAreaInput("uniprot_ids", "Enter Target Uniprot IDs", width = '100%',
                                              placeholder = "Separate by commas. Ex: P01308,P01317")
                              )
                     ),
                     actionButton(inputId = "align", label = "Compare Structures")
                   )
  ),
  dashboardBody(
    fluidPage(
      useShinyjs(),
      tabsetPanel(
        tabPanel("Home",
                 h3("Instructions"),
                 p("Use the sidebar to upload sequences and click 'Compare Structures'. Then navigate to view results in the tabs above."),
                 br(),
                 h4("🧪 Example UniProt IDs to try:"),
                 p("Reference ID: ", code("O00206")),
                 p("Target IDs: ", code("Q9GL65,Q9QUK6,G1SH24")),
                 actionButton("demo_button", "Load Demo Example")
        ),
        tabPanel("Alignment Results",
                 titlePanel("Alignment Table and Plot"),
                 fluidRow(column(12, tableOutput("Aln_table"))),
                 fluidRow(column(12, plotOutput("plot_out"))),
                 fluidRow(column(2, downloadButton("download", "Download Plot")))
        ),
        tabPanel("Detailed Results", 
                 h3("Per-target Alignment Details"), 
                 uiOutput("detailed_output")
        ),
        tabPanel("3D Structure",
                 h3("3D Structure Viewer"),
                 uiOutput("ngl_ui"),
                 p("Rendering 3D structure of proteins.")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  aln_data <- reactiveVal(NULL)
  ref_details <- reactiveVal(list())
  alignment_html <- reactiveVal(list())
  target_meta <- reactiveVal(NULL)
  ref_meta <- reactiveVal(NULL)
  
  ### Demo button observer ###
  observeEvent(input$demo_button, {
    
    updateSelectInput(session, "ref_inputType", selected = "Uniprot Accession ID")
    updateTextInput(session, "uniprot", value = "O00206")
    updateSelectInput(session, "Target_inputType", selected = "Uniprot Accession IDs")
    updateTextAreaInput(session, "uniprot_ids", value = "Q9GL65,Q9QUK6,G1SH24")
    updateTextInput(session, "epitope_input", value = "FLPDIFTEL")
    runjs("setTimeout(function() { $('#align').prop('disabled', false); }, 1000);")
  })
  
  observeEvent(input$align, {
    withProgress(message = "Processing...", value = 0.1, {
      
      # --- Validate and clean reference input ---
      ref_id <- NULL
      if (input$ref_inputType == "Uniprot Accession ID") {
        ref_id <- trimws(input$uniprot)
        if (is.null(ref_id) || ref_id == "" || grepl("[^A-Za-z0-9]", ref_id)) {
          showNotification("Invalid or missing reference UniProt ID.", type = "error")
          return()
        }
      }
      
      # --- Validate and clean target IDs ---
      raw_target_ids <- trimws(input$uniprot_ids %||% "")
      target_ids_split <- if (raw_target_ids != "") strsplit(raw_target_ids, ",")[[1]] %>% trimws() else character(0)
      if (any(grepl("[^A-Za-z0-9]", target_ids_split))) {
        showNotification("Invalid characters in target UniProt IDs.", type = "error")
        return()
      }
      
      Sys.sleep(0.2)  # give inputs time to settle
      incProgress(0.2, detail = "Getting targets...")
      # --- Fetch reference sequence and metadata ---
      if (input$ref_inputType == "File") {
        req(input$fileInput)
        ref_seq <- readAAStringSet(input$fileInput$datapath)
        merged_ref_df <- data.frame(Gene.Names = names(ref_seq), Organism = "Uploaded File")
        ref_name <- names(ref_seq)[1]
        ref_org <- "Uploaded File"
        ref_fullname <- "NA"
      } else {
        refseq_df <- tryCatch({
          GetSequences(ref_id)
        }, error = function(e) {
          showNotification("Failed to retrieve reference sequence.", type = "error"); NULL })
        if (is.null(refseq_df)) return()
        
        reftaxa_df <- tryCatch({
          GetNamesTaxa(ref_id)
        }, error = function(e) {
          showNotification("Failed to retrieve reference taxonomy.", type = "error"); NULL })
        if (is.null(reftaxa_df)) return()
        
        merged_ref_df <- merge(refseq_df, reftaxa_df, by = 'row.names', all = TRUE)
        ref_seq <- AAStringSet(merged_ref_df$Sequence)
        names(ref_seq) <- paste0(merged_ref_df$Entry, "_", merged_ref_df$Organism)
        ref_name <- merged_ref_df$Gene.Names[1]
        ref_org <- merged_ref_df$Organism[1]
        ref_fullname <- merged_ref_df$Protein.names[1]
      }
      
      ref_details(list(id = ref_id, name = ref_name,
                       fullname = ref_fullname, len = width(ref_seq), org = ref_org))
      ref_meta(merged_ref_df)
      
      # --- Fetch target sequences and metadata ---
      target_seqs <- AAStringSet()
      meta_df <- data.frame()
      
      if (input$Target_inputType %in% c("Fasta files", "Both")) {
        req(input$filesinput)
        fasta_seqs <- readAAStringSet(input$filesinput$datapath)
        target_seqs <- c(target_seqs, fasta_seqs)
        meta_df <- bind_rows(meta_df,
                             data.frame(Entry = names(fasta_seqs),
                                        Gene.Names = names(fasta_seqs),
                                        Organism = "Uploaded File",
                                        Protein.names = names(fasta_seqs)))
      }
      
      if (length(target_ids_split) > 0) {
        seq_df <- tryCatch({
          GetSequences(target_ids_split)
        }, error = function(e) {
          showNotification("Failed to retrieve one or more target sequences.", type = "error"); NULL })
        if (is.null(seq_df)) return()
        
        taxa_df <- tryCatch({
          GetNamesTaxa(target_ids_split)
        }, error = function(e) {
          showNotification("Failed to retrieve one or more target taxonomy records.", type = "error"); NULL })
        if (is.null(taxa_df)) return()
        
        merged_df <- merge(seq_df, taxa_df, by = 'row.names', all = TRUE)
        uni_seqs <- AAStringSet(merged_df$Sequence)
        names(uni_seqs) <- paste0(merged_df$Entry, "_", merged_df$Organism)
        target_seqs <- c(target_seqs, uni_seqs)
        target_meta_df <- merged_df %>%
          select(Entry, Gene.Names, Organism, Protein.names)
        meta_df <- bind_rows(meta_df, target_meta_df)
      }
      
      target_meta(meta_df)
      
      # Pairwise Sequence Alignment---------------------------------------------
      incProgress(0.2, detail = "Aligning sequences...")
      local_aln <- sapply(target_seqs, pairwiseAlignment, ref_seq, type = "local", substitutionMatrix = "BLOSUM50", scoreOnly = FALSE)
      global_aln <- sapply(target_seqs, pairwiseAlignment, ref_seq, type = "global", substitutionMatrix = "BLOSUM50", scoreOnly = FALSE)
      
      # Epitope Coverage Score -------------------------------------------------
      epitope_seq <- gsub("\\s+", "", input$epitope_input)
      epitope_scores <- rep(NA, length(target_seqs))
      if (epitope_seq != "") {
        for (i in seq_along(target_seqs)) {
          t_seq <- as.character(target_seqs[[i]])
          w_size <- nchar(epitope_seq)
          dist_vals <- sapply(1:(nchar(t_seq) - w_size + 1), function(pos) {
            stringdist(substr(t_seq, pos, pos + w_size - 1), epitope_seq, method = "lv")
          })
          best_match <- min(dist_vals)
          coverage <- 100 * (1 - (best_match / w_size))
          epitope_scores[i] <- ifelse(coverage < 50, 0, round(coverage, 2))
        }
      }
      
      # 3D protein structure RMSD evaluation -------------------------------------
      incProgress(0.2, detail = "Computing RMSD...")
      rmsd_results <- list()
      fixed_id <- input$uniprot
      ids <- unlist(strsplit(input$uniprot_ids, ","))
      for (uniprot_id2 in ids) {
        url1 <- paste0("https://alphafold.ebi.ac.uk/files/AF-", fixed_id, "-F1-model_v4.pdb")
        url2 <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot_id2, "-F1-model_v4.pdb")
        if (!url.exists(url1) || !url.exists(url2)) next
        download.file(url1, destfile = paste0(fixed_id, "_AF.pdb"), quiet = TRUE)
        download.file(url2, destfile = paste0(uniprot_id2, "_AF.pdb"), quiet = TRUE)
        pdb1 <- read.pdb(paste0(fixed_id, "_AF.pdb"))
        pdb2 <- read.pdb(paste0(uniprot_id2, "_AF.pdb"))
        pdb1.inds <- atom.select(pdb1, elety = c("CA", "CB", "CG", "CG1", "CG2"))
        pdb2.inds <- atom.select(pdb2, elety = c("CA", "CB", "CG", "CG1", "CG2"))
        aln <- struct.aln(pdb1, pdb2, pdb1.inds, pdb2.inds)
        rmsd_val <- aln$rmsd[length(aln$rmsd)]
        taxa <- GetNamesTaxa(uniprot_id2)
        pair_name <- paste0(uniprot_id2, "_", taxa$Organism)
        rmsd_results[[pair_name]] <- rmsd_val
      }
      
      rmsd_numeric <- unlist(rmsd_results)
      rmsd_scaled <- round(100 * (1 - ((rmsd_numeric - 0)/(5))), 2)
      
      # 3D protein structure visualization ---------------------------------------
      incProgress(0.2, detail = "Rendering 3D structure...")
      output$ngl_ui <-  renderUI({
        req(input$uniprot, input$uniprot_ids, rmsd_scaled, aln_data())
        
        # Reference protein visualization
        reference = input$uniprot
        
        # check for the pdb file in directory path
        file_path <- file.path( paste0(reference, "_AF.pdb"))
        
        if (!file.exists(file_path)) 
        { return(NULL)  }
        taxas <- GetNamesTaxa(reference)
        refrence_box <- box(
          title = paste( "Reference protein :", reference,"-", taxas$Organism,",",taxas$Gene.Names), width = 12, solidHeader = TRUE, status = "success",
          NGLVieweROutput(outputId = paste0("ngl_plot_", reference), height = "400px")
        )
        
        
        # target proteins ranked based on table
        tabledata <- aln_data()
        # use the Rank in table to sort target proteins
        ranked_ids <- tabledata %>% arrange(Rank) %>% pull(Accession_ID_Species_name)
        
        target_boxes <- lapply(ranked_ids,function(acession_id_name) {
          
          # split Accesion ID name to get Uniprot ID only
          id <- strsplit(acession_id_name,"_")[[1]][1]
          
          rank <- tabledata %>% 
            filter(Accession_ID_Species_name == acession_id_name) %>% 
            pull(Rank)
          
          
          # check for the pdb file in directory path
          file_path <- file.path( paste0(id, "_AF.pdb"))
          
          if (!file.exists(file_path)) 
          { return(NULL)  }
          taxas <- GetNamesTaxa(id)
          box(
            title = paste("Target Protein Rank ", rank,': ',id,"-", taxas$Organism,",",taxas$Gene.Names,".   Scaled RMSD: ",rmsd_scaled[paste0(id, "_", taxas$Organism)]),
            width = 12, solidHeader = TRUE, status = "success",
            NGLVieweROutput(outputId = paste0("ngl_plot_", id), height = "400px")
          )
        })
        # this is important to list all boxes in tab
        tagList(refrence_box,target_boxes)
      })
      
      observe({
        req(input$uniprot, input$uniprot_ids)
        all_ids <- unique(c(input$uniprot, unlist(strsplit(input$uniprot_ids, ","))))
        
        for (id in all_ids) {
          local_path <- file.path(paste0(id, "_AF.pdb"))
          
          local({
            my_id <- id
            my_path <- local_path
            
            output[[paste0("ngl_plot_", my_id)]] <- renderNGLVieweR({
              req(file.exists(my_path))
              NGLVieweR(my_path) %>%
                stageParameters(backgroundColor = "white") %>%
                setQuality("high") %>%
                setSpin(FALSE) %>%
                addRepresentation("cartoon",
                                  param = list(
                                    name = "cartoon",
                                    colorScheme = "residueindex"
                                  ))
            })
          })
        }
      })
      
      
      
      
      # Table of alignment scores and percent Identity ---------------------------
      incProgress(0.2, detail = "Generating results...")
      df <- data.frame(
        Global_Aln_Score = sapply(global_aln, score),
        Global_PID = round(sapply(global_aln, pid, type = "PID3")),
        Local_Aln_Score = sapply(local_aln, score),
        Local_PID = round(sapply(local_aln, pid, type = "PID3"), 0),
        Epitope_Coverage = epitope_scores
      ) %>% rownames_to_column("Accession_ID_Species_name")
      
      df$RMSD <- rmsd_numeric[df$Accession_ID_Species_name]
      df$RMSD_Scaled <- rmsd_scaled[df$Accession_ID_Species_name]
      df <- df %>% arrange(desc(Epitope_Coverage)) %>% mutate(Rank = row_number())
      aln_data(df)
      alignment_html(lapply(names(target_seqs), function(nm) {
        alignment <- pairwiseAlignment(ref_seq[[1]], target_seqs[[nm]], type = "global", substitutionMatrix = "BLOSUM50")
        format_alignment_multiline(as.character(pattern(alignment)), as.character(subject(alignment)))
      }))
    })
  })
  
  output$Aln_table <- renderTable({ req(aln_data()); aln_data() })
  output$plot_out <- renderPlot({
    req(aln_data())
    reference_taxa <- GetNamesTaxa(input$uniprot)
    aln_data() %>%
      pivot_longer(cols = c("Global_PID", "Local_PID", "RMSD_Scaled", "Epitope_Coverage"),
                   names_to = "Metrics", values_to = "Percentage") %>%
      mutate(Accession_ID_Species_name = fct_reorder(Accession_ID_Species_name, -Percentage)) %>%
      ggplot(aes(x = Accession_ID_Species_name, y = Percentage / 100, fill = Metrics)) +
      geom_col(position = position_dodge2()) +
      coord_flip() +
      geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(0.9), hjust = 0) +
      scale_fill_brewer(palette = "RdYlBu") +
      theme_classic() +
      labs(title = paste0(" Species Antibody evaluation scores for ", reference_taxa$Gene.Names," of ", reference_taxa$Organism),
           y = "Scores",
           x = "Target Proteins")
  })
  
  output$download <- downloadHandler(
    filename = function () { paste0("AbBridge_Alignment_Plot.png") },
    content = function(file) {
      ggsave(file, plot = output$plot_out(), device = "png", width = 18, height = 12)
    }
  )
  #Detailed Alignment Results
  output$detailed_output <- renderUI({
    req(aln_data())
    ref <- ref_details()
    df <- aln_data()
    meta <- target_meta()
    html_output <- list(
      h4(paste("Reference Protein ID:", ref$id)),
      h4(paste("Gene Name:", ref$name)),
      h4(paste("Full Name:", ref$fullname)),  
      h4(paste("Length:", ref$len)),
      h4(paste("Organism:", ref$org)),
      hr(),
      h4("Alignment Color Legend"),
      HTML("<ul>
           <li><span style='color:green'><b>Green</b></span>: Match</li>
           <li><span style='color:red'><b>Red</b></span>: Mismatch</li>
           <li><span style='color:grey'><b>Gray</b></span>: Gap</li>
         </ul>"),
      hr()
    )
    alignments <- alignment_html()
    names(alignments) <- df$Accession_ID_Species_name
    target_index <- 1
    for (nm in names(alignments)) {
      row <- df %>% filter(Accession_ID_Species_name == nm)
      entry <- strsplit(nm, "_")[[1]][1]
      organism <- strsplit(nm, "_")[[1]][2]
      gene_name <- meta$Gene.Names[meta$Entry == entry]
      full_name <- meta$Protein.names[meta$Entry == entry]  
      target_len <- nchar(gsub("-", "", gsub("<[^>]+>", "", gsub("<br>", "", alignments[[nm]]))))
      html_output <- append(html_output, list(
        wellPanel(
          h4(paste0("Target ", target_index, ": ", entry)),
          p(paste("Gene Name:", gene_name)),
          p(paste("Full Name:", full_name)),  
          p(paste("Organism:", organism)),
          p(paste("Target Length:", target_len)),
          p(paste("Local Alignment Score:", row$Local_Aln_Score)),
          p(paste("Global PID:", row$Global_PID, "%")),
          p(paste("Rank:", row$Rank)),
          HTML(alignments[[nm]])
        )
      ))
      target_index <- target_index + 1
    }
    do.call(tagList, html_output)
  })
}

shinyApp(ui = ui, server = server)

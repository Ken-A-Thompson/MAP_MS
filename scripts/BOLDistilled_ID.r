# ---- Setup ----
packages <- c(
  "shiny",
  "shinyFiles",
  "seqinr",
  "data.table"
)

installed <- rownames(installed.packages())

for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Allow access to entire filesystem
roots <- c(home = normalizePath("~"))

# ---- Helper Functions ----
clean_fasta <- function(fasta_file, cleaned_file) {
  seqs <- read.fasta(fasta_file, as.string = TRUE, set.attributes = FALSE)
  cleaned_seqs <- lapply(seqs, function(s) {
    s <- toupper(s)
    s <- gsub("-", "", s)
    s <- gsub("[IRYSWKMBDHV]", "N", s)
    return(s)
  })
  write.fasta(sequences = cleaned_seqs, names = names(seqs), file.out = cleaned_file)
}

run_vsearch <- function(fasta, db, method, threads = 2) {
  out_file <- tempfile(fileext = ".txt")
  if (method == "usearch_global") {
    cmd <- paste(
      "vsearch",
      "--usearch_global", shQuote(fasta),
      "--db", shQuote(db),
      "--blast6out", shQuote(out_file),
      "--id 0.75 --maxhits 5 --maxaccepts 5",
      "--threads", threads
    )
  } else if (method == "sintax") {
    cmd <- paste(
      "vsearch",
      "--sintax", shQuote(fasta),
      "--db", shQuote(db),
      "-tabbedout", shQuote(out_file),
      "-strand plus -sintax_cutoff 0.6",
      "--threads", threads
    )
  } else stop("Unknown method")
  system(cmd)
  return(out_file)
}

parse_usearch_global <- function(file, taxonomy_file) {
  df <- fread(file, sep="\t", header=FALSE, data.table = FALSE)
  if (ncol(df) < 4) return(data.frame(Error = "Unexpected VSEARCH output"))
  colnames(df)[1:4] <- c("Query","Hit","Pct_ID","Overlap_bp")
  df <- df[df$Overlap_bp >= 500, ]
  if (nrow(df) == 0) { return(data.frame()) }
  df <- df[order(df$Query, -df$Pct_ID, -df$Overlap_bp), ]
  df <- df[!duplicated(df$Query), ]
  df$BIN <- sapply(strsplit(df$Hit, "\\|"), function(x) x[2])
  tax <- fread(taxonomy_file, sep="\t", header=TRUE, data.table = FALSE)
  df <- merge(df, tax, by.x="BIN", by.y="bin", all.x=TRUE)
  df$BIN_Match <- ifelse(df$Pct_ID >= 97.7, "BIN_MATCH", "NO_MATCH")
  df <- df[,c(2,26,1,4,5,14:23)]
  return(df)
}

parse_sintax <- function(file) {
  df <- read.table(file, sep="\t", header=FALSE, stringsAsFactors = FALSE, quote = "")
  if (ncol(df) < 4) return(data.frame(Error = "Unexpected SINTAX output"))
  query <- df$V1
  raw_tax <- df$V4
  tax_df <- data.frame(
    Query = query,
    kingdom = "", phylum = "", class = "", order = "",
    family = "", genus = "", species = "", stringsAsFactors = FALSE
  )
  for (i in seq_along(raw_tax)) {
    ranks <- unlist(strsplit(raw_tax[i], ","))
    for (r in ranks) {
      r_clean <- sub("^[a-z]:", "", r)
      if (grepl("^k:", r)) tax_df$kingdom[i] <- r_clean
      if (grepl("^p:", r)) tax_df$phylum[i] <- r_clean
      if (grepl("^c:", r)) tax_df$class[i] <- r_clean
      if (grepl("^o:", r)) tax_df$order[i] <- r_clean
      if (grepl("^f:", r)) tax_df$family[i] <- r_clean
      if (grepl("^g:", r)) tax_df$genus[i] <- r_clean
      if (grepl("^s:", r)) tax_df$species[i] <- r_clean
    }
  }
  return(tax_df)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("BOLDistilled ID (local computer)"),
  
  # Use fluidRow to allow panels to take full width
  fluidRow(
    column(width = 12,
           wellPanel(
             tags$h3("1. FASTA Input"),
             radioButtons("fasta_input_type", strong("How will you provide the FASTA?"),
                          choices = c("Select FASTA file" = "file",
                                      "Paste FASTA text" = "paste"),
                          selected = "file"),
             conditionalPanel(
               condition = "input.fasta_input_type == 'file'",
               shinyFilesButton("fasta_file", "📁 Select FASTA file", "Select FASTA", multiple = FALSE),
               tags$small(em(textOutput("fasta_path")))
             ),
             conditionalPanel(
               condition = "input.fasta_input_type == 'paste'",
               textAreaInput("fasta_text", "Paste FASTA here:", rows = 10, width = "100%")
             )
           )
    )
  ),
  
  fluidRow(
    column(width = 12,
           wellPanel(
             tags$h3("2. Select Query Method"),
             radioButtons("method", strong("Choose method"),
                          choices = c("Assign taxonomy (SINTAX)" = "sintax",
                                      "Match to BINs (usearch_global)" = "usearch_global"))
           )
    )
  ),
  
  fluidRow(
    column(width = 12,
           wellPanel(
             tags$h3("3. Reference Database Files"),
             conditionalPanel(
               condition = "input.method == 'sintax'",
               shinyFilesButton("sintax_db", "📁 Select _sintax.fasta reference file", "Select SINTAX Reference", multiple = FALSE),
               tags$small(em(textOutput("sintax_path")))
             ),
             conditionalPanel(
               condition = "input.method == 'usearch_global'",
               shinyFilesButton("vsearch_db", "📁 Select _SEQUENCES_vsearch file", "Select BIN Reference", multiple = FALSE),
               tags$small(em(textOutput("vsearch_path"))),
               tags$br(), tags$br(),
               shinyFilesButton("taxonomy_file", "📁 Select _TAXONOMY.tsv file", "Select taxonomy file", multiple = FALSE),
               tags$small(em(textOutput("taxonomy_path")))
             )
           )
    )
  ),
  
  fluidRow(
    column(width = 12,
           wellPanel(
             tags$h3("4. Output Directory"),
             shinyDirButton("outdir", "📁 Select Output Folder", "Select Folder"),
             tags$small(em(textOutput("outdir_path"))),
             tags$hr(),
             actionButton("run_btn", "🚀 Run Query", class = "btn btn-primary btn-lg", width = "100%")
           )
    )
  ),
  
  fluidRow(
    column(width = 12,
           h3("Status"),
           verbatimTextOutput("status")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  shinyFileChoose(input, "fasta_file", roots = roots, session = session)
  shinyFileChoose(input, "sintax_db", roots = roots, session = session)
  shinyFileChoose(input, "vsearch_db", roots = roots, session = session)
  shinyFileChoose(input, "taxonomy_file", roots = roots, session = session)
  shinyDirChoose(input, "outdir", roots = roots, session = session)
  
  outdir <- reactive({
    req(input$outdir)
    parseDirPath(roots, input$outdir)
  })
  
  output$outdir_path <- renderText({ if(is.null(input$outdir)) "No output directory selected" else outdir() })
  
  fasta_path <- reactive({
    if (input$fasta_input_type == "file") {
      if (is.null(input$fasta_file)) return(NULL)
      parseFilePaths(roots, input$fasta_file)$datapath
    } else {
      req(input$fasta_text)
      tmp <- tempfile(fileext = ".fasta")
      writeLines(input$fasta_text, tmp)
      tmp
    }
  })
  
  sintax_path <- reactive({ if(is.null(input$sintax_db)) return(NULL); parseFilePaths(roots, input$sintax_db)$datapath })
  vsearch_path <- reactive({ if(is.null(input$vsearch_db)) return(NULL); parseFilePaths(roots, input$vsearch_db)$datapath })
  taxonomy_path <- reactive({ if(is.null(input$taxonomy_file)) return(NULL); parseFilePaths(roots, input$taxonomy_file)$datapath })
  
  output$fasta_path <- renderText({ fasta_path() })
  output$sintax_path <- renderText({ sintax_path() })
  output$vsearch_path <- renderText({ vsearch_path() })
  output$taxonomy_path <- renderText({ taxonomy_path() })
  
  observeEvent(input$run_btn, {
    req(fasta_path(), outdir())
    withProgress(message = "Processing...", value = 0, {
      incProgress(0.1, detail = "Cleaning sequences...")
      fasta <- fasta_path()
      cleaned <- tempfile(fileext = ".fasta")
      clean_fasta(fasta, cleaned)
      
      if(input$method == "sintax") {
        req(sintax_path())
        db <- sintax_path()
        incProgress(0.3, detail = "Running SINTAX...")
        out_file <- run_vsearch(cleaned, db, "sintax", threads = 6)
        incProgress(0.7, detail = "Parsing results...")
        df <- parse_sintax(out_file)
        out_name <- paste0("TAX_ASSIGN_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
      } else {
        req(vsearch_path(), taxonomy_path())
        db <- vsearch_path(); tax_file <- taxonomy_path()
        incProgress(0.3, detail = "Running usearch_global...")
        out_file <- run_vsearch(cleaned, db, "usearch_global", threads = 6)
        incProgress(0.7, detail = "Parsing results...")
        df <- parse_usearch_global(out_file, tax_file)
        out_name <- paste0("BIN_MATCH_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
      }
      
      incProgress(0.9, detail = "Saving output...")
      out_path <- file.path(outdir(), out_name)
      write.table(df, out_path, sep="\t", row.names = FALSE, quote = FALSE)
      incProgress(1, detail = "Done!")
    })
    output$status <- renderText(paste("✅ Done! Output saved to:\n", out_path))
  })
}

# ---- Run App ----
shinyApp(ui, server)

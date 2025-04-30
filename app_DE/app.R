#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(tidyverse)
library(bslib)
library(ggplot2)
library(here)
library(viridisLite)
library(gprofiler2)
library(plotly)
Masterfile <- readRDS(here::here("data/Masterfile.rda")) %>% filter(!is.na(padj)) %>% dplyr::select(-funct.class)
Masterfile_padj <- readRDS(here::here("data/Masterfile_padj.rda"))

de <- Masterfile %>% filter(!is.na(padj)) %>% subset(padj < 0.01 & abs(log2FoldChange) >= 1)
de_4hr <- subset(de, time == "4hr")
de_24hr <- subset(de, time == "24hr")

bind_with_names <- function(df, name) {
  df <- df %>% mutate(condition = name)
}

replaceNames <- function(df) {
  names(df)[names(df) == "C_Albicans"] <- "CAlb"
  names(df)[names(df) == "PolyIC"] <- "PolIC"
  names(df)[names(df) == "S_Aureus"] <- "SAur"
  return(df)
}

# Combine all data frames into one, with the list names as a column
path_4hr <- readRDS(here::here("data/pathway_res_4.rda")) %>% replaceNames()
path_4hr <- map2_df(path_4hr, names(path_4hr), bind_with_names) %>% mutate(timepoint = "4hr")
path_24hr <- readRDS(here::here("data/pathway_res_24.rda")) %>% replaceNames()
path_24hr <- map2_df(path_24hr, names(path_24hr), bind_with_names) %>% mutate(timepoint = "24hr")
path <- rbind(path_4hr, path_24hr) %>% mutate(set = ifelse(NES > 0, "up", "down"))

ui <- fluidPage( 
  page_sidebar(
    sidebar = sidebar(
      # Select variable for y-axis
      selectInput(
        inputId = "condition",
        label = "condition",
        choices = c("LPS", "CAlb", "PolIC", "SAur"),
        selected = "LPS"
      ),
      
      selectInput(
        inputId = "timepoint",
        label = "stimulation time point",
        choices = c("4hr", "24hr"),
        selected = "4hr"
      ),
      
      # sliderInput(
      #   inputId = "n_pathways",
      #   label = "Number of pathways to display",
      #   min = 1,
      #   max = 20,
      #   value = 5
      # ),
      
      sliderInput(
        inputId = "padj",
        label = "DEG Padj cutoff",
        min = 0,
        max = 0.05,
        value = 0.01
      ),
      
      sliderInput(
        inputId = "lfc",
        label = "DEG LFC cutoff",
        min = 0,
        max = 8,
        value = 1,
        step = 0.1
      ),
      downloadButton('download_degs_1', 'Download selected DEGs'),
      downloadButton('download_degs_2', 'Download all DEGs'),
      downloadButton('download_path_1', 'Download selected pathways'),
      downloadButton('download_path_2', 'Download all pathways'),
      downloadButton('download_overlap_4hr', 'Download overlapping DEGs @ 4hr'),
      downloadButton('download_overlap_24hr', 'Download overlapping DEGs @ 24hr'),
      downloadButton('download_unchar', 'Download uncharacterized DEGs')
      
    ),
    plotlyOutput(outputId = "plot"),
    DT::DTOutput("pathway_table")
  )
)

# Define server

server <- function(input, output, session) {
  observe({
    updateSelectInput(session, "timepoint", choices = unique(Masterfile$`time-point`))
    updateSelectInput(session, "condition", choices = unique(Masterfile$condition))
    
    updateSelectInput(session, "n_pathways", choices = seq(1,15))
  })
  
  filteredData <- reactive({
    req(input$timepoint)
    req(input$condition)
    Masterfile %>%
      filter(!is.na(padj)) %>% 
      mutate(padj_log10 = -log10(padj)) %>% 
      filter(`time-point` == input$timepoint, condition == input$condition) %>% 
      mutate(Significance = ifelse(padj > input$padj, "NS",
                                   ifelse(padj <= input$padj & log2FoldChange >= input$lfc, "DEG",
                                   ifelse(padj <= input$padj & log2FoldChange <= -input$lfc, "DEG", "Non-DEG"))))
  })
  
  degs <- reactive({
    req(input$timepoint)
    req(input$condition)
    Masterfile %>%
      filter(padj <= 0.01 & abs(log2FoldChange) >= 1)
  })
  
  pathwayData <- reactive({
    req(input$timepoint)
    req(input$condition)

    path %>% subset(padj < 0.01) %>% mutate(padj_log10 = -log10(padj)) %>% filter(timepoint == input$timepoint, condition == input$condition) %>%
            arrange(-NES, padj) %>% 
            mutate(pval = signif(pval, 2), padj = signif(padj, 2), NES = round(NES, 2)) %>% 
            dplyr::select(pathway, pval, padj, NES, size)
  })
  
  # output$scatterplot <- renderPlot({
  #   data <- filteredData()
  #   ggplot(data = data, aes(x = log2FoldChange, y = padj_log10, color = Significance)) +
  #     theme_bw() +
  #     geom_point()
  # })
  
  output$plot <- renderPlotly({ 
    d <- filteredData() 
    # p <-
    #   ggplot(data = d, aes(x = log2FoldChange, y = padj_log10, color = Significance)) +
    #   ggtitle(paste0(input$condition, " at ", input$timepoint)) +
    #   theme_bw() +
    #   geom_point()
    
    plot_ly(type = 'scatter', mode = 'markers', data = as_tibble(d), x = ~log2FoldChange, y = ~padj_log10, color = ~Significance, text = ~geneName,
            marker = list(size = 5)) 
  }) 
  
  output$pathway_table <- DT::renderDataTable({
    data <- pathwayData()
    data
  })
  
  output$download_degs_1 <- downloadHandler(

    filename = function() {
      paste("quantseq_degs_filtered_", input$condition, "_", input$timepoint, "_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(Masterfile %>%
                  filter(`time-point` == input$timepoint, condition == input$condition) %>%
                  mutate(Significance = ifelse(padj > input$padj, "NS",
                                               ifelse(padj <= input$padj & log2FoldChange >= input$lfc, "DEG",
                                                      ifelse(padj <= input$padj & log2FoldChange <= -input$lfc, "DEG", "Non-DEG")))),
                file, row.names = FALSE)
    })
  
  output$download_degs_2 <- downloadHandler(
    filename = function() { 
      paste("quantseq_degs_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(Masterfile, file, row.names = FALSE)
    })
  
  output$download_path_1 <- downloadHandler(
    filename = function() { 
      paste("quantseq_pathways_filtered_", input$condition, "_", input$timepoint, "_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(path %>% filter(timepoint == input$timepoint, condition == input$condition) %>%
                  arrange(-NES, padj) %>% 
                  dplyr::select(pathway, pval, padj, NES, size), 
                file, row.names = FALSE)
    })
  
  output$download_path_2 <- downloadHandler(
    filename = function() { 
      paste("quantseq_pathways_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(path, file, row.names = FALSE)
    })
  
  output$download_overlap_4hr <- downloadHandler(
    filename = function() { 
      paste("quantseq_overlap_4hr_", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(Reduce(intersect, lapply(split(de_4hr, de_4hr$condition), function(t) return(unique(t$geneName)))), file, row.names = FALSE, col.names = FALSE)
    })
  
  output$download_overlap_24hr <- downloadHandler(
    filename = function() { 
      paste("quantseq_overlap_24hr_", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(Reduce(intersect, lapply(split(de_24hr, de_24hr$condition), function(t) return(unique(t$geneName)))), file, row.names = FALSE, col.names = FALSE)
    })
  
  output$download_unchar <- downloadHandler(
    filename = function() { 
      paste("quantseq_unchar_", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(unique(subset(Masterfile_padj, uniprot_rev_score %in% c(1,2))$geneName), file, row.names = FALSE, col.names = FALSE)
    })

  
  # output$pathways <- renderPlot({
  #   data <- pathwayData()
  #   
  #   ggplot(data, aes(x = NES, y = reorder(pathway, -NES), fill = padj_log10)) +
  #     theme_bw() +
  #     geom_point() +
  #     facet_grid(set~., scales = "free_y") +
  #     labs(x = "Normalized enrichment score (NES)", y = "GO:BP pathway") +
  #     theme_bw()
  # })
}

# Create a Shiny app object
shinyApp(ui = ui, server = server)

#To deploy
#rsconnect::deployApp("/Users/emilvorsteveld/work/healthy_quantseq/main/app_DE")
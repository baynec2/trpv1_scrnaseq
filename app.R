#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(dplyr)
library(shiny)
library(shinythemes)
library(org.Hs.eg.db)

#testing 
#input = list(TRPV1_RPM_threshold = 100, n_greater_than_0 = 32, GO_term = "GO:0006955",max_scale = 100)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("darkly"),
                # Application title
                titlePanel("SCRNAseq of Dorsal Root Ganglion Neurons Expressing TRPV1 by GO Term"),
                
                # Sidebar with a slider input for number of bins 
                sidebarLayout(
                  sidebarPanel(
                    
                    textInput("GO_term", "GO Term", value = "GO:0006955", width = NULL, placeholder = NULL),
                    sliderInput("TRPV1_RPM_threshold",
                                "TRPV1 must have this many or more reads per million to include cell",
                                min = 100,
                                max = 1000,
                                value = 100),
                    sliderInput("n_greater_than_0",
                                "including only genes that have zero read counts for this number of cells or less",
                                min = 0,
                                max = 32,
                                value = 32),
                    textInput("max_scale", "max number on scale", value = "auto", width = NULL, placeholder = NULL),
                    downloadButton("download_all_data", "download all of the data"),
                    downloadButton("download_selection", "download data from heatmap")
                  ),
                  # Show a plot of the generated distribution
                  mainPanel(
                    plotly::plotlyOutput("Heatmap")
                  )
                ))


# Define server logic required to draw a histogram
server <- function(input,output){
  #Reading in data 
  data = readRDS("cleaned_data_100_RPM")
  datasetInput <- reactive({
    data
  })
    output$Heatmap <- plotly::renderPlotly({
      
      #filtering for just TRPV1
      f = data %>% 
        dplyr::filter(gene == "TRPV1")
      #defining samples we are interested in
      soi = f %>% 
        dplyr::filter(value > input$TRPV1_RPM_threshold) %>% 
        dplyr::pull(sample_id)
      
      #filtering the data to fit our criteria
      filtered = data %>% 
        dplyr::filter(sample_id %in% soi) %>% 
        dplyr::mutate(gene = toupper(gene)) %>% 
        dplyr::group_by(gene) %>% 
        dplyr::mutate(n = sum(value ==0)) %>% 
        #arbitrary filter here to cut noise
        dplyr::filter(n <= input$n_greater_than_0) %>% 
        dplyr::select(-n) %>% 
        dplyr::ungroup()
      
      #Selecting the columns that we want reported
      cols = c("SYMBOL","GENENAME")
      #Actually pulling the wanted data
      annots = AnnotationDbi::select(org.Hs.eg.db, keys = input$GO_term,columns = cols, keytype="GO") %>% 
        dplyr::select(SYMBOL,GENENAME,GO) %>% 
        dplyr::distinct()
      
      #adding the receptor we are interested in
      TRPV1 = data.frame(SYMBOL = "TRPV1", GENENAME = "transient receptor potential cation channel subfamily V member 1",GO = "GO:0034220")
      
      annots = dplyr::bind_rows(annots,TRPV1)
      
      #merging our data and annotation
      resultTable = dplyr::inner_join(filtered, annots, by = c("gene"= "SYMBOL")) 
      
      heatmap = resultTable %>% 
        #mutate(value = log10(value)) %>% 
        tidyr::pivot_wider(names_from = "sample_id",
                           values_from = value) %>% 
        dplyr::mutate(hover_info = paste0(gene,": ",GENENAME))
      
      heatmap_f = heatmap %>% 
        dplyr::select(-gene,-GENENAME,-hover_info)
      
      hover_text = matrix(rep(heatmap$hover_info,ncol(heatmap_f)),ncol = ncol(heatmap_f), byrow = FALSE)
      
      #Auto setting the color scale unless specifically asked for
      if(input$max_scale == "auto"){
      heatmaply::heatmaply(heatmap_f,labRow = heatmap$gene, custom_hovertext = hover_text) %>% 
        plotly::layout(width=1500,height = 1000)
      }else{
        heatmap_2 = heatmap_f %>% 
          dplyr::select(-GO) %>% 
          as.matrix()
        #Just setting max value to specified max 
        heatmap_2[heatmap_2 > as.numeric(input$max_scale)] = as.numeric(input$max_scale)
        
        heatmaply::heatmaply(heatmap_2,labRow = heatmap$gene, custom_hovertext = hover_text,row_side_colors = heatmap_f$GO) %>% 
        plotly::layout(width=1500,height = 1000)
      }  
    })
    # Downloadable csv of selected dataset ----
    output$download_all_data <- downloadHandler(
      filename = function() {
        paste("all_data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(data, file)
      })
    # Downloadable csv of selected dataset ----
    output$download_selection <- downloadHandler(
      filename = function() {
        paste("selected_data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(resultTable, file)
      })
}

# Run the application 
shinyApp(ui = ui, server = server)

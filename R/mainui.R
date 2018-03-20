#' Graphical User Interface for Single-Cell Visualization Analysis Toolkit
#'
#' start shiny app UI for scVAT
#'
#' @param data ...
#'
#' @return
#'
#' @author thbin
#'
#' @examples
#' 
startVATGUI<-function(data_name){
  #public data pools
  #globalObjects = ls(.GlobalEnv)
  #if(".assay.dataset" %in% globalObjects){
  #  oldData <- .GlobalEnv$.assay.dataset
  #}
  vat <- .GlobalEnv[[data_name]]
  avail.group <- c("Choose one"="",getPropMeta(vat,"Group"))
  avail.analysis.key <- c("Choose one" = "", getAnalysisKey(vat))
  diffData <- data.frame()
  diff.method <- c("t","wilcox","DESeq2")
  app <- shinyApp(
    ui <- dashboardPage(
      dashboardHeader(title = "Single-Cell Visualization Analysis Toolkit",titleWidth = 450), skin="black",
      dashboardSidebar(
        sidebarMenu(
          menuItem("tSNE Plot",tabName="tsnePage",icon=icon("spinner")),
          menuItem("Manual Clustering (tSNE)", tabName="tsneClusterPage",icon=icon("certificate")),
          menuItem("Differential Analysis", tabName="diffPage",icon=icon("compass")),
          menuItem("Analysis Visualization", tabName="analysisPage",icon=icon("spinner"))
        )
      ),
      dashboardBody(
        tabItems(
          tabItem(tabName="tsnePage",
                  fluidRow(
                    box(width=12, status="success",
                        HTML("<h4>tSNE Plot</h4>")
                    ),
                    fluidRow(
                      column(width = 8,
                             plotlyOutput("tsnePlot",height=600)
                      ),
                      column(width = 4,
                             selectInput("tsneGroup","Group:",avail.group,selected=avail.group[2])
                             #,selectInput("tsneShape","Shape:",avail.group,selected=avail.group[1])
                             #,sliderInput("tsne", "Number of observations:", 1, 100, 50)
                      )
                    )
                  )
          ),
          tabItem(tabName="tsneClusterPage",
                   fluidRow(
                     box(width = 12, status="success",
                         HTML("<h4>Manual Clustering the cells manually using tSNE Plot</h4>")
                     ),
                     fluidRow(
                       column(width = 8,
                              plotlyOutput("tsneClusterPlot",height=600)
                       ),
                       column(width = 4,
                              selectInput("tsneClusterGroup","Group:",avail.group,selected="0")
                              ,selectInput("tsneClusterShape","Shape:",avail.group,selected="0")
                              #,sliderInput("tsne", "Number of observations:", 1, 100, 50)
                              ,numericInput("tsneSetGroup","Group Id:","0")
                              ,wellPanel(
                                actionButton("saveTsne","Save")
                                ,actionButton("refreshTsne","Refresh")
                                ,actionButton("resetTsne","Reset")
                              )
                       )
                     )
                   )
                   
          ),
          tabItem(tabName="diffPage",
                  fluidRow(
                    box(width = 12, status="success",
                        HTML("<h4>Differentializing Analysis</h4>")
                    ),
                    fluidRow(
                      box(width=12,
                        column(width=10,
                          div(style = "display:inline-block",selectInput("diffGroup"," Cluster: ",avail.group,selected=avail.group[2]))
                          ,div(style = "display:inline-block",selectInput("diffMethod"," Method: ",diff.method,selected=diff.method[1]))
                          ,div(style = "display:inline-block",uiOutput("chooseDiffGroup1"))
                          ,div(style = "display:inline-block",uiOutput("chooseDiffGroup2"))
                          ,div(style = "display:inline-block",numericInput("diffLogFC","Min. LogFC:",0.25,min=0, max=1,step=0.01))
                          ,div(style = "display:inline-block",numericInput("diffAvg","Min. AVG:",0.1,min=0, max=1,step=0.01))
                          ,div(style = "display:inline-block",checkboxInput("diffOnlyPos","Only Postive",FALSE))
                        ),
                        column(width=2,
                               actionButton("diffDo","Calculate",icon=icon("play")),
                               downloadButton("diffDownload","Download")
                               )
                      ),
                      wellPanel(width=12,
                        DTOutput("diffDataTable", width="100%")  
                      )
                    )
                  )

          ),
          tabItem(tabName="analysisPage",
                  fluidRow(
                    box(width=12, status="success",
                        HTML("<h4>More Analysis Data Visualization</h4>")
                    ),
                    fluidRow(
                      column(width = 8,
                             plotlyOutput("analysisPlot",height=600)
                      ),
                      column(width = 4,
                             uiOutput("chooseAnalysisKey"),
                             uiOutput("choosePlotDims"),
                             textInput("analysisGene1","Gene1:",""),
                             textInput("analysisGene2","Gene2:",""),
                             textInput("analysisColor","Color:","lightgrey,blue,orange,red"),
                             checkboxInput("analysisGradient", "Gradient", TRUE)
                             ,checkboxInput("analysisMultiPlot","Plot Multiple Genes", FALSE)
                             ,numericInput("analysisPlotRows","Number of Rows:",2, min = 1,step = 1)
                      )
                    )
                  )
          )
        )
      )
    ),
    server <- function(input, output, session){
      
      notify.id <- NULL
      
      data <- reactive(
        eval(parse(text=data_name))
      )
      
      output$tsnePlot <- renderPlotly({
        input$refreshTsne
        plotTSNE(vat,group.id = input$tsneGroup, source="tsneP")
      })
      
      #Start Manual Clustering (tSNE)
      output$tsneClusterPlot <- renderPlotly({
        input$refreshTsne
        plotTSNE(vat, group.id="manual.cluster", gradient = TRUE, source="tsneCP")
      })
      observeEvent(input$saveTsne,{
        selected_data <- event_data("plotly_selected",source="tsneCP")
        if(length(selected_data)>0){
          val <- input$tsneSetGroup
          select_id <- selected_data$pointNumber + 1
          vat@cell.props$manual.cluster[select_id] <<- as.numeric(val)
        }
      })
      observeEvent(input$resetTsne,{
        vat@cell.props$manual.cluster <<- 0
      })
      #End Manual Clustering (tSNE)
      
      #Start Analysis Visualization
      output$analysisPlot <- renderPlotly({
        key <- input$analysisKey
        if(isEmpty(key)) {
          return ()
        }
        dims <- input$plotDims
        if(is.null(dims)) return()
        if(length(dims)<2 || length(dims)>3){
          return()
        }
        
        gene1 <- input$analysisGene1
        gene2 <- input$analysisGene2
        colors <- input$analysisColor
        colors <- strsplit(colors,split=",")[[1]]
        size <- NULL
        sizes <- NULL
        if(length(dims)==3){
          size <- 1
          sizes <- c(1:5)
        }
        if(isEmpty(gene1)){
          plotAnalysis(vat,key = key, dims=dims, source="analysisP",size=size, sizes=sizes)
        }else{
          gene1 <- strsplit(gene1,split=",")[[1]]
          if(input$analysisMultiPlot && length(dims)==2){#3d plot.ly doesn't work!
            plotMultipleGenes(vat, gene1, nrows = input$analysisPlotRows, dims=dims, 
                              key = key, gradient = input$analysisGradient, 
                              colors = colors[1:2],size=size, sizes=sizes)
          }else{
            if(isEmpty(gene2)){
              plotGene(vat, gene1, dims=dims, 
                       key = key, gradient = input$analysisGradient,
                       colors = colors[1:2],size=size, sizes=sizes)
            }else{
              plotTwoGenes(vat, gene1, gene2, dims=dims, 
                           key = key,colors=colors,size=size, sizes=sizes)
            }
          }
        }
      })
      output$chooseAnalysisKey <- renderUI({
        selectInput("analysisKey","Analysis Data:",avail.analysis.key,selected=avail.analysis.key[1])
      })
      output$choosePlotDims <- renderUI({
        key <- input$analysisKey
        if(isEmpty(key)) return()
        dims <- getAnalysisColnames(vat = vat,key=key)
        selectInput("plotDims", "Plot Dims",dims, multiple=TRUE)
      })
      #End Analysis Visualization
      
      #Start Defferential Analysis
      output$chooseDiffGroup1 <- renderUI({
        group.key <- input$diffGroup
        all.groups <- sort(unique(vat@cell.props[,group.key]))
        all.groups <- c(as.character(all.groups),"All")
        selectInput("diffGroup1"," Group 1: ",all.groups,selected=all.groups[1],multiple = TRUE)
      })
      output$chooseDiffGroup2 <- renderUI({
        group.key <- input$diffGroup
        all.groups <- sort(unique(vat@cell.props[,group.key]))
        all.groups <- c("Choose one"="",as.character(all.groups))
        selectInput("diffGroup2"," Group 2: ",all.groups,multiple = TRUE)
      })
      observeEvent(input$diffDo,{
        notify.id <<- showNotification("Calculating differences, Please waiting...",duration=NULL)
        group.key <- input$diffGroup
        group1 <- input$diffGroup1
        group2 <- input$diffGroup2
        if(isEmpty(group2)) group2 <- NULL
        method <- input$diffMethod
        diffLogFC <- input$diffLogFC
        diffAvg <- input$diffAvg
        only.pos <- input$diffOnlyPos
        if("All" %in% group1){
          diffData <<- tryCatch({
            doAllDiffAnalysis(
              vat, group.key = group.key,method = method,
              min.logfc = diffLogFC, min.avg = diffAvg, only.pos = only.pos
            )
          }, error = function(cond){
            removeNotification(notify.id)
            showNotification(sprintf("Run Error: %s",cond),type="error", duration=10)
            return(NULL)
          }
          )
        }else{
          diffData <<- tryCatch({
            doDiffAnalysis(
              vat, group1 = group1, group2 = group2, group.key = group.key,
              method = method,
              min.logfc = diffLogFC, min.avg = diffAvg, only.pos = only.pos
            )
          }, error = function(cond){
            removeNotification(notify.id)
            showNotification(sprintf("Run Error: %s",cond),type="error", duration=10)
            return(NULL)
          }
          )
        }
      })
      output$diffDownload <- downloadHandler(
        filename = function(){
          diffLogFC <- input$diffLogFC
          paste0("diff_", diffLogFC, ".csv")
        },
        content = function(file){
          write.csv(diffData, file, row.names = FALSE)
        }
      )
      output$diffDataTable <- renderDT({
        input$diffDo
        if (!is.null(notify.id)) removeNotification(notify.id)
        notify.id <<- NULL
        diffData
      }#,options=list(searching=FALSE,pageLength = 20)
      )
      #End Defferential Analysis
      
      session$onSessionEnded(function(){
        .GlobalEnv[[data_name]] <- vat
        stopApp()
        })
    }
  )
  #shinyApp(ui, server, options = list(launch.browser=TRUE))
  runApp(app, launch.browser = TRUE)
}

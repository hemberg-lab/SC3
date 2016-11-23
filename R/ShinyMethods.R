#' Opens \code{SC3} results in an interactive session in a web browser.
#'
#' Runs interactive \code{shiny} session of \code{SC3} based on precomputed clusterings.
#'
#' @param object an object of \code{SCESet} class
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize
#' all precomputed clusterings.
#' 
#' @name sc3_interactive
#' @aliases sc3_interactive, sc3_interactive,SCESet-method
#'
#' @importFrom shiny HTML actionButton animationOptions checkboxGroupInput column div downloadHandler downloadLink eventReactive fluidPage fluidRow h4 headerPanel htmlOutput need observe observeEvent p plotOutput reactiveValues renderPlot renderUI selectInput shinyApp sliderInput stopApp tabPanel tabsetPanel uiOutput updateSelectInput validate wellPanel withProgress conditionalPanel reactive outputOptions tags radioButtons downloadButton
#' @importFrom utils head
#' @importFrom stats median
#' @importFrom graphics plot
#' @export
#'
sc3_interactive.SCESet <- function(object) {
    consensus <- object@sc3$consensus
    if (is.null(consensus)) {
        warning(paste0("Please run sc3_calc_consens() first!"))
        return()
    }
    
    ks <- as.numeric(names(consensus))
    dataset <- get_processed_dataset(object)
    
    ## define UI parameters
    plot.height <- 600
    plot.height.small <- 300
    
    ## define server global variables
    values <- reactiveValues()
    
    svm <- FALSE
    if (!is.null(object@sc3$svm_train_inds)) {
        svm <- TRUE
        dataset <- dataset[,object@sc3$svm_train_inds]
    }
    
    biology <- FALSE
    if (!is.null(object@sc3$biology)) {
        if(object@sc3$biology) {
            biology <- TRUE
        }
    }
    
    pdata <- colnames(make_col_ann_for_heatmaps(object, colnames(object@phenoData@data)))
    
    shinyApp(
        ui = fluidPage(
            tags$head(tags$style(
                HTML(".shiny-output-error-validation {
                     color: red;
}")
                )),
            fluidRow(column(12,
                            HTML("<br>"))),
            fluidRow(
                column(
                    3,
                    h4("Parameters"),
                    wellPanel(fluidRow(
                        column(
                            12,
                            conditionalPanel(
                                "output.ks_length!='1'",
                                sliderInput(
                                    "clusters",
                                    label = "Number of clusters k",
                                    min = min(ks),
                                    max = max(ks),
                                    value = stats::median(ks),
                                    step = 1,
                                    animate = animationOptions(interval = 2000,
                                                               loop = FALSE)
                                )
                            )
                        ),
                        column(12,
                               conditionalPanel(
                                   "output.ks_length=='1'",
                                   HTML(paste0("<h4>k = ", unique(ks), "<h4>"))
                               ))
                    ),
                    fluidRow(column(
                        12,
                        conditionalPanel("output.is_svm",
                                         p(HTML(
                                             paste0(
                                                 "<font color = 'red'>Your data was clustered based on ",
                                                 length(object@sc3$svm_train_inds),
                                                 " selected cells.</font>"
                                             )
                                         )))
                    ))),
                    
                    conditionalPanel(
                        "input.main_panel == 'Consensus' || input.main_panel == 'Expression' || input.main_panel == 'DE' || input.main_panel == 'Markers'",
                    wellPanel(fluidRow(
                        column(12,
                               HTML("<b>Show phenoData</b>"),
                               HTML("<font size = 1>"),
                               checkboxGroupInput("pDataColumns", label = NULL, pdata, selected = NULL),
                               HTML("</font>")
                        ))))
                ),
                column(6,
                       uiOutput('mytabs')),
                column(3,
                       fluidRow(
                           column(
                               12,
                               h4("Panel description"),
                               htmlOutput("explanation"),
                               conditionalPanel(
                                   "input.main_panel == 'Markers' && output.is_biology",
                                   wellPanel(
                                       sliderInput(
                                           "auroc.threshold",
                                           label = "AUROC threshold",
                                           min = 0.6,
                                           max = 0.95,
                                           value = 0.85,
                                           step = 0.05
                                       ),
                                       radioButtons(
                                           "p.val.mark",
                                           label = "p-value threshold",
                                           choices = c(
                                               "0.01" = 0.01,
                                               "0.05" = 0.05,
                                               "0.1" = 0.1
                                           ),
                                           selected = "0.01",
                                           inline = TRUE
                                       )
                                   )
                               ),
                               conditionalPanel(
                                   "input.main_panel == 'DE' && output.is_biology",
                                   wellPanel(
                                       radioButtons(
                                           "p.val.de",
                                           label = "p-value threshold",
                                           choices = c(
                                               "0.01" = 0.01,
                                               "0.05" = 0.05,
                                               "0.1" = 0.1
                                           ),
                                           selected = "0.01",
                                           inline = TRUE
                                       )
                                   )
                               )
                           )
                       ))
            )
                ),
        server = function(input, output, session) {
            # render tabpanel
            output$mytabs = renderUI({
                myTabs <- list(
                    tabPanel(
                        "Consensus",
                        plotOutput('consensus',
                                   height = plot.height,
                                   width = "100%")
                    ),
                    tabPanel(
                        "Silhouette",
                        plotOutput('silh',
                                   height = plot.height,
                                   width = "100%")
                    ),
                    tabPanel(
                        "Stability",
                        plotOutput(
                            'StabilityPlot',
                            height = plot.height.small,
                            width = "100%"
                        )
                    ),
                    tabPanel(
                        "Expression",
                        plotOutput('matrix',
                                   height = plot.height,
                                   width = "100%")
                    ),
                    tabPanel("DE",
                             uiOutput('plot_de_genes')),
                    tabPanel("Markers",
                             uiOutput('plot_mark_genes'))
                )
                do.call(tabsetPanel,
                        c(myTabs,
                          id = "main_panel",
                          type = "pills"))
            })
            # observer for marker genes
            observe({
                if (biology) {
                    # get all marker genes
                    markers <- organise_marker_genes(object, input$clusters, as.numeric(input$p.val.mark), as.numeric(input$auroc.threshold))
                    values$n.markers <- nrow(markers)
                    # get top 10 marker genes of each cluster
                    markers <- markers_for_heatmap(markers)
                    clusts <- unique(markers[,1])
                    if (is.null(clusts))
                        clusts <- "None"
                    values$mark.res <- markers
                    updateSelectInput(session, "cluster", choices = clusts)
                } else {
                    values$n.markers <- 0
                }
            })
            # observer for DE genes
            observe({
                if (biology) {
                    de_genes <- organise_de_genes(object, input$clusters, as.numeric(input$p.val.de))
                    values$n.de.genes <- length(de_genes)
                    values$n.de.plot.height <- length(head(de_genes, 50))
                } else {
                    values$n.de.genes <- 0
                }
            })
            ## REACTIVE PANEL DESCRIPTIONS
            output$explanation <- renderUI({
                res <- ""
                if (length(input$main_panel) > 0) {
                    if (grepl("Consensus", input$main_panel)) {
                        res <- HTML(
                            "The consensus matrix is a <i>N</i>x<i>N</i>
                            matrix, where <i>N</i> is the number of cells.
                            It represents similarity between the cells based
                            on the averaging of clustering results from all
                            combinations of clustering parameters. Similarity 0
                            (<font color = 'blue'>blue</font>) means that
                            the two cells are always assigned to different clusters.
                            In contrast, similarity 1 (<font color = 'red'>red</font>)
                            means that the two cells are always assigned
                            to the same cluster. The consensus matrix is
                            clustered by hierarchical clustering and has
                            a diagonal-block structure. Intuitively, the
                            perfect clustering is achieved when all diagonal
                            blocks are completely <font color = 'red'>red</font>
                            and all off-diagonal elements are completely <font color = 'blue'>blue</font>.
                            For a more objective clustering validation please visit the <b>Silhouette</b> panel."
                        )
                    }
                    if (grepl("Silhouette", input$main_panel)) {
                        res <- HTML(
                            "As opposed to visual exploration of the
                            consensus matrix, a silhouette is a quantitative
                            measure of the diagonality of the consensus
                            matrix. An average silhouette width
                            (shown at the bottom left of the silhouette
                            plot) varies from 0 to 1, where 1 represents
                            a perfectly block-diagonal consensus matrix
                            and 0 represents a situation where there is
                            no block-diagonal structure. The best
                            clustering is achieved when the average
                            silhouette width is close to 1."
                        )
                    }
                    if (grepl("Stability", input$main_panel)) {
                        res <- HTML(
                            "Stability index shows how stable each cluster
                            is accross the selected range of <i>k</i>s.
                            The stability index varies between 0 and 1, where
                            1 means that the same cluster appears in every
                            solution for different <i>k</i>."
                        )
                    }
                    if (grepl("Expression", input$main_panel)) {
                        res <- HTML(
                            "The expression panel represents the original
                            input expression matrix (cells in columns and
                            genes in rows) after cell and gene filters.
                            Genes are clustered by <i>k</i>-means with <i>k</i> = 100
                            (dendrogram on the left) and the heatmap
                            represents the expression levels of the gene
                            cluster centers after <i>log2</i>-scaling."
                        )
                    }
                    if (grepl("DE", input$main_panel)) {
                        res <-
                            HTML(
                                paste0(
                                    "SC3 found <b>",
                                    values$n.de.genes,
                                    "</b> differentially expressed genes based
                                    on the obtained clustering.<br>",
                                    "Differential expression is calculated using
                                    the non-parametric Kruskal-Wallis test.
                                    A significant <i>p</i>-value indicates that gene
                                    expression in at least one cluster
                                    stochastically dominates one other cluster.
                                    SC3 provides a list of all differentially
                                    expressed genes with adjusted <i>p</i>-values < 0.01
                                    and plots gene expression profiles of the
                                    50 genes with the lowest <i>p</i>-values. Note
                                    that the calculation of differential
                                    expression after clustering can introduce
                                    a bias in the distribution of <i>p</i>-values,
                                    and thus we advise to use the <i>p</i>-values for
                                    ranking the genes only.<br>The <i>p</i>-value threshold
                                    can be changed using the radio buttons below.<br><br>"
                                )
                                )
                    }
                    if (grepl("Markers", input$main_panel)) {
                        res <-
                            HTML(
                                paste0(
                                    "SC3 found <b>",
                                    values$n.markers,
                                    "</b> marker genes based
                                    on the obtained clustering.<br>",
                                    "To find marker genes, for each gene a binary classifier is constructed
                                    based on the mean cluster expression values.
                                    The classifier prediction is then calculated
                                    using the gene expression ranks. The area
                                    under the receiver operating characteristic
                                    (ROC) curve is used to quantify the accuracy
                                    of the prediction. A <i>p</i>-value is assigned to
                                    each gene by using the Wilcoxon signed rank
                                    test. The genes with the area under the ROC
                                    curve (AUROC) > 0.85 and with the <i>p</i>-value
                                    < 0.01 are defined as marker genes and the
                                    top 10 marker genes of each cluster are
                                    visualized in this panel. The AUROC and the
                                    p-value thresholds can be changed using the
                                    slider and radio buttons below.<br><br>"
                                )
                                )
                    }
                    }
                return(res)
                    })
            ## REACTIVE PANELS
            # plot consensus matrix
            output$consensus <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    sc3_plot_consensus(object, as.numeric(input$clusters), show_pdata = input$pDataColumns)
                })
            })
            # plot silhouette
            output$silh <- renderPlot({
                sc3_plot_silhouette(object, as.numeric(input$clusters))
            })
            # plot stability
            output$StabilityPlot <- renderPlot({
                validate(need(
                    length(object@sc3$consensus) > 1,
                    "\nStability cannot be calculated for a single k value!"
                ))
                sc3_plot_cluster_stability(object, input$clusters)
            })
            # plot expression matrix
            output$matrix <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    set.seed(1234567)
                    sc3_plot_expression(object, as.numeric(input$clusters), show_pdata = input$pDataColumns)
                })
            })
            # plot marker genes
            output$mark_genes <- renderPlot({
                validate(need(
                    biology,
                    "\nPlease run sc3_calc_biology() first!"
                ))
                validate(
                    need(
                        !is.null(values$mark.res),
                        "\nThere were no marker genes found! Try change the parameters."
                    )
                )
                validate(need(
                    try(!is.null(rownames(dataset)))
                    ,
                    "\nNo gene names provided in the input expression matrix!"
                ))
                sc3_plot_markers(
                    object,
                    as.numeric(input$clusters),
                    as.numeric(input$auroc.threshold),
                    as.numeric(input$p.val.mark),
                    show_pdata = input$pDataColumns
                )
            })
            plotHeightMark <- function() {
                return(150 + 10.8 * nrow(values$mark.res))
            }
            output$plot_mark_genes <- renderUI({
                validate(need(
                    biology,
                    "\nPlease run sc3_calc_biology() first!"
                ))
                validate(
                    need(
                        !is.null(values$mark.res),
                        "\nThere were no marker genes found! Try change the parameters."
                    )
                )
                validate(need(
                    try(!is.null(rownames(dataset)))
                    ,
                    "\nNo gene names provided in the input expression matrix!"
                ))
                plotOutput("mark_genes",
                           height = plotHeightMark(),
                           width = "100%")
            })
            # plot DE genes
            output$de_genes <- renderPlot({
                validate(need(
                    biology,
                    "\nPlease run sc3_calc_biology() first!"
                ))
                validate(
                    need(
                        values$n.de.genes > 0,
                        "\nThere were no DE genes found! Try change the parameters."
                    )
                )
                validate(need(
                    try(!is.null(rownames(dataset)))
                    ,
                    "\nNo gene names provided in the input expression matrix!"
                ))
                sc3_plot_de_genes(object,
                                  as.numeric(input$clusters),
                                  as.numeric(input$p.val.de),
                                  show_pdata = input$pDataColumns)
            })
            plotHeightDE <- function() {
                return(150 + 10.8 * values$n.de.plot.height)
            }
            output$plot_de_genes <- renderUI({
                validate(need(
                    biology,
                    "\nPlease run sc3_calc_biology() first!"
                ))
                validate(
                    need(
                        values$n.de.genes > 0,
                        "\nThere were no DE genes found! Try change the parameters."
                    )
                )
                validate(need(
                    try(!is.null(rownames(dataset)))
                    ,
                    "\nNo gene names provided in the input expression matrix!"
                ))
                plotOutput("de_genes",
                           height = plotHeightDE(),
                           width = "100%")
            })
            
            # REACTIVE BUTTONS
            
            ## OUTPUTS USED FOR CONDITIONAL PANELS
            ks_length <- reactive({
                return(as.character(length(ks)))
            })
            
            output$ks_length <- reactive({
                return(ks_length())
            })
            
            is_svm <- reactive({
                return(svm)
            })
            
            output$is_svm <- reactive({
                return(is_svm())
            })
            
            is_biology <- reactive({
                return(biology)
            })
            
            output$is_biology <- reactive({
                return(is_biology())
            })
            
            # stop App on closing the browser
            session$onSessionEnded(function() {
                stopApp()
            })
            outputOptions(output, 'ks_length', suspendWhenHidden = FALSE)
            outputOptions(output, 'is_svm', suspendWhenHidden = FALSE)
            outputOptions(output, 'is_biology', suspendWhenHidden = FALSE)
                    },
        
        # launch App in a browser
        options = list(launch.browser = TRUE)
                        )
                    }

#' @rdname sc3_interactive
#' @aliases sc3_interactive
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_interactive", signature(object = "SCESet"),
          function(object) {
              sc3_interactive.SCESet(object)
          })

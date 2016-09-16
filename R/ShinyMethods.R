#' Open SC3 results in an interactive session in a web browser
#'
#' Runs interactive session of SC3 based on precomputed objects
#'
#' @param object an object of "SCESet" class
#'
#' @return Opens a browser window with an interactive shiny app and visualize
#' all precomputed clusterings.
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
    dataset <- object@sc3$processed_dataset
    
    ## define UI parameters
    plot.height <- 600
    plot.height.small <- 300
    
    ## define server global variables
    values <- reactiveValues()
    
    svm <- FALSE
    if (!is.null(object@sc3$svm_train_inds)) {
        svm <- TRUE
    }
    
    biology <- FALSE
    if (!is.null(object@sc3$biology)) {
        biology <- TRUE
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
                                                 "<font color = 'red'>Your data was clustered in the SVM regime based on ",
                                                 length(object@sc3$svm_train_inds),
                                                 " cells. When you have found the best number of clusters <b><em>k</em></b>, go back to your terminal session and run sc3_run_svm() to predict the labels of the other cells.</font>"
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
                        )))),
                    
                    fluidRow(
                        column(
                            12,
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
                            ),
                            conditionalPanel("input.main_panel == 'tSNE'",
                                             wellPanel(
                                                 sliderInput(
                                                     "perplexity",
                                                     label = "Perplexity",
                                                     min = floor(0.7 * ncol(dataset) / 5),
                                                     max = floor(1.3 * ncol(dataset) / 5),
                                                     value = floor(ncol(dataset) / 5)
                                                 )
                                             ))
                        )
                    )
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
                                       selectInput(
                                           "cluster",
                                           "Choose a cluster for GO analysis:",
                                           c("None" = "NULL")
                                       ),
                                       if (object@sc3$rselenium) {
                                           actionButton("GO", label = "Analyze!")
                                       },
                                       if (object@sc3$rselenium) {
                                           HTML(
                                               "<br>(opens <a href = 'http://biit.cs.ut.ee/gprofiler' target='_blank'>g:Profiler</a> in Firefox)"
                                           )
                                       } else {
                                           HTML(
                                               "<font color='red'>To be able to run GO analysis you need to install RSelenium library. You can do that by closing this window and then running 'RSelenium::checkForServer()' command in your R session. This will download the required library. After that please rerun SC3 again. More details are available <a href = 'https://cran.r-project.org/web/packages/RSelenium/vignettes/RSelenium-basics.html' target='_blank'>here</a></font>."
                                           )
                                       }
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
                        "Labels",
                        div(htmlOutput('labels'),
                            style = "font-size:80%")
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
                    tabPanel(
                        "tSNE",
                        plotOutput('tSNEplot',
                                   height = plot.height,
                                   width = "100%")
                    ),
                    tabPanel("DE",
                             uiOutput('plot_de_genes')),
                    tabPanel("Markers",
                             uiOutput('plot_mark_genes')),
                    tabPanel(
                        "Outliers",
                        plotOutput('outliers',
                                   height = plot.height.small,
                                   width = "100%")
                    )
                )
                do.call(tabsetPanel,
                        c(myTabs,
                          id = "main_panel",
                          type = "pills"))
            })
            # observer for labels
            observe({
                res <- prepare_output(object, input$clusters)
                values$labels <- res$labels
                values$labels1 <- res$labels1
            })
            # observer for marker genes
            observe({
                if (biology) {
                    markers <-
                        object@sc3$biology[[as.character(input$clusters)]]$markers
                    markers <-
                        markers[markers$AUC >= as.numeric(input$auroc.threshold) &
                                    markers$p.value < as.numeric(input$p.val.mark),]
                    values$n.markers <- nrow(markers)
                    mark.res.plot <-
                        mark_gene_heatmap_param(markers)
                    clusts <- unique(mark.res.plot$sc3_clusters)
                    if (is.null(clusts))
                        clusts <- "None"
                    values$mark.res <- mark.res.plot
                    updateSelectInput(session, "cluster", choices = clusts)
                } else {
                    values$n.markers <- 0
                }
            })
            # observer for DE genes
            observe({
                if (biology) {
                    de.genes <-
                        object@sc3$biology[[as.character(input$clusters)]]$de.genes
                    de.genes <-
                        de.genes[de.genes$p.value < as.numeric(input$p.val.de), , drop = FALSE]
                    values$n.de.genes <- nrow(de.genes)
                    values$n.de.plot.height <-
                        nrow(head(de.genes, 50))
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
                    if (grepl("Labels", input$main_panel)) {
                        res <- HTML(
                            "The labels panel shows how cells are distributed
                            in the clusters that were obtained using the
                            original cell indexes from the expression matrix.
                            To help visualize how the clustering changes
                            when changing from <i>k</i> - 1 to <i>k</i>, the labels are
                            colour-coded corresponding to the custering
                            results for <i>k</i> - 1."
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
                    if (grepl("tSNE", input$main_panel)) {
                        res <-
                            HTML(
                                "<a href = 'https://lvdmaaten.github.io/tsne/' target='_blank'>tSNE</a> (t-Distributed
                                Stochastic Neighbor Embedding) method is used to
                                map high-dimensional data to a 2D
                                space while preserving local distances between
                                cells. tSNE has become a very popular visualisation
                                tool. SC3 imports the Rtsne function from the
                                <a href='https://cran.r-project.org/web/packages/Rtsne/index.html' target='_blank'>
                                Rtsne package</a> to perform the tSNE analysis.
                                The colors on the plot correspond to the clusters
                                identified by SC3.<br>
                                One of the most sensitive parameters in tSNE analysis is the
                                so-called <i>perplexity</i>. SC3 defines the default
                                <i>perplexity</i> as <i>N</i>/5, where <i>N</i> is
                                the number of cells, but also allows to change it
                                in a reasonable interval using the slider in the <b>Parameters</b> panel."
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
                                    can be changed using the radio buttons in the <b>Parameters</b> panel."
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
                                    slider and radio buttons in the <b>Parameters</b> panel.<br>
                                    In addition, SC3 allows you to run a Gene
                                    Ontology (GO) and Pathway Enrichment analysis using a panel below.<br><br>"
                                )
                                )
                    }
                    if (grepl("Outliers", input$main_panel)) {
                        res <- HTML(
                            "Outlier cells in each cluster are detected
                            using robust distances, calculated using
                            the minimum covariance determinant (MCD).
                            The outlier score shows how different a
                            cell is from all other cells in the cluster
                            and it is defined as the differences between
                            the square root of the robust distance and
                            the square root of the 99.99% quantile of
                            the Chi-squared distribution.<br><br>"
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
            # output cell labels
            output$labels <- renderUI({
                labs <- "<br/><br/>"
                if (!is.null(values$labels1)) {
                    labs1 <- list()
                    cols <- iwanthue(input$clusters - 1)
                    for (i in 1:(input$clusters - 1)) {
                        col <- cols[i]
                        ind <-
                            unlist(strsplit(as.character(values$labels1[i,]),
                                            " "))
                        for (j in ind) {
                            labs1[[j]] <-
                                paste0("<font color=\"",
                                       col,
                                       "\">",
                                       j,
                                       "</font>")
                        }
                    }
                    for (i in 1:input$clusters) {
                        ind <- unlist(strsplit(as.character(values$labels[i,]),
                                               " "))
                        for (j in ind) {
                            labs <- c(labs, labs1[[j]])
                        }
                        labs <- c(labs, c("<br/>", "<hr>"))
                    }
                } else {
                    for (i in 1:input$clusters) {
                        ind <- unlist(strsplit(as.character(values$labels[i,]),
                                               " "))
                        labs <- c(labs, ind, c("<br/>", "<hr>"))
                    }
                }
                
                HTML(paste0(labs))
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
            # make tSNE plot
            output$tSNEplot <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    sc3_plot_tsne(object,
                                  input$clusters,
                                  input$perplexity)
                })
            })
            # plot marker genes
            output$mark_genes <- renderPlot({
                validate(need(
                    !is.null(object@sc3$biology),
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
                    !is.null(object@sc3$biology),
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
                    !is.null(object@sc3$biology),
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
                    !is.null(object@sc3$biology),
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
            # plot outliers
            output$outliers <- renderPlot({
                validate(need(
                    !is.null(object@sc3$biology),
                    "\nPlease run sc3_calc_biology() first!"
                ))
                sc3_plot_cell_outliers(object,
                                       as.numeric(input$clusters))
            })
            
            # REACTIVE BUTTONS
            
            # GO button
            observeEvent(input$GO, {
                if (!is.null(values$mark.res)) {
                    open_gprofiler(rownames(values$mark.res[values$mark.res$clusts == as.numeric(input$cluster),]))
                }
            })
            
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

#' @rdname sc3_interactive.SCESet
#' @aliases sc3_interactive
#' @importClassesFrom scater SCESet
#' @export
setMethod("sc3_interactive", signature(object = "SCESet"),
          function(object) {
              sc3_interactive.SCESet(object)
          })

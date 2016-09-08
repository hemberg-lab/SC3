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
    
    consensus <- object@consensus$sc3_consensus
    if ( is.null(consensus) ) {
        warning(paste0("Please run sc3_calc_consens() first!"))
        return()
    }
    
    ks <- as.numeric(names(consensus))
    dataset <- object@consensus$sc3_processed_dataset
    show.original.labels <- FALSE
    
    ## define UI parameters
    plot.height <- 600
    plot.height.small <- 300
    
    ## define server global variables
    values <- reactiveValues()
    
    svm <- FALSE
    if(!is.null(object@consensus$svm_train_inds)) {
        svm <- TRUE
    }
    
    shinyApp(
        ui = fluidPage(
            tags$head(
                tags$style(
                    HTML(
                        ".shiny-output-error-validation {
                            color: red;
                        }"
                    )
                )
            ),
            fluidRow(
                column(
                    12,
                    HTML(
                        "<br>"
                    )
                )
            ),
            fluidRow(
                column(
                    3,
                    h4("Parameters"),
                    wellPanel(
                    fluidRow(
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
                                    animate = animationOptions(
                                        interval = 2000,
                                        loop = FALSE
                                    )
                                )
                            )
                        ),
                        column(
                            12,
                            conditionalPanel(
                                "output.ks_length=='1'",
                                HTML(paste0("<h4>k = ", unique(ks), "<h4>")
                                )
                            )
                        )
                    ),
                    fluidRow(
                        column(
                            12,
                            conditionalPanel(
                                "output.is_svm",
                                h4("SVM prediction"),
                                p(
                                    paste0(
                                        "Your data was clustered based on ",
                                        length(object@consensus$svm_train_inds),
                                        " cells. When you have found the best clustering parameters, go back to your terminal session and run XXX to predict the labels of the other cells."
                                    )
                                )
                            )
                        )
                    )),
                    
                    fluidRow(
                        column(
                            12,
                            conditionalPanel(
                                "input.main_panel == 'Markers'",
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
                                    ),
                                    selectInput(
                                        "cluster",
                                        "Choose a cluster for GO analysis:",
                                        c("None" = "NULL")
                                    ),
                                    if(object@consensus$rselenium) {
                                        actionButton("GO", label = "Analyze!")
                                    },
                                    if(object@consensus$rselenium) {
                                        HTML("<br>(opens <a href = 'http://bioinfo.vanderbilt.edu/webgestalt/' target='_blank'>WebGestalt</a> in Firefox)")
                                    } else {
                                        HTML("<font color='red'>To be able to run GO analysis you need to install RSelenium library. You can do that by closing this window and then running 'RSelenium::checkForServer()' command in your R session. This will download the required library. After that please rerun SC3 again. More details are available <a href = 'https://cran.r-project.org/web/packages/RSelenium/vignettes/RSelenium-basics.html' target='_blank'>here</a></font>.")
                                    }
                                )
                            ),
                            conditionalPanel(
                                "input.main_panel == 'DE'",
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
                            conditionalPanel(
                                "input.main_panel == 'tSNE'",
                                wellPanel(
                                    sliderInput(
                                        "perplexity",
                                        label = "Perplexity",
                                        min = floor(0.7 * ncol(dataset) / 5),
                                        max = floor(1.3 * ncol(dataset) / 5),
                                        value = floor(ncol(dataset) / 5)
                                    )
                                )
                            )
                        )
                    )
                ),
                column(
                    6,
                    uiOutput('mytabs')
                ),
                column(
                    3,
                    fluidRow(
                        column(
                            12,
                            h4("Panel description"),
                            htmlOutput("explanation")
                        )
                    )
                )
            )
        ),
        
        server = function(input, output, session) {
            # render tabpanel
            output$mytabs = renderUI({
                myTabs <- list(
                    tabPanel(
                        "Consensus",
                        plotOutput(
                            'consensus',
                            height = plot.height,
                            width = "100%"
                        )
                    ),
                    tabPanel(
                        "Silhouette",
                        plotOutput(
                            'silh', 
                            height = plot.height, 
                            width = "100%"
                        )
                    ),
                    tabPanel(
                        "Labels",
                        div(
                            htmlOutput('labels'),
                            style = "font-size:80%"
                        )
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
                        plotOutput(
                            'matrix', 
                            height = plot.height, 
                            width = "100%"
                        )
                    ),
                    tabPanel(
                        "tSNE",
                        plotOutput(
                            'tSNEplot', 
                            height = plot.height, 
                            width = "100%"
                        )
                    ),
                    tabPanel(
                        "DE",
                        uiOutput('plot_de_genes')
                    ),
                    tabPanel(
                        "Markers",
                        uiOutput('plot_mark_genes')
                    ),
                    tabPanel(
                        "Outliers",
                        plotOutput(
                            'outliers',
                            height = plot.height.small,
                            width = "100%"
                        )
                    )
                )
                do.call(
                    tabsetPanel,
                    c(
                        myTabs,
                        id = "main_panel",
                        type = "pills"
                    )
                )
            })
            
            # MAIN OBSERVER
            # used for extraction of precalculated variables and
            # updating the reactive variables
            # this observer executes whenever any parameter in the left-side
            # parameter panel is changed
            observe({
                res <- prepare_output(object, input$clusters)
                
                # assign results to reactive variables
                values$consensus <- object@consensus$sc3_consensus[[as.character(input$clusters)]]$consensus
                values$labels <- res$labels
                values$labels1 <- res$labels1
                values$hc <- res$hc
                values$silh <- res$silh
                values$new.labels <- res$new.labels
                values$original.labels <- colnames(dataset)[res$cell.order]

                # reindex the new clusters in ascending order
                d <- dataset
                d <- d[ , res$cell.order]
                colnames(d) <- values$new.labels
                values$dataset <- d
                
                # reorder new cell labels in the same order as cells in the 
                # input matrix
                values$original.labels1 <- values$original.labels[res$cell.order]
                values$new.labels1 <- values$new.labels[res$cell.order]
            })
            # this observer updates the GO drop down menu with current clusters
            # after calculating marker genes
            observe({
                markers <- 
                    object@consensus$sc3_biology[[as.character(input$clusters)]]$markers
                mark.res.plot <- mark_gene_heatmap_param(markers)
                clusts <- unique(mark.res.plot$clusts)
                updateSelectInput(session, "cluster", choices = clusts)
            })
            
            ## REACTIVE PANEL DESCRIPTIONS
            output$explanation <- renderUI({
                res <- ""
                if(length(input$main_panel) > 0) {
                    if(grepl("Consensus", input$main_panel)) {
                        res <- HTML("The consensus matrix is a <i>N</i>x<i>N</i> 
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
                                    For a more objective clustering validation please visit the <b>Silhouette</b> panel.")
                    }
                    if(grepl("Silhouette", input$main_panel)) {
                        res <- HTML("As opposed to visual exploration of the 
                                    consensus matrix, a silhouette is a quantitative 
                                    measure of the diagonality of the consensus 
                                    matrix. An average silhouette width 
                                    (shown at the bottom left of the silhouette 
                                    plot) varies from 0 to 1, where 1 represents 
                                    a perfectly block-diagonal consensus matrix 
                                    and 0 represents a situation where there is 
                                    no block-diagonal structure. The best 
                                    clustering is achieved when the average 
                                    silhouette width is close to 1.")
                    }
                    if(grepl("Labels", input$main_panel)) {
                        res <- HTML("The labels panel shows how cells are distributed 
                                    in the clusters that were obtained using the 
                                    original cell indexes from the expression matrix. 
                                    To help visualize how the clustering changes 
                                    when changing from <i>k</i> - 1 to <i>k</i>, the labels are 
                                    colour-coded corresponding to the custering 
                                    results for <i>k</i> - 1.")
                    }
                    if(grepl("Stability", input$main_panel)) {
                        res <- HTML("Stability index shows how stable each cluster
                                    is accross the selected range of <i>k</i>s.
                                    The stability index varies between 0 and 1, where
                                    1 means that the same cluster appears in every
                                    solution for different <i>k</i>.")
                    }
                    if(grepl("Expression", input$main_panel)) {
                        res <- HTML("The expression panel represents the original 
                                    input expression matrix (cells in columns and 
                                    genes in rows) after cell and gene filters. 
                                    Genes are clustered by <i>k</i>-means with <i>k</i> = 100 
                                    (dendrogram on the left) and the heatmap 
                                    represents the expression levels of the gene 
                                    cluster centers after <i>log2</i>-scaling.")
                    }
                    if(grepl("tSNE", input$main_panel)) {
                        res <- HTML("<a href = 'https://lvdmaaten.github.io/tsne/' target='_blank'>tSNE</a> (t-Distributed
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
                                    in a reasonable interval using the slider below.<br><br>")
                    }
                    if(grepl("DE", input$main_panel)) {
                        res <- HTML(paste0("SC3 found <b>", length(values$de.res), "</b> differentially expressed genes based
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
                                    can be changed using the radio buttons below.<br><br>"))
                    }
                    if(grepl("Markers", input$main_panel)) {
                        res <- HTML(paste0("SC3 found <b>", nrow(values$mark.res), "</b> marker genes based
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
                                    slider and radio buttons below.<br>
                                    In addition, SC3 allows you to run a Gene 
                                    Ontology (GO) analysis using a panel below.<br><br>"))
                    }
                    if(grepl("Outliers", input$main_panel)) {
                        res <- HTML("Outlier cells in each cluster are detected 
                                    using robust distances, calculated using 
                                    the minimum covariance determinant (MCD). 
                                    The outlier score shows how different a 
                                    cell is from all other cells in the cluster 
                                    and it is defined as the differences between 
                                    the square root of the robust distance and 
                                    the square root of the 99.99% quantile of 
                                    the Chi-squared distribution (this parameter
                                    can be controlled using a panel below).<br><br>")
                    }
                }
                return(res)
            })
            
            ## REACTIVE PANELS
            
            # plot consensus matrix
            output$consensus <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    if(show.original.labels) {
                        ann <- data.frame(
                            Input.Labels = factor(colnames(dataset))
                        )
                        pheatmap::pheatmap(
                            values$consensus,
                            cluster_rows = values$hc,
                            cluster_cols = values$hc,
                            cutree_rows = input$clusters,
                            cutree_cols = input$clusters,
                            annotation_col = ann,
                            show_rownames = FALSE,
                            show_colnames = FALSE
                        )
                    } else {
                        sc3_plot_consensus(
                            object, 
                            as.numeric(input$clusters)
                        )
                    }
                })
            })
            
            # plot silhouette
            output$silh <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    plot(values$silh, col = "black")
                })
            })
            
            # output cell labels
            output$labels <- renderUI({
                labs <- "<br/><br/>"
                if(!is.null(values$labels1)) {
                    labs1 <- list()
                    cols <- iwanthue(input$clusters - 1)
                    for(i in 1:(input$clusters - 1)) {
                        col <- cols[i]
                        ind <- unlist(strsplit(as.character(values$labels1[i, ]),
                                               " "))
                        for(j in ind) {
                            labs1[[j]] <-
                                paste0("<font color=\"", col, "\">", j, "</font>")
                        }
                    }
                    for(i in 1:input$clusters) {
                        ind <- unlist(
                            strsplit(
                                as.character(values$labels[i, ]),
                                " "
                            )
                        )
                        for(j in ind) {
                            labs <- c(labs, labs1[[j]])
                        }
                        labs <- c(labs, c("<br/>", "<hr>"))
                    }
                } else {
                    for(i in 1:input$clusters) {
                        ind <- unlist(
                            strsplit(
                                as.character(values$labels[i, ]),
                                " "
                            )
                        )
                        labs <- c(labs, ind, c("<br/>", "<hr>"))
                    }
                }
                
                HTML(paste0(labs))
            })
            
            # plot stability
            output$StabilityPlot <- renderPlot({
                validate(
                    need(
                        length(object@consensus$sc3_consensus) > 1,
                        "\nStability cannot be calculated for a single k value!"
                    )
                )
                sc3_plot_cluster_stability(object, input$clusters)
            })
            
            # plot expression matrix
            output$matrix <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    set.seed(1234567)
                    if(show.original.labels) {
                        ann <- data.frame(
                            Input.Labels = factor(values$original.labels)
                        )
                        t <- values$dataset
                        colnames(t) <- rownames(ann)
                        pheatmap::pheatmap(
                            t,
                            kmeans_k = 100,
                            cluster_cols = FALSE,
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            annotation_col = ann,
                            gaps_col = values$col.gaps
                        )
                    } else {
                        sc3_plot_expression(object, as.numeric(input$clusters))
                    }
                })
            })
            
            # make tSNE plot
            output$tSNEplot <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    sc3_plot_tsne(
                        object,
                        input$clusters,
                        input$perplexity
                    )
                })
            })
            
            # plot marker genes
            output$mark_genes <- renderPlot({
                sc3_plot_markers(
                    object, 
                    as.numeric(input$clusters), 
                    as.numeric(input$auroc.threshold), 
                    as.numeric(input$p.val.mark)
                )
            })
            plotHeightMark <- function() {
                validate(
                    need(
                        !is.null(object@consensus$sc3_biology),
                        "\nPlease run sc3_calc_biology() first!"
                    )
                )
                validate(
                    need(
                        try(!is.null(rownames(dataset))),
                        "\nNo gene names provided in the input expression matrix!"
                    )
                )
                markers <- 
                    object@consensus$sc3_biology[[as.character(input$clusters)]]$markers
                mark.res.plot <- mark_gene_heatmap_param(markers)
                return(50 + 10.8*nrow(mark.res.plot))
            }
            
            output$plot_mark_genes <- renderUI({
                plotOutput(
                    "mark_genes", 
                    height = plotHeightMark(), 
                    width = "100%"
                )
            })
            
            # plot DE genes
            output$de_genes <- renderPlot({
                sc3_plot_de_genes(
                    object, 
                    as.numeric(input$clusters), 
                    as.numeric(input$p.val.de)
                )
            })
            plotHeightDE <- function() {
                validate(
                    need(
                        !is.null(object@consensus$sc3_biology),
                        "\nPlease run sc3_calc_biology() first!"
                    )
                )
                validate(
                    need(
                        try(!is.null(rownames(dataset))),
                        "\nNo gene names provided in the input expression matrix!"
                    )
                )
                de.genes <- 
                    object@consensus$sc3_biology[[as.character(input$clusters)]]$de.genes
                return(50 + 10.8*nrow(head(de.genes, 50)))
            }
            output$plot_de_genes <- renderUI({
                plotOutput(
                    "de_genes", 
                    height = plotHeightDE(), 
                    width = "100%"
                )
            })
            
            # plot outliers
            output$outliers <- renderPlot({
                validate(
                    need(
                        !is.null(object@consensus$sc3_biology),
                        "\nPlease run sc3_calc_biology() first!"
                    )
                )
                sc3_plot_cell_outliers(
                    object,
                    as.numeric(input$clusters)
                )
            })

            # REACTIVE BUTTONS
            
            # GO button
            observeEvent(input$GO, {
                open_webgestalt_go(
                    rownames(
                        values$mark.res[values$mark.res[,2] == input$cluster, ]
                    )
                )
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
            
            # stop App on closing the browser
            session$onSessionEnded(function() {
                stopApp()
            })
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
          function(
              object
          ) {
              sc3_interactive.SCESet(
                  object
              )
          })


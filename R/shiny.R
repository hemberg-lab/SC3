#' SC3 interactive function
#'
#' Runs interactive session of SC3 based on precomputed objects
#'
#' @param input.param parameters precomputed by sc3() with interactivity = FALSE (sc3.interactive.arg).
#'
#' @return Opens a browser window with an interactive shiny app and visualize
#' all precomputed clusterings.
#'
#' @importFrom shiny HTML actionButton animationOptions checkboxGroupInput column div downloadHandler downloadLink eventReactive fluidPage fluidRow h4 headerPanel htmlOutput need observe observeEvent p plotOutput reactiveValues renderPlot renderUI selectInput shinyApp sliderInput stopApp tabPanel tabsetPanel uiOutput updateSelectInput validate wellPanel withProgress conditionalPanel reactive outputOptions tags radioButtons downloadButton
#' @importFrom ggplot2 ggplot aes geom_bar geom_point scale_fill_manual scale_color_manual guides theme_bw labs aes_string xlab ylab geom_rug ylim
#' @importFrom utils head write.table
#' @importFrom stats as.dendrogram order.dendrogram cutree median
#' @importFrom graphics plot
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom WriteXLS WriteXLS
#' @importFrom Rtsne Rtsne
#'
#' @export
#'
sc3_interactive <- function(input.param) {
    ## define UI parameters
    plot.height <- 600
    plot.height.small <- 300
    
    ## define server global variables
    values <- reactiveValues()
    
    if(!is.null(input.param$svm.num.cells)) {
        with_svm <- TRUE
        values$svm <- FALSE
    } else {
        with_svm <- FALSE
        values$svm <- TRUE
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
                        paste0("<h3>SC3 clustering of ", 
                               input.param$filename, 
                               "</h3>")
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
                            sliderInput(
                                "clusters",
                                label = "Number of clusters k",
                                min = min(as.numeric(unlist(input.param$cons.table[,3]))),
                                max = max(as.numeric(unlist(input.param$cons.table[,3]))),
                                value = median(as.numeric(unlist(input.param$cons.table[,3]))),
                                step = 1,
                                animate = animationOptions(
                                    interval = 2000,
                                    loop = FALSE
                                )
                            )
                        )
                    ),
                    fluidRow(
                        column(
                            12,
                            conditionalPanel(
                                "!output.is_svm",
                                h4("SVM prediction"),
                                p(
                                    paste0(
                                        "Your data was clustered based on ",
                                        input.param$svm.num.cells,
                                        " cells. When you have found the best clustering parameters, press this button to predict labels of the other cells and to perform biological interpretation:\n\n"
                                    )
                                ),
                                actionButton(
                                    "svm",
                                    label = "Run SVM"
                                )
                            )
                        )
                    )),
                    h4("Export results"),
                    wellPanel(
                        conditionalPanel(
                            "output.is_svm",
                            fluidRow(
                                column(6,
                                       downloadButton(
                                           'excel',
                                           label = "To Excel")
                                ),
                                column(6,
                                       actionButton(
                                           'save',
                                           label = "To R session")
                                )
                            ),
                            fluidRow(
                                column(12,
                                       HTML("<br>"),
                                       downloadButton(
                                           'consens_download',
                                           label = "Consensus")
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
                    ),
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
                                    conditionalPanel(
                                        "output.is_mark",
                                        selectInput(
                                            "cluster",
                                            "Choose a cluster for GO analysis:",
                                            c("None" = "NULL")
                                        ),
                                        if(input.param$rselenium.installed) {
                                            actionButton("GO", label = "Analyze!")
                                        },
                                        if(input.param$rselenium.installed) {
                                            HTML("<br>(opens <a href = 'http://bioinfo.vanderbilt.edu/webgestalt/' target='_blank'>WebGestalt</a> in Firefox)")
                                        } else {
                                            HTML("<font color='red'>To be able to run GO analysis you need to install RSelenium library. You can do that by closing this window and then running 'RSelenium::checkForServer()' command in your R session. This will download the required library. After that please rerun SC3 again. More details are available <a href = 'https://cran.r-project.org/web/packages/RSelenium/vignettes/RSelenium-basics.html' target='_blank'>here</a></font>.")
                                        }
                                    )
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
                                        min = floor(0.7 * ncol(input.param$dataset) / 5),
                                        max = floor(1.3 * ncol(input.param$dataset) / 5),
                                        value = floor(ncol(input.param$dataset) / 5)
                                    )
                                )
                            ),
                            conditionalPanel(
                                "input.main_panel == 'Outliers'",
                                wellPanel(
                                    radioButtons(
                                        "chisq.quantile",
                                        label = "% of the Chi-squared quantile",
                                        choices = c(
                                            "95%" = 0.95,
                                            "99%" = 0.99,
                                            "99.99%" = 0.9999
                                        ),
                                        selected = "0.9999",
                                        inline = TRUE
                                    )
                                )
                            )
                        )
                    )
                )
            )
        ),
        
        server = function(input, output, session) {
            # render tabpanel
            output$mytabs = renderUI({
                if(values$svm) {
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
                } else {
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
                        )
                    )
                }
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
                # get all results for k
                res <- input.param$cons.table[
                    as.numeric(input.param$cons.table[ , 3]) == input$clusters, 
                    4][[1]]
                # get all results for k-1
                values$labels1 <- NULL
                if((input$clusters - 1) %in% as.numeric(input.param$cons.table[ , 3])) {
                    res1 <- input.param$cons.table[
                        as.numeric(input.param$cons.table[ , 3]) == (input$clusters - 1), 
                        4][[1]]
                    values$labels1 <- res1[[2]]
                }
                # assign results to reactive variables
                values$consensus <- res[[1]]
                values$labels <- res[[2]]
                values$hc <- res[[3]]
                values$silh <- res[[4]]
                # order cells according to hierarchical clustering of the
                # consensus matrix
                clusts <- cutree(values$hc, input$clusters)
                cell.order <- order.dendrogram(as.dendrogram(values$hc))
                d <- input.param$dataset
                colnames(d) <- clusts
                d <- d[ , cell.order]
                values$original.labels <- colnames(input.param$dataset)[cell.order]
                # prepare a consensus matrix for downloading
                tmp <- data.frame(values$consensus[cell.order, cell.order])
                colnames(tmp) <- values$original.labels
                rownames(tmp) <- seq(length=nrow(tmp))
                values$consensus.download <- tmp
                # reindex the new clusters in ascending order
                values$new.labels <- reindex_clusters(colnames(d))
                colnames(d) <- values$new.labels
                values$dataset <- d
                
                # calculate stability of the clusters
                
                # check if there are more than 1 k value in ks range
                values$stability <- NULL
                if(length(unique(as.numeric(input.param$cons.table[ , 3]))) > 1) {
                    stab.res <- input.param$cons.table[ , 3:4]
                    values$stability <- StabilityIndex(
                        stab.res,
                        input$clusters
                    )
                }
                
                # reorder new cell labels in the same order as cells in the 
                # input matrix
                ord <- order(cell.order)
                values$original.labels1 <- values$original.labels[ord]
                values$new.labels1 <- values$new.labels[ord]
                # define gaps between clusters for the heatmap
                values$col.gaps <- which(diff(as.numeric(colnames(d))) != 0)
                # update reactive variables used for hiding some panels on the
                # webpage
                values$mark <- FALSE
                # update the svm reactive variable
                if(with_svm) {
                    values$svm <- FALSE
                }
                # update output variables
                values$cells.outliers <- data.frame()
                values$de.genes <- data.frame()
                values$marker.genes <- data.frame()
                SC3.results <- list()
                # update the output cell labels
                values$cell.labels <- 
                    data.frame(new.labels = as.numeric(values$new.labels),
                               original.labels = values$original.labels,
                               stringsAsFactors = FALSE)
                values$cell.labels1 <- 
                    data.frame(original.labels = values$original.labels1,
                               new.labels = as.numeric(values$new.labels1),
                               stringsAsFactors = FALSE)
            })
            # this observer updates the GO drop down menu with current clusters
            # after calculating marker genes
            observe({
                if(values$mark) {
                    clusts <- unique(values$mark.res$clusts)
                    updateSelectInput(session, "cluster", choices = clusts)
                }
                
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
                                    combinations of clustering parameters 
                                    (checkboxes in the left panel). Similarity 0  
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
            output$consensus <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    if(input.param$show.original.labels) {
                        ann <- data.frame(
                            Input.Labels = factor(colnames(input.param$dataset))
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
                        pheatmap::pheatmap(
                            values$consensus,
                            cluster_rows = values$hc,
                            cluster_cols = values$hc,
                            cutree_rows = input$clusters,
                            cutree_cols = input$clusters,
                            show_rownames = FALSE,
                            show_colnames = FALSE
                        )
                    }
                })
            })
            
            output$silh <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    plot(values$silh, col = "black")
                })
            })
            
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
            
            output$StabilityPlot <- renderPlot({
                validate(
                    need(
                        !is.null(values$stability),
                        "\nStability cannot be calculated for a single k value!"
                    )
                )
                withProgress(message = 'Plotting...', value = 0, {
                    d <- data.frame(
                        Cluster = factor(1:length(values$stability)),
                        Stability = values$stability)
                    ggplot(d, aes(x = d$Cluster, y = d$Stability)) +
                        geom_bar(stat = "identity") +
                        ylim(0, 1) +
                        labs(x = "Cluster", y = "Stability Index") +
                        theme_bw()
                })
            })
            
            output$matrix <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    set.seed(1234567)
                    if(input.param$show.original.labels) {
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
                        pheatmap::pheatmap(
                            values$dataset,
                            kmeans_k = 100,
                            cluster_cols = FALSE,
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            gaps_col = values$col.gaps
                        )
                    }
                })
            })
            
            output$tSNEplot <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    set.seed(1234567)
                    tsne_out <- Rtsne::Rtsne(
                        t(values$dataset),
                        perplexity = input$perplexity
                    )
                    df_to_plot <- as.data.frame(tsne_out$Y)
                    df_to_plot$Cluster <- factor(
                        colnames(values$dataset),
                        levels = unique(colnames(values$dataset))
                    )
                    comps <- colnames(df_to_plot)[1:2]
                    ggplot(df_to_plot, aes_string(x = comps[1],
                                                  y = comps[2],
                                                  color = "Cluster")) +
                        geom_point() +
                        xlab("Dimension 1") +
                        ylab("Dimension 2") +
                        theme_bw()
                })
            })
            output$mark_genes <- renderPlot({
                d <- values$plot.mark[[1]]
                d.param <- values$plot.mark[[2]]
                col.gaps <- values$plot.mark[[3]]
                pheatmap::pheatmap(
                    d[rownames(d.param$mark.res.plot), , drop = FALSE],
                    show_colnames = FALSE,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    annotation_row = d.param$row.ann,
                    annotation_names_row = FALSE,
                    gaps_row = d.param$row.gaps,
                    gaps_col = col.gaps,
                    cellheight = 10
                )
            })
            plotHeightMark <- function() {
                validate(
                    need(
                        try(!is.null(rownames(input.param$dataset))),
                        "\nNo gene names provided in the input expression matrix!"
                    )
                )
                withProgress(message = 'Calculating Marker genes...',
                             value = 0,
                             {
                                 # prepare dataset for plotting
                                 if(with_svm) {
                                     d <- values$dataset.svm
                                     col.gaps <- values$col.gaps.svm
                                 } else {
                                     d <- values$dataset
                                     col.gaps <- values$col.gaps
                                 }
                                 
                                 values$mark <- FALSE
                                 
                                 # define marker genes
                                 values$mark.res <- get_marker_genes(
                                     d,
                                     as.numeric(colnames(d)),
                                     as.numeric(input$auroc.threshold),
                                     as.numeric(input$p.val.mark)
                                 )
                                 # check the results of mark_genes_main:
                                 # global variable mark.res
                                 validate(
                                     need(
                                         try(dim(values$mark.res)[1] != 0),
                                         "\nUnable to find significant marker genes from obtained clusters! Try to change either the number of clusters k, the AUROC threshold or the p-value threshold."
                                     )
                                 )
                                 colnames(values$mark.res) <- 
                                     c("AUC","clusts","p.value")
                                 d.param <- mark_gene_heatmap_param(
                                     values$mark.res,
                                     unique(colnames(d))
                                 )
                                 values$mark <- TRUE
                                 values$marker.genes <- data.frame(
                                     new.labels = as.numeric(values$mark.res[,2]),
                                     gene = rownames(values$mark.res),
                                     AUROC = as.numeric(values$mark.res[,1]),
                                     p.value = as.numeric(values$mark.res[,3]),
                                     stringsAsFactors = FALSE
                                 )
                            }
                        )
                values$plot.mark <- list(d, d.param, col.gaps)
                return(50 + 10.8*nrow(d.param$mark.res.plot))
            }
            
            output$plot_mark_genes <- renderUI({
                plotOutput(
                    "mark_genes", 
                    height = plotHeightMark(), 
                    width = "100%"
                )
            })
            output$de_genes <- renderPlot({
                d <- values$plot.de[[1]]
                d.param <- values$plot.de[[2]]
                col.gaps <- values$plot.de[[3]]
                pheatmap::pheatmap(
                    d[names(head(values$de.res, 50)), , drop = FALSE],
                    show_colnames = FALSE,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    annotation_row = d.param$row.ann,
                    annotation_names_row = FALSE,
                    gaps_col = col.gaps,
                    cellheight = 10
                )
            })
            plotHeightDE <- function() {
                validate(
                    need(try(!is.null(rownames(input.param$dataset))),
                         "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating DE genes...', value = 0, {
                    # prepare dataset for plotting
                    if(with_svm) {
                        d <- values$dataset.svm
                        col.gaps <- values$col.gaps.svm
                    } else {
                        d <- values$dataset
                        col.gaps <- values$col.gaps
                    }
                    # define de genes
                    values$de.res <- get_de_genes(
                        d, 
                        colnames(d), 
                        as.numeric(input$p.val.de)
                    )
                    # check the results of de_genes_main:
                    # global variable de.res
                    validate(
                        need(
                            try(length(values$de.res) != 0),
                            "\nUnable to find significant differentially expressed genes from obtained clusters! Try to change either the number of clusters k or the p-value threshold."
                        )
                    )
                    d.param <- de_gene_heatmap_param(head(values$de.res, 50))
                    res <- data.frame(
                        gene = names(values$de.res),
                        p.value = as.numeric(values$de.res),
                        stringsAsFactors = FALSE
                    )
                    rownames(res) <- NULL
                    values$de.genes <- res
                })
                values$plot.de <- list(d, d.param, col.gaps)
                return(50 + 10.8*length(head(values$de.res, 50)))
            }
            output$plot_de_genes <- renderUI({
                plotOutput("de_genes", height = plotHeightDE(), width = "100%")
            })
            
            
            output$outliers <- renderPlot({
                withProgress(
                    message = 'Calculating cell outliers...',
                    value = 0, 
                    {
                        # prepare dataset for plotting
                        if(with_svm) {
                            d <- values$dataset.svm
                        } else {
                            d <- values$dataset
                        }
                        
                        # compute outlier cells
                        values$outl.res <- get_outl_cells(
                            d,
                            colnames(d),
                            as.numeric(input$chisq.quantile)
                        )
                        
                        t <- as.data.frame(values$outl.res)
                        colnames(t)[1] <- "outl"
                        t$Cluster <- names(values$outl.res)
                        t$Cells <- 1:dim(t)[1]
                        t$Cluster <-
                        factor(
                            t$Cluster,
                            levels =
                                unique(
                                    as.character(
                                        sort(as.numeric(t$Cluster))
                                    )
                                )
                        )
                        cols <- iwanthue(length(unique(t$Cluster)))
                        Cells <- outl <- Cluster <- NULL
                        
                        values$cells.outliers <- if(with_svm) {
                            data.frame(
                                new.labels = as.numeric(names(values$outl.res)),
                                original.labels = values$original.labels.svm,
                                MCD.dist = as.numeric(values$outl.res),
                                stringsAsFactors = FALSE
                            )
                        } else {
                            data.frame(
                                new.labels = as.numeric(names(values$outl.res)),
                                original.labels = values$original.labels,
                                MCD.dist = values$outl.res,
                                stringsAsFactors = FALSE
                            )
                        }
                        
                        ggplot(t, aes(x = t$Cells,
                                  y = t$outl,
                                  fill = t$Cluster, 
                                  color = t$Cluster)) +
                        geom_bar(stat = "identity") +
                        geom_point() +
                        scale_fill_manual(values = cols) +
                        scale_color_manual(values = cols) +
                        guides(color = FALSE, fill = FALSE) +
                        labs(x = "Cells", y = "Outlier score") +
                        theme_bw()
                    }
                )
            })

            # REACTIVE BUTTONS
            
            # SVM button
            observeEvent(input$svm, {
                withProgress(message = 'Running SVM...', value = 0, {
                    d <- input.param$dataset
                    colnames(d) <- values$new.labels1
                    original.labels <-
                        c(
                            values$original.labels1,
                            colnames(input.param$study.dataset)
                         )
                    
                    values$dataset.svm <- input.param$study.dataset
                    colnames(values$dataset.svm) <-
                        support_vector_machines(
                            d,
                            input.param$study.dataset,
                            "linear"
                        )
                    
                    tmp <- cbind(d, values$dataset.svm)
                    cols <- colnames(tmp)
                    inds <- NULL
                    for(i in unique(colnames(tmp))) {
                        inds <- c(inds, which(cols == i))
                    }
                    tmp <- tmp[ , inds]
                    
                    values$original.labels.svm <- original.labels[inds]
                    colnames(tmp) <- reindex_clusters(colnames(tmp))
                    values$new.labels.svm <- colnames(tmp)
                    values$dataset.svm <- tmp
                    
                    values$col.gaps.svm <-
                        which(
                            diff(as.numeric(colnames(values$dataset.svm))) != 0
                        )
                    
                    values$cell.labels <- data.frame(
                        new.labels = values$new.labels.svm,
                        original.labels = values$original.labels.svm
                    )
                    # original ordering of the cells
                    ord <- order(input.param$svm.inds)
                    values$cell.labels1 <- data.frame(
                        original.labels = original.labels[ord],
                        new.labels = values$new.labels.svm[order(inds)][ord]
                    )
                    
                    values$svm <- TRUE
                    values$svm.clusters <- paste(input$clusters, collapse = "_")
                    values$svm.distance <- paste(input$distance, collapse = "_")
                    values$svm.dimRed <- paste(input$dimRed, collapse = "_")
                })
            })
            
            # SAVE TO R session button
            observeEvent(input$save, {
                SC3.results <<- list(
                    cell.labels = values$cell.labels,
                    cell.labels.original.order = values$cell.labels1,
                    consensus = values$consensus.download,
                    de.genes = values$de.genes,
                    marker.genes = values$marker.genes,
                    cells.outliers = values$cells.outliers
                )
            })
            
            # GO button
            observeEvent(input$GO, {
                open_webgestalt_go(
                    rownames(
                        values$mark.res[values$mark.res[,2] == input$cluster, ]
                    )
                )
            })
            
            ## OUTPUTS USED FOR CONDITIONAL PANELS
            output$is_svm <- reactive({
                return(values$svm)
            })
            output$is_mark <- reactive({
                return(values$mark)
            })

            # DOWNLOAD buttons
            
            # Save to xls button
            output$excel <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-", input.param$filename, ".xls")
                },
                
                content = function(file) {
                    WriteXLS(list(values$cell.labels,
                                  values$cell.labels1,
                                  values$de.genes,
                                  values$marker.genes,
                                  values$cells.outliers),
                             ExcelFileName = file,
                             SheetNames = c("Cell labels",
                                            "Cell labels (original order)",
                                            "DE genes",
                                            "Marker genes",
                                            "Cell outliers"))
                }
            )
            
            # Save consensus matrix button
            output$consens_download <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-consensus-", input.param$filename, ".txt")
                },
                
                content = function(file) {
                    write.table(values$consensus.download,
                                file = file,
                                row.names = FALSE,
                                quote = FALSE,
                                sep = "\t")
                }
            )
            
            # stop App on closing the browser
            session$onSessionEnded(function() {
                stopApp()
            })
            
            # hide panels
            outputOptions(output, 'is_mark', suspendWhenHidden = FALSE)
            outputOptions(output, 'is_svm', suspendWhenHidden = FALSE)
        },
        # launch App in a browser
        options = list(launch.browser = TRUE)
    )
}

library(gtools)
library(shiny)
library(shinyFiles)
library(rhdf5)
library(SparseM)
library(ggplot2)
library(reshape2)

shinyServer(function(input, output, session) {
  
  # --------------------
  # Setup/Load Data Page
  # --------------------
  # Use shinyFiles library for selecting local files
  # The built-in fileInput copies the file -- not necessary
  shinyFileChoose(input, 'timeseries_database', root=list(home='~',root='/'), filetypes=c('h5'), session=session)
  
  output$tst_filename <- renderPrint({
    if (is.null(input$timeseries_database))
      return("Please select your time-series database file")
    as.character(parseFilePaths(c(home='~'), input$timeseries_database)$datapath[1])
  })

  # Declare some variables here so we can viciously abuse R's scoping
  # In the future, I should figure out if Shiny can handle this properly
  # because this has gotten a bit out of hand
  plot_width <- 12
  plot_height <- 4
  max_cluster <- NULL
  min_param <- NULL
  step_size <- NULL
  taxonomic_ids <- NULL
  sequence_ids <- NULL
  sequencecluster_labels <- NULL
  tp <- NULL
  modCRTLIST <- NULL
  tsdatabase_path <- NULL
  nsequences <- NULL
  nsamples <- NULL
  nparams <- NULL
  nobs <- NULL
  timeseries_table <- NULL
  timeseries_totals <- NULL
  cluster_names <- NULL
  file_version <- NULL
  mask <- NULL
  values <- reactiveValues()
  
  loadFilesResult <- observeEvent(input$loadFiles, {
    withProgress(message = 'Loading database file, please wait...', value = 0, {
      tsdatabase_path <<- normalizePath(as.character(parseFilePaths(c(home='~'), input$timeseries_database)$datapath[1]))
      #First, load the cluster file
      
      h5data <- h5ls(tsdatabase_path)
      nobs <<- as.numeric(h5data[h5data$name=="sequenceids",]$dim)
      nsamples <<- as.numeric(h5data[h5data$name=="names",]$dim)
      nparams <<- as.numeric(strsplit(h5data[h5data$name=="clusters",]$dim,' x')[[1]][1])
      incProgress(0.1, detail = "Reticulating splines...")
      #Time cluster labels loaded on the fly
      #Time-series table
      tslength<-as.numeric(h5data[h5data$name=="data",]$dim)
      #Small enough to read in one go
      tsindptr<-h5read(tsdatabase_path,"timeseries/indptr")
      incProgress(0.25, detail = "De-chunking time dimension...")
      #The time-series data must be read in chunks to prevent huge memory
      #usage by the HDF5 library
      tsdata<-vector(length=tslength)
      tsindices<-vector(length=tslength)
      chunks<-seq(1,tslength,10000)
      #Make sure we get all the chunks
      if (chunks[length(chunks)] != tslength) {
        chunks<-c(chunks,tslength)
      }
      for (i in 1:(length(chunks)-1)){
        tsdata[chunks[i]:chunks[i+1]] <- h5read(tsdatabase_path,
             "timeseries/data",index=list(chunks[i]:chunks[i+1]))
        tsindices[chunks[i]:chunks[i+1]] <- h5read(tsdatabase_path,
             "timeseries/indices",index=list(chunks[i]:chunks[i+1]))
      }
      #Correct for R's 1-based indexing
      tsindices <- tsindices+1
      tsindptr <- tsindptr+1
      #Make a sparse matrix, immediately converting it to a full matrix
      #This is necessary because the sparse classes can't do everything
      #that the full class can do
      timeseries_table <<- as.matrix(new("matrix.csr",ra=as.numeric(tsdata),
                 ja=as.integer(tsindices),ia=as.integer(tsindptr),
                 dimension=as.integer(c(nobs,nsamples))))
      timeseries_totals <<- colSums(timeseries_table)
      nsequences <<- sum(timeseries_totals)
      incProgress(0.5, detail="Identifying life-forms...")
      #Read in auxiliary information
      taxonomic_ids <<- h5read(tsdatabase_path, "genes/taxonomy")
      sequencecluster_labels <<- h5read(tsdatabase_path, "genes/sequenceclusters")
      sequence_ids <<- h5read(tsdatabase_path, "genes/sequenceids")
      #Grab some of the required numbers
      incProgress(0.75, detail="Uploading mission parameters...")
      tp <<- h5read(tsdatabase_path,"samples/time")
      cluster_params <- h5readAttributes(tsdatabase_path,"genes/clusters")
      min_param <<- cluster_params$param_min
      max_param <<- cluster_params$param_max
      step_size <<- cluster_params$param_step
      file_version <- h5readAttributes(tsdatabase_path,"/")$origin_version
      if (! is.null(file_version)) {
          file_version <- strsplit(file_version, ".", fixed=TRUE)[[1]]
          if ((file_version[2] >= 1) & (file_version[2] >= 1)) {
              #This dataset exists in version >= 0.1.1
              mask <<- h5read(tsdatabase_path,"samples/mask")
          } else {
              mask <<- rep(1, length(tp))
          }
      } else {
          mask <<- rep(1, length(tp))
      }
      incProgress(1, detail="Triangulation complete.")
    })
  })
  
  output$tsdb_filename <- renderText({
    if (is.null(input$timeseries_database)) {
        "Please select a file"
    } else {
      normalizePath(as.character(parseFilePaths(c(home='~'), input$timeseries_database)$datapath[1]))
    }
  })
  
  # --------------------
  # Data Summary Page
  # --------------------
  
  #Plot the number of clusters vs epsilon
  nclusters <- eventReactive(input$plotnclust, {
    withProgress(message = 'Retrieving clusters from database file...', value = 0, {
      nclusts <- c()
      noise_size <- c()
      for (i in 1:nparams) {
        incProgress(1/nparams)
        clusts <- h5read(tsdatabase_path,"genes/clusters",index=list(i,NULL))
        noise_size <- c(noise_size, sum(timeseries_table[which(clusts == "-1"),]))
        nclusts <- c(nclusts, length(unique(t(clusts))))
      }
      list(nclusts=nclusts,noise_size=noise_size)
    })
  })
  
  output$clusterVsEps <- renderPlot({
      layout(as.matrix(cbind(1,2,3,4)))
      params <- 0:(nparams-1)*step_size+min_param
      res <- nclusters()
      nclusts <- res$nclusts
      noise_size <- res$noise_size
      plot(params, nclusts, type='l', xlab="epsilon", ylab="Number of time-series clusters", cex.lab=1.5, cex.axis=1.5)
      plot(params, noise_size/nsequences, type='l', xlab="epsilon", ylab="Total sequences in noise bin", cex.lab=1.5, cex.axis=1.5)
      output$maxCluster<<-renderText({
        as.numeric(max(nclusts))
      })
      max_eps_ind <- which(nclusts == max(nclusts))[1]
      output$maxEps<<-renderText({
        as.numeric(params[max_eps_ind])
      })
      output$epsNoise<<-renderText({
        as.numeric(noise_size[max_eps_ind]/nsequences)
      })
      clusts <- as.vector(h5read(tsdatabase_path,"genes/clusters",index=list(max_eps_ind,NULL)))
      clust_totals <- aggregate(rowSums(timeseries_table), by=list(clusts), FUN=sum)
      hist(log10(clust_totals$x), breaks=100, xlab="Log10 cluster size (number of sequences)", cex.lab=1.5, cex.axis=1.5)
      unique_seqs <- table(clusts)
      hist(log10(unique_seqs), breaks=100, xlab="Log10 cluster size (number of unique sequences)", cex.lab=1.5, cex.axis=1.5)
      totSeqsAtMax <- mean(clust_totals$x)
      uniqSeqsAtMax <- mean(unique_seqs)
      output$clustSize <- renderText({
          paste("Mean Total Sequences: ", totSeqsAtMax, ", Mean Unique Sequences: ", uniqSeqsAtMax, sep="")
      })
  }, height=300, width=1200)
  
  #Simpson Index (list items must be proportions)
  Simpson <- function(x) Reduce("+", x^2)
  
  tax_consistency <- eventReactive(input$plottemptaxconsist, {
      tax_levels <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
      #Requires GreenGenes style classifications
      short_code <- c("k__","p__","c__","o__","f__","g__","s__")
      taxdf <- data.frame()
      for (level in 2:7) {
      simpson_indices <- c()
      for (clusterid in unique(values$current_clusters)) {
      #Exclude noise from calculation
        if (clusterid != "-1") {
          tax_ids <- taxonomic_ids[values$current_clusters==clusterid]
          split_tax<- sapply(tax_ids,function(x) strsplit(x,";")[[1]][level])
          taxtable <- table(as.factor(split_tax))
          taxtable <- taxtable[rownames(taxtable!=short_code[level])]
          proptaxtable <- taxtable/sum(taxtable)
          simpson_index <- Simpson(proptaxtable)
          #Simpson index may be NULL if table is empty
          #(i.e., no classifications at level)
          if (!is.null(simpson_index))
            simpson_indices <- c(simpson_indices, simpson_index)
        }
      }
          taxdf <- rbind(taxdf,
                         data.frame(values=simpson_indices,
                         level=rep(tax_levels[level], length(simpson_indices))))
      }
      taxdf
  })
      
  output$temptaxconsist <- renderPlot({
      #Plotting temporal taxonomic consistency and output table
      output$temptaxconsisttab <- renderTable({
        taxdf <- tax_consistency()
        aggregate(taxdf$values, by=list(taxdf$level), FUN=mean)
      })
      p<-ggplot(tax_consistency(), aes(y=values, x=level))                                                     
      p+geom_boxplot()
  })
  

  
  # --------------------
  # Explore Clusters Page
  # --------------------
  
  filterClusters <- function() {
    if (input$selectTaxon != "") {
      selectIndices <- NULL
      try({
          selectIndices <- grep(input$selectTaxon, taxonomic_ids, ignore.case=TRUE)
      })
      if (length(selectIndices) > 0) {
          values$selected_clusters <<- values$current_clusters[selectIndices]
      } else {
        values$selected_clusters <<- NULL
      }
    } else {
      values$selected_clusters <<- values$current_clusters
    }
  }
  
  filterSeqClusters <- function() {
    if (input$selectTaxonSeqClust != "") {
      selectIndices <- NULL
      try({
        selectIndices <- grep(input$selectTaxonSeqClust, taxonomic_ids, ignore.case=TRUE)
      })
      if (length(selectIndices) > 0) {
        values$selected_seqclusters <<- sequencecluster_labels[selectIndices]
      } else {
        values$selected_seqclusters <<- NULL
      }
    } else {
      values$selected_seqclusters <<- sequencecluster_labels
    }
  }
  
  #Populate the values$current_clusters variable with the cluster names for
  #the given clustering parameter
  changeClusterResult <- observeEvent(input$cluster_param, {
    cluster_index <- (input$cluster_param-min_param)/step_size+1
    values$current_clusters <<- as.vector(h5read(tsdatabase_path,"genes/clusters",index=list(cluster_index,1:nobs)))
    values$selected_clusters <<- values$current_clusters
    #If we have a taxa filter string, filter out clusters that don't contain that string
    filterClusters()
  }, ignoreNULL = TRUE)
  
  filterTaxa <- observeEvent(input$selectTaxon, {
      filterClusters()
  })
  
  filterTaxaSeq <- observeEvent(input$selectTaxonSeqClust, {
      filterSeqClusters()
  })
  
  #Slider for selecting epsilon parameter
  output$clusterParamSelector <- renderUI({
    numericInput("cluster_param", "Cluster parameter (eps):", 
                 min = min_param, max = max_param, step = step_size, value = min_param)
  })
  
  output$plotConsistButton <- renderUI({
      if (is.null(input$cluster_param))
        return()
      actionButton("plottemptaxconsist",paste("Plot Temporal/Taxonomic Consistency for Epsilon =",input$cluster_param))
  })
  
  #Select time series cluster widget
  output$clusterSelector <- renderUI({
    if (is.null(input$cluster_param))
      return()
    cluster_abunds <- aggregate(rowSums(timeseries_table), by=list(values$current_clusters), FUN=sum)
    cluster_names <- as.numeric(cluster_abunds[order(cluster_abunds[,2], decreasing=TRUE),1])
    cluster_names <- cluster_names[cluster_names %in% values$selected_clusters]
    if (length(cluster_names) == 0) {
      stop("No clusters matching taxa filter.")
      return()
    }
    if ((length(cluster_names) > 1) & (cluster_names[1] == -1)) {
        select <- cluster_names[2]
    } else {
        select <- cluster_names[1]
    }
    selectInput('cluster_number', 'Cluster number:', choices = cluster_names,
                multiple = FALSE, selectize = FALSE, selected = select)
  })
  
  #Select OTU number widget
  output$OTUSelector <- renderUI({
    filterSeqClusters()
    cluster_abunds <- aggregate(rowSums(timeseries_table), by=list(sequencecluster_labels), FUN=sum)
    cluster_names <- cluster_abunds[order(cluster_abunds[,2], decreasing=TRUE),1]
    cluster_names <- cluster_names[cluster_names %in% values$selected_seqclusters]
    if (length(cluster_names) == 0) {
      stop("No clusters matching taxa filter.")
      return()
    }
    select <- cluster_names[1]
    selectInput('otu_number', 'OTU Cluster number:', choices = cluster_names,
                multiple = FALSE, selectize = FALSE, selected = select)
  })

  #Main time series cluster plot
  output$main_plot <- renderPlot({
    if (is.null(input$cluster_param))
      return()
    if (is.null(input$cluster_number))
      return()
    if (as.numeric(input$cluster_number) == -1) {
       stop("The \"-1\" cluster contains sequences DBSCAN has labelled as noise. This plot is often too large to produce, so it is disabled. See below for a list of sequences that were classified as noise at the given epsilon parameter.")
       return()
    }
    # Take only the time-series that are in the current cluster
    subset_table <- timeseries_table[values$current_clusters==input$cluster_number, ]
    if (dim(subset_table)[1]>0) {
        plotWrapper(subset_table)
    } else {
        stop("The selected subset contains no data.")
        return()
    }
  }, height=plot_height*100, width=plot_width*100)
  
  # Plot by OTU number
  output$otu_plot <- renderPlot({
    input$cluster_param
    if (is.null(input$otu_number))
      return()
    # Take only the time-series that are in the current cluster
    subset_table <- timeseries_table[sequencecluster_labels==input$otu_number, ]
    subset_table <- as.matrix(subset_table)
    # Fix behaviour when there's only one row, it makes it
    # column-wise, so we have to transpose it
    if (dim(subset_table)[2] == 1) {
      subset_table <- t(subset_table)
    }
    plotWrapper(subset_table, otu_plot=TRUE)
  }, height=plot_height*100, width=plot_width*100)
  
  #Generates a table with the abundance, taxonomy, and cluster information
  output$infotable <- renderDataTable({
    if (is.null(input$cluster_param))
      return()
    if (is.null(input$cluster_number))
      return()
    subset_ids <- sequence_ids[values$current_clusters==input$cluster_number]
    subset_table <- timeseries_table[values$current_clusters==input$cluster_number, ]
    tax_ids <- taxonomic_ids[values$current_clusters==input$cluster_number]
    phylo_clusters <- sequencecluster_labels[values$current_clusters==input$cluster_number]
    data.frame(Abundance=rowSums(subset_table),
               PhyloClusterNumber=phylo_clusters,
               TaxonomicID=tax_ids,
               SequenceID=subset_ids)
  })
  
  #Generates a table with the abundance, taxonomy, and cluster information
  output$otutable <- renderDataTable({
    if (is.null(input$otu_number))
      return()
    subset_ids <- sequence_ids[sequencecluster_labels==input$otu_number]
    subset_table <- timeseries_table[sequencecluster_labels==input$otu_number, ]
    subset_table <- as.matrix(subset_table)
    if (dim(subset_table)[2] == 1) {
      subset_table <- t(subset_table)
    }
    tax_ids <- taxonomic_ids[sequencecluster_labels==input$otu_number]
    time_clusts <- values$current_clusters[sequencecluster_labels==input$otu_number]
    df <- data.frame(Abundance=rowSums(subset_table),
               TimeClustNumber=time_clusts,
               TaxonomicID=tax_ids,
               SequenceID=subset_ids)
    if (input$excludeNoise) {
      df <- df[df$TimeClustNumber!=-1,]
    }
    df
  })
  
  output$saveTable <- downloadHandler(
    filename <- function() {
      paste('timeseries_cluster-eps_', input$cluster_param, "-num_", input$cluster_number, '.csv', sep="")
    },
    content <- function(file) {
        subset_ids <- sequence_ids[values$current_clusters==input$cluster_number]
        subset_table <- timeseries_table[values$current_clusters==input$cluster_number, ]
        tax_ids <- taxonomic_ids[values$current_clusters==input$cluster_number]
        phylo_clusters <- sequencecluster_labels[values$current_clusters==input$cluster_number]
        save_df <- data.frame(SequenceID=subset_ids,
               Abundance=rowSums(subset_table),
               TaxonomicID=tax_ids,
               PhyloClusterNumber=phylo_clusters)
        time_series <- as.data.frame(timeseries_table[values$current_clusters==input$cluster_number, ])
        colnames(time_series) <- tp
        save_df <- cbind(save_df, time_series)
        write.csv(save_df, file)
    }
  )
  
  output$saveAll <- downloadHandler(
    filename <- function() {
      paste('all_clusters-eps_', input$cluster_param, '.csv', sep="")
    },
    content <- function(file) {
        ids <- sequence_ids
        full_table <- timeseries_table
        tax_ids <- taxonomic_ids
        phylo_clusters <- sequencecluster_labels
        time_clusts <- values$current_clusters
        save_df <- data.frame(SequenceID=ids,
               Abundance=rowSums(full_table),
               TaxonomicID=tax_ids,
               PhyloClusterNumber=phylo_clusters,
               TimeCluster=time_clusts)
        time_series <- as.data.frame(timeseries_table)
        colnames(time_series) <- tp
        save_df <- cbind(save_df, time_series)
        write.csv(save_df, file)
    }
  )
  
  # Save time cluster plot as SVG
  output$saveMainPlot <- downloadHandler(
  filename <- function() {
      paste('timeseries_cluster-eps_', input$cluster_param, "-num_", input$cluster_number, '.svg', sep="")
    },
    content <- function(file) {
      if (is.null(input$cluster_param))
        return()
      if (is.null(input$cluster_number))
        return()
      # Take only the time-series that are in the current cluster
      subset_table <- timeseries_table[values$current_clusters==input$cluster_number, ]   
      svg("tsplot.svg",height=plot_height,width=plot_width)
      # Plotting code for main time-series plot
      layout(as.matrix(cbind(c(1,3),c(2,4))))
      plotWrapper(subset_table)
      dev.off()
      file.copy("tsplot.svg", file)
    }
  )
  
  output$saveOTUTable <- downloadHandler(
    filename <- function() {
      paste("otu_num_", input$otu_number, ".csv", sep="")
    },
    content <- function(file) {
      subset_ids <- sequence_ids[sequencecluster_labels==input$otu_number]
      subset_table <- timeseries_table[sequencecluster_labels==input$otu_number, ]
      subset_table <- as.matrix(subset_table)
      if (dim(subset_table)[2] == 1) {
        subset_table <- t(subset_table)
      }
      tax_ids <- taxonomic_ids[sequencecluster_labels==input$otu_number]
      time_clusts <- values$current_clusters[sequencecluster_labels==input$otu_number]
      
      save_df <- data.frame(Abundance=rowSums(subset_table),
                   TimeClustNumber=time_clusts,
                   TaxonomicID=tax_ids,
                   SequenceID=subset_ids)
      time_series <- as.data.frame(subset_table)
      colnames(time_series) <- tp
      save_df <- cbind(save_df, time_series)
      write.csv(save_df, file)
    }
  )
  
  # Save OTU plot as SVG
  output$saveOTUPlot <- downloadHandler(
    filename <- function() {
      paste('otu_num_', input$otu_number, '.svg', sep="")
    },
    content <- function(file) {
      if (is.null(input$otu_number))
      return()
      # Take only the time-series that are in the current cluster
      subset_table <- timeseries_table[sequencecluster_labels==input$otu_number, ]
      subset_table <- as.matrix(subset_table)
      # Fix behaviour when there's only one row, it makes it
      # column-wise, so we have to transpose it
      if (dim(subset_table)[2] == 1) {
        subset_table <- t(subset_table)
      }
      # Plotting code for main OTU plot
      svg("otuplot.svg",height=plot_height,width=plot_width)
      layout(as.matrix(cbind(c(1,3),c(2,4))))
      plotWrapper(subset_table, otu_plot=TRUE)
      dev.off()
      file.copy("otuplot.svg", file)
    }
  )
  
  plotWrapper <- function(subset_table, otu_plot=FALSE) {
    # Make a normalized version
    col_norm_subset_table <- t(apply(subset_table,1,function(x) x/timeseries_totals))
    # If a column has been removed by filter step, normalization
    # returns NaNs and for some reason Inf
    # set it as zero to remove gaps in plot
    col_norm_subset_table[is.nan(col_norm_subset_table)] <- 0
    col_norm_subset_table[is.infinite(col_norm_subset_table)] <- 0
    double_norm <- col_norm_subset_table/rowSums(col_norm_subset_table)
    # Plotting code for main time-series plot
    p1 <- plotTimeSeriesMatrix(col_norm_subset_table, "Normalized by Reads per Time Point", "Time", "Relative Abundance", otu_plot)
    p2 <- plotTimeSeriesMatrix(double_norm, "Normalized by Reads Per Time Point and Within Sequence", "Time", "Relative Abundance", otu_plot)
    multiplot(p1, p2)
  }
  
  plotTimeSeriesMatrix<-function(mat, plot_title, x_label, y_label, otu_plot=FALSE) {
      #Add the time points, melt the matrix
      colnames(mat) <- tp
      mmat <- melt(t(mat))
      mmat$mask <- mask
      colnames(mmat) <- c("time","sequence","value","mask")
      if (otu_plot) {
          time_clusts <- values$current_clusters[sequencecluster_labels==input$otu_number]
          mmat$time_clust <- time_clusts[mmat$sequence]
          if (input$excludeNoise) {
              mmat <- mmat[mmat$time_clust != -1,]
              if (nrow(mmat) == 0) {
                stop("All sequences in this OTU were labelled by DBSCAN as noise. Nothing to plot.")
                return()
              }
          }
          p <- ggplot(mmat,aes(x=time,y=value,group=sequence,colour=time_clust))+geom_line(alpha=0.3, size=1.5)
      } else {
          p <- ggplot(mmat,aes(x=time,y=value,group=sequence))+geom_line(alpha=0.3, size=1.5)
      }
      if (length(unique(mask)) > 1) {
          p <- p+facet_wrap(~mask, nrow=1, scales="free")
      }
      p <- p+ggtitle(plot_title)+xlab(x_label)+ylab(y_label)
      return(p)
  }
  
})

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

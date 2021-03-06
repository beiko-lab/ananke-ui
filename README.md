# ananke-ui
Shiny-based interface for exploring ananke results.

###Screenshot:
<img src="http://kiwi.cs.dal.ca/~beikolab/AnankeUI.png" width="600">

### Pre-requisites
Required R libraries are: [shiny](https://cran.r-project.org/web/packages/shiny/index.html), [gtools](https://cran.r-project.org/web/packages/gtools/index.html), [shinyFiles](https://cran.r-project.org/web/packages/shinyFiles/index.html), [SparseM](https://cran.r-project.org/web/packages/SparseM/index.html), [ggplot2](https://cran.r-project.org/web/packages/SparseM/index.html), [gtools](https://cran.r-project.org/web/packages/gtools/index.html) and [rhdf5](http://bioconductor.org/packages/release/bioc/html/rhdf5.html). With the exception of rhdf5, these can be installed from CRAN with the usual `install.packages()` command. rhdf5 can be installed via Bioconductor:
```
install.packages(c("shiny","shinyFiles","SparseM","ggplot2","gtools"))
source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5")
```

A .h5 file generated by timeclust is required as input. More information on this is available in the [Ananke GitHub repository](https://github.com/beiko-lab/ananke/).

### Launching the UI
The UI can be launched from an R session:
```
library(shiny)
runGitHub("beiko-lab/ananke-ui")
```

### Loading A Data File
On the initial page, click the "Time-series database file" button and use the file selector to locate your .h5 database file created by ananke.

Once your file is selected, click the "Load data!" button. It may take a minute to load the file. Progress will be indicated in the top-right corner of the page. **Known issue:** Do not navigate to another pane until the file has finished loading, or the pages will not render properly.

### Data Summary
On this pane, there is an option to plot the number of clusters versus the epsilon clustering parameter. This can help decide where to begin exploring your data. If *ananke* was set to cluster over a large range of epsilon values, this operation could take a significant amount of time.

A good parameter for beginning exploration is the epsilon that returns the largest number of clusters. This is shown as a peak in the nclust vs epsilon parameters plot.

Once an epsilon parameter has been selected on the "Explore Clusters" pane, the Data Summary pane will have an option to plot the taxonomic consistency (via Simpson's index) of the time clusters.

### Explore Clusters
This pane allows you to explore your sequence clusters from two angles: by time clusters and by OTU clusters (if that data has been input into the database file). 

###### Exploring Time Clusters
Cluster IDs can be selected from the dropdowns, and are sorted by total sequence abundance. For time clusters, the clustering parameter epsilon can be changed using the spinner input. As this value is increased, the time series clusters become more "coarse" and will contain more dissimilar time series.

Four plots based on the time series cluster are generated. In each plot, each line represents a unique sequence as it is tracked across time. The top-left plot contains the time series plots of the raw sequence counts. In the top-right, counts have been normalized by sequence depth. In the bottom-right, counts have been normalized only within each sequence's time series (the counts for each sequence is divided by the total count of the sequence). In the bottom-right, the sequences have been doubly normalized by sequence depth and within each sequence's time series.

Below the plots is a data table containing the taxonomic information and sequence identity cluster information (if provided), as well as the abundance of the sequence and the sequence hash. The sequence hash is unique and allows a sequence of interest to be tracked down in the unique sequence output file from ananke.

Both the time series plot and the data table can be saved using links in the left sidepanel. Plots are saved as SVG files, and data tables as CSV files. The data table contains the information present in the table, as well as the sequence counts at each time point (i.e., you can download this table and plot your time series in Excel if that's what you're into).

###### Exploring OTUs
OTUs are presented in an identical way to the time series clusters. The time series plots are coloured by time clusters (at the epsilon currently selected in the time cluster tab). This allows the user to spot OTUs with temporally discordant sequences at a glance. OTUs with multiple line colours contain sequences which demonstrate different temporal dynamics. The data table lists the time cluster IDs, allowing you to explore these on the previous tab (you can type them into the cluster ID dropdown to get to them quickly).

Example of an OTU displaying multiple temporal dynamics. Lines are coloured according to time-series cluster.
<img src="http://kiwi.cs.dal.ca/~beikolab/DiscordantOTU.png" width="600">

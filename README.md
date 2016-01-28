### Description

__SC3__ (<b>S</b>ingle-<b>C</b>ell <b>C</b>onsensus <b>C</b>lustering) is an interactive tool for the unsupervised clustering of cells from single cell RNA-Seq experiments.

The latest stable version of __SC3__ can be installed from BioConductor:  
http://bioconductor.org/packages/devel/bioc/html/SC3.html

This repository is a developmental version of __SC3__. To install the developmental version please use:

```{R}
install.packages("devtools")
devtools::install_github("hemberg-lab/SC3")
library(SC3)
RSelenium::checkForServer()
```

__SC3__ manuscript is under review, but also available on bioRxiv:  
[http://biorxiv.org/content/early/2016/01/13/036558](http://biorxiv.org/content/early/2016/01/13/036558)

Please report any bugs, comments, issues or suggestions here:  
[https://github.com/hemberg-lab/SC3/issues](https://github.com/hemberg-lab/SC3/issues)

If you have any other questions please contact [Vladimir Kiselev](mailto:vk6@sanger.ac.uk).

### Test run

To test that the package has been installed successfully please run the following command:

```{R}
library(SC3)
sc3(treutlein, ks = 3:7, cell.filter = TRUE)
```

It should open __SC3__ in a browser window.

### "Built-in" datasets

There is one built-in dataset that is automatically loaded with __SC3__:

| Dataset | Source | __N__ cells | __k__ clusters |
--- | --- | --- | --- |
| [treutlein](http://www.nature.com/nature/journal/v509/n7500/full/nature13173.html) | Distal lung epithelium | 80 | 5 |

One can explore clusterings of this dataset by running the following commands (__ks__ parameter defines a region of __k__ needed to be investigated - see the next section):

```{R}
sc3(treutlein, ks = 3:7)
sc3(deng, ks = 8:12)
```

### Running __SC3__

If you would like to check the clustering of your own __dataset__ for __k__ (number of clusters) from 2 to 5, then you need to run the following:

```{R}
sc3(dataset, ks = 2:5)                        # without filtering of lowly expressed cells
sc3(dataset, ks = 2:5, cell.filter = TRUE)    # with filtering of lowly expressed cells
sc3(dataset, ks = 2:5, interactivity = FALSE) # without interactive visualisation
```

For more details please read the documentation by typing ?sc3

### Input file format

To run __SC3__ on an input file containing an expression matrix one need to preprocess the input file so that it looks as follows:


|  | cell1 | cell2 | cell3 | cell4 | cell5 
--- | --- | --- | --- | --- | ---
| __gene1__ | 1 | 2 | 3 | 4 | 5 
| __gene2__ | 1 | 2 | 3 | 4 | 5 
| __gene3__ | 1 | 2 | 3 | 4 | 5 


The first row of the expression matrix (with cell labels, e.g. __cell1__, __cell2__, etc.) should contain one fewer field than all other rows. Separators should be either spaces or tabs. If separators are commas (,) then the extension of the file must be .csv. If a path to your input file is "/path/to/input/file/expression-matrix.txt", to run it:

```{R}
sc3("/path/to/input/file/expression-matrix.txt", ks = 2:5)
```

### License

GPL-3

## specific the columns of names 
library(plyr)

suppressPackageStartupMessages(require(argparse)) # don't say "Loading required package: optparse"

parser <- ArgumentParser(description='Process some files into Sankey plots')


parser$add_argument('--name_col', dest='name_col', type="integer",
  default=1,
  help='Which row store the names')

parser$add_argument('--label_col', dest='label_col', type="integer",
  default=2,
  help='Which row store the labels')

parser$add_argument('files', metavar='FILE', type="character", nargs='+',
  help='Files that contain the cluster assignment')

#parser$print_help()
args = parser$parse_args()
#args <- parser$parse_args(c("--label_col", "3","heart.R5.statH","heart.R6.statH","heart.R7.statH"))
#print(commandArgs(TRUE))
print(args)

file_list = list()
for ( file in args$files) {
  tmp = read.delim(file)[,c(args$name_col,args$label_col)]
  file_list[[file]] = tmp
}

source = NULL
target = NULL
value  = NULL

for (idx in 1:(length(file_list)-1)){
dat = merge(file_list[[idx]],file_list[[idx+1]],by=1)
trans = count(dat[,-1])

source = c(source,paste(idx,trans[,1],sep=".C"))
target = c(target,paste(idx+1,trans[,2],sep=".C"))
value = c(value,trans[,3])
}

label = unique(c(source,target))
source2 = match(source,label)-1
target2 = match(target,label)-1


library(plotly)
p <- plot_ly(
    type = "sankey",
    orientation = "h",

    node = list(
      label =  label, #c("A1", "A2", "B1", "B2", "C1", "C2"),
#      color = c("blue", "blue", "blue", "blue", "blue", "blue"),
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),

    link = list(
      source = source2, #c(0,1,0,2,3,3),
      target = target2, #c(2,3,3,4,4,5),
      value =  value  #c(8,4,2,8,4,2)
    )
  ) %>% 
  layout(
    title = "Basic Sankey Diagram",
    font = list(
      size = 10
    )
)

htmlwidgets::saveWidget(as_widget(p), "sankey_plot.html")


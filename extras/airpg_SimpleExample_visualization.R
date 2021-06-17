#!/usr/bin/R
## author: Michael Gruenstaeudl, email: m.gruenstaeudl@fu-berlin.de

########################################################################

library(ggplot2)
library(tcltk)
library(tools)
library(dplyr)
library(forcats)

########################################################################

start_date = "2003-01-01"
end_date = "2021-12-31"

########################################################################

AvailTableFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
AvailTableData = read.csv(AvailTableFn, sep = "\t")
IRTableFn = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
out_fn = dirname(IRTableFn)
IRTableData = read.csv(IRTableFn, sep = "\t")
combinedDF = merge(AvailTableData, IRTableData, by="ACCESSION")

########################################################################

combinedDF = transform(combinedDF,
IRa_REPORTED_START = as.integer(as.character(IRa_REPORTED_START)),
IRb_REPORTED_START = as.integer(as.character(IRb_REPORTED_START)),
IRa_REPORTED_END = as.integer(as.character(IRa_REPORTED_END)),
IRb_REPORTED_END = as.integer(as.character(IRb_REPORTED_END)),
IRa_REPORTED_LENGTH = as.integer(as.character(IRa_REPORTED_LENGTH)),
IRb_REPORTED_LENGTH = as.integer(as.character(IRb_REPORTED_LENGTH))
)

########################################################################

posMatch = combinedDF[which(combinedDF[,"IRa_REPORTED"]=="yes"),]
negMatch = combinedDF[which(combinedDF[,"IRa_REPORTED"]=="no"),]
posDatesData = as.Date(posMatch$CREATE_DATE, format="%Y-%m-%d")
negDatesData = as.Date(negMatch$CREATE_DATE, format="%Y-%m-%d")
posTab = table(cut(posDatesData, 'year'))
negTab = table(cut(negDatesData, 'year'))
posPlotData = data.frame(DATE=as.Date(names(posTab)), FREQ_RECORDS=as.vector(posTab), CRITERION="positive")
negPlotData = data.frame(DATE=as.Date(names(negTab)), FREQ_RECORDS=as.vector(negTab), CRITERION="negative")
plotData = rbind(posPlotData, negPlotData)
plotData = plotData[order(plotData$DATE),]
plotData = plotData %>% group_by(CRITERION) %>% mutate(CUMFREQ=cumsum(FREQ_RECORDS), .keep = "all")

########################################################################

base_plot = ggplot(data=plotData, aes(x=DATE, y=CUMFREQ, fill=forcats::fct_rev(CRITERION)), width=1) +  #forcats::fct_rev inverts the order
    geom_bar(stat="identity", position="stack", alpha=0.5)

myPlot = base_plot + 
    xlab("Year") + 
    ylab("Cumulative Number of Records\n") + 
    scale_x_date(
        limits=c(as.Date(start_date), as.Date(end_date)),
        date_breaks="1 year",
        minor_breaks=NULL,
        expand=expansion(0),
        date_labels="%Y"
    ) + 
    scale_fill_manual(values=c("grey50", "grey0"), name="Complete IR annotations", labels=c("Yes", "No")) + 
    theme_minimal() + 
    theme(plot.title = element_text(size=16),
          plot.subtitle = element_text(size=16, face="italic"),
          axis.text=element_text(size=16),
          axis.text.x = element_text(size=16, angle=45, hjust=0, vjust=0.5),
          axis.title=element_text(size=18, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),
          legend.key.width=unit(1,"cm"),
          legend.position = "bottom",
          legend.text = element_text(size=14),
          legend.title = element_text(size=14))

ggsave(file = "./airpg_SimpleExample_visualization.pdf", plot=myPlot)

########################################################################



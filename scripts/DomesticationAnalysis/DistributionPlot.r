library(ggplot2)
library(patchwork)

argv<-commandArgs(TRUE)   

data <- read.table(argv[1], header = TRUE, sep="\t")

data$Chr <- factor(data$Chr, levels=c("C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09"))
chrColor <- c("#FFFFFF", "#9933FF", "#FF33FF", "#FF3399", "#FF3333", "#FF9933", "#99FF33", "#33FF99", "#3399FF")

p1 <- ggplot(data, aes(x=(Start+End)/2/1000000, y=Count, color=Chr, fill=Chr)) +
      geom_point(size=0.8, shape=16, alpha=0.8) +
      facet_grid(. ~ Chr, as.table=TRUE, scales="free_x", space="free_x")  +
      geom_hline(aes(yintercept=187), colour="red", linetype="dashed")  +       ##- top 5%
      geom_hline(aes(yintercept=187), colour="blue", linetype="dashed")  +      ##- top 10%
      theme(panel.background = element_blank(),axis.line = element_line(), axis.text.x=element_text(size=6, angle=0), legend.key.size = unit(0.4,'cm'), legend.text = element_text(size=6,angle=0), legend.position = 'none') + 
      xlab("Chromosome") + ylab("Candidate k-mers")



outpdf <-  paste(argv[1], ".distribution_density.pdf", sep="")


#mergeP <- p1 + p2 + plot_layout(widths = c(5, 1))


ggsave(file=outpdf, plot=p1, width = 7, height = 2.8)



## * load R packages
library(data.table)
library(ggplot2)

## * Load data
df.sim <- read.table("resSim.txt") ## load here the file containing the simulation results

## * Reshape results to the long format
## probably need to be updated 
dtL <- melt(as.data.table(df.sim), id.vars = c("n","cens","rho_1","rho_2"),
            value.name = "rEfficiency", variable.name = "estimator")
## dtL[, rho := paste0("\rho[1]=",rho_1," \rho[2]=",rho_2,"")]
dtL[, rho := paste0("rho[1]==",rho_1,"~~~~rho[2]==",rho_2)]
dtL[, censoring := paste0(cens,"*\'%\'~","censoring")]

## * Graphical display
figure1 <- ggplot(dtL, aes(x=n, y = rEfficiency, group = estimator, color = estimator, shape = estimator))
figure1 <- figure1 + geom_abline(intercept = 1, slope = 0, color = "red")
figure1 <- figure1 + geom_line(size=1.25) + geom_point(size=2.5)
figure1 <- figure1 + ylab("Efficiency relative to the Fine-Gray estimator") + xlab("Sample size")
figure1 <- figure1 + facet_grid(censoring~rho, labeller = label_parsed)
figure1 <- figure1 + scale_shape_manual(values = rep(16:17,3))
figure1 <- figure1 + scale_color_manual(values = c("red","orange","blue","purple","forestgreen","darkolivegreen2"))
figure1 <- figure1 + theme(text = element_text(size=15), legend.position = "bottom")
figure1


## * Export as pdf
## ggsave(figure1, filename = "figure1.pdf")




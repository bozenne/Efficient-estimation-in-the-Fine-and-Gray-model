## * load R packages
library(mets) ## version >= 1.2.9

## * Load data
if(file.exists("kumar.txt")){ ## load original data
    kumar <- read.table("kumar.txt")
}else{ ## simulated data
    set.seed(10)
    n <- 864
    kumar <- kumarsim(n,depcens=1)
}

## * Process data
kumar$cause <- kumar$status
kumar$ttt24 <- kumar[,6]

## * Run estimation procedures

## relapse FG-model
c2 <- cifreg(Event(time,cause) ~ gp + dnr + preauto + ttt24 , data = kumar, propodds = NULL, cause = 2)
## relapse FG-model with stratified Kaplan-Meier
cc2 <- cifreg(Event(time,cause)~ gp + dnr + preauto + ttt24, data = kumar, propodds = NULL, cause = 2, cens.model =~ strata(gp,dnr,preauto))

## augmenting using Aalen-Johansen for working models, iteratively until convergence, here two steps
fgaugS <- FG_AugmentCifstrata(Event(time,cause) ~ gp + dnr + preauto + ttt24 + strata(gp,dnr,preauto,ttt24),
                              data = kumar, cause = 2,
                              E = c2$E)
fgaugS2 <- FG_AugmentCifstrata(Event(time,cause) ~ gp + dnr + preauto + ttt24 + strata(gp,dnr,preauto,ttt24),
                               data = kumar, cause = 2, E = fgaugS$E)


## augmenting using Aalen-Johansen for working models and with censoring weights, iteratively until convergence 
fgauglS <- FG_AugmentCifstrata(Event(time,cause) ~ gp + dnr + preauto + ttt24 + strataC(gp,dnr,preauto) + strata(gp,dnr,preauto,ttt24), data = kumar , cause = 2, E = c2$E)
fgauglS2 <- FG_AugmentCifstrata(Event(time,cause) ~ gp + dnr + preauto + ttt24 + strataC(gp,dnr,preauto) + strata(gp,dnr,preauto,ttt24), data = kumar, cause = 2, E = fgauglS$E)

## * Gather results into a table
table1 <- cbind(
    round(summary(c2)$coef[,1:2],3),
    round(summary(fgaugS2)$coef[,1:2],3),
    round(summary(cc2)$coef[,1:2],3),
    round(summary(fgauglS2)$coef[,1:2],3)
)
rownames(table1) <- c("gp","dnr","preauto","ttt24")

table1 <- rbind(colnames(table1),table1)
colnames(table1) <- c("FG","","AUG","","FG-CM","","AUG-CM","")
print(table1, quote = FALSE)



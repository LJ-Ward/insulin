############################################################
## TWO-GROUP ROC
############################################################

library(pROC)

############################################
## LOAD
############################################

df <- read.csv("df1.csv",
               sep=";",
               stringsAsFactors=FALSE)

df$trueID <- factor(df$trueID,
                    levels=c("control","insulin"))

df$insulin <- as.numeric(df$insulin_score)

############################################
## COLOUR 
############################################

cols <- c(insulin="#FF7F0E")

############################################
## EMPTY PLOT 
############################################

plot(NA, NA,
     xlim=c(0,1), ylim=c(0,1),
     xlab="1 - Specificity",
     ylab="Sensitivity",
     asp=1)

abline(0,1,lty=2,col="grey70")

############################################
## DRAW ROC 
############################################

r_insulin <- roc(df$trueID=="insulin",
                 df$insulin,
                 direction="<",
                 quiet=TRUE)

# STEP LINE
lines(1 - r_insulin$specificities,
      r_insulin$sensitivities,
      type="s",
      lwd=3,
      col=cols["insulin"])

# CI ribbon
ci <- ci.se(r_insulin, specificities=seq(0,1,length=100))
x  <- 1 - as.numeric(rownames(ci))

polygon(c(x, rev(x)),
        c(ci[,1], rev(ci[,3])),
        col=adjustcolor(cols["insulin"],0.18),
        border=NA)

############################################
## TEXT 
############################################

ciA <- ci.auc(r_insulin)

text(0.55, 0.12,
     paste0("insulin  AUC=",
            sprintf("%.2f",auc(r_insulin)),
            " (",
            sprintf("%.2f",ciA[1]),"–",
            sprintf("%.2f",ciA[3]),")"),
     col=cols["insulin"],
     pos=4)

############################################
## LEGEND
############################################

legend("bottomright",
       legend="insulin",
       col=cols["insulin"],
       lwd=3,
       bty="n")

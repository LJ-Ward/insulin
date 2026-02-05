############################################
## MULTI-CLASS ROC
############################################

library(pROC)

############################################
## LOAD
############################################

df <- read.csv("df1.csv",
               sep=";",
               stringsAsFactors=FALSE)

df$trueID <- factor(df$trueID,
                    levels=c("control","hyper","insulin"))

df$control <- as.numeric(df$control_score)
df$hyper   <- as.numeric(df$hyper_score)
df$insulin <- as.numeric(df$insulin_score)

############################################
## COLOURS
############################################

cols <- c(
  control="#1F77B4",
  hyper="#009E73",
  insulin="#FF7F0E"
)

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
## FUNCTION TO DRAW ONE ROC
############################################

drawROC <- function(truth, score, col){
  
  r <- roc(truth, score, direction="<", quiet=TRUE)
  
  # ROC line (STEP STYLE — clean)
  lines(1 - r$specificities,
        r$sensitivities,
        type="s",
        lwd=3,
        col=col)
  
  # CI ribbon (THIS is the only fix)
  ci <- ci.se(r, specificities=seq(0,1,length=100))
  
  x  <- 1 - as.numeric(rownames(ci))
  
  polygon(c(x, rev(x)),
          c(ci[,1], rev(ci[,3])),   # LOWER & UPPER ONLY
          col=adjustcolor(col,0.18),
          border=NA)
  
  return(r)
}

############################################
## DRAW ALL THREE
############################################

r_control <- drawROC(df$trueID=="control", df$control, cols["control"])
r_hyper   <- drawROC(df$trueID=="hyper",   df$hyper,   cols["hyper"])
r_insulin <- drawROC(df$trueID=="insulin", df$insulin, cols["insulin"])

############################################
## TEXT
############################################

addText <- function(r, label, col, y){
  ciA <- ci.auc(r)
  text(0.55, y,
       paste0(label,"  AUC=",
              sprintf("%.2f",auc(r)),
              " (",
              sprintf("%.2f",ciA[1]),"–",
              sprintf("%.2f",ciA[3]),")"),
       col=col, pos=4)
}

addText(r_control,"control",cols["control"],0.25)
addText(r_hyper,"hyper",cols["hyper"],0.18)
addText(r_insulin,"insulin",cols["insulin"],0.11)

############################################
## LEGEND
############################################

legend("bottomright",
       legend=c("control","hyper","insulin"),
       col=cols,
       lwd=3,
       bty="n")

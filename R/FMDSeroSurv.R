
###############################################################################
# Program for estimation of state-level proportions of DIVA positivity rates  #
# Estimation of national level DIVA-positivity rates FMDV infection           #
###############################################################################

#' This package is developed for the estimation of state & national level estimators and the related measures for sero-surveilance.
#' @param SeroSurvData NSP-Ab antibody n X 3 data frame (n: number of states in the sample) obtained from sero-surveillance, where row represents the states, first column represents the total sample collected from each state and third column represents the number of DIVA positive samples.
#' @param Census_Data Animal census N X 1 data frame (state wise bovine population), where rows are states and column is the bovine population (cattle + buffalo) (e.g. N = 28 before 2014 & N= 29 after 2014)
#'
#' @return
#' A data frame containing the state-level and national level estimates of DIVA positivity rates and related parameters of FMDV infection:
#' \itemize{
#'   \item1 DIVA_Positive, Var_DIVA_Prop., StError_Prop., MarginErr_Prop., CI_Lower, CI_Upper, Pred_Total are the DIVA positive rates (proportion as estimator), variance of the estimator, standard error, margin of error, 95 percentage CI and predicted total of animals having history FMD infection at the individual state-level.
#'   \item2 DIVA_Pos_Prop, Pos Rate, StError, CV, MarError, Lower_CI, Upper_CI, PredNatTotal, SE_PredTotal are the DIVA positive proportion, rates, its standard error, margin of error, co-efficient of variation, 95 percentage CI, predicted total, standard error of predicted total at the country/national level.
#'  }
#'
#'  @author Samarendra Das
#'
#'  @example
#'  set.seed(23)
#'  SeroSurvData <- data.frame(round(runif(5, 100, 200)), rnbinom(5, mu = 40, size = 1), row.names = c("AP", "KN", "MH", "OD", "Guj"))
#'  colnames(SeroSurvData) <- c("Samples", "Positive")
#'  Census_Data <- data.frame(round(runif(10, 400, 600)), row.names=c("AP", "KN", "MH", "OD", "Guj", "Goa", "Raj", "Assam", "Miz", "Meg"))
#'  colnames(Census_Data) <- "Total"
#'  #Do not run
#'  result <- FMDSeroSurv (SeroSurvData, Census_Data)
#'  state_estimate <- result$StateEstimates
#'  Nat_estimate <- result$NationalEstimates
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#'
#'

FMDSeroSurv <- function (SeroSurvData, Census_Data) {


  if(!is.matrix(SeroSurvData) & !is.data.frame(SeroSurvData) & class(SeroSurvData)[1] != "dgCMatrix")
    stop("Wrong input data type of 'SeroSurvData'")
  if(sum(is.na(SeroSurvData)) > 0)
    stop("NAs detected in input 'SeroSurvData'");gc();
  if(sum(SeroSurvData < 0) > 0)
    stop("Negative values detected in input 'SeroSurvData'");gc();
  if(all(SeroSurvData == 0))
    stop("All elements of input 'SeroSurvData' are zeros");gc();
  if(any(colSums(SeroSurvData) == 0))
    warning("Library size of zeros detected in 'SeroSurveillance Data'");gc();

  if(!is.matrix(Census_Data) & !is.data.frame(Census_Data) & class(Census_Data)[1] != "dgCMatrix")
    stop("Wrong input data type of 'Census_Data'")
  if(sum(is.na(Census_Data)) > 0)
    warning("NAs detected in input 'Census_Data'");gc();
  if(sum(na.omit(Census_Data) < 0) > 0)
    stop("Negative values detected in input 'Census_Data'");gc();
  if(all(Census_Data == 0))
    stop("All elements of input 'Census_Data' are zeros");gc();

  if (nrow(Census_Data) < nrow(SeroSurvData))
    stop("No. of states in sample must be less than the total states")

###inputs from SeroSurvData
#N <- 29                 ####total number of states
N <- nrow(na.omit(Census_Data))
n <- nrow (SeroSurvData)         #######number of states considered in sampling

mi <- as.numeric(SeroSurvData[,1])  #####number random samples from states
M.all <- as.numeric(Census_Data [, 1])     # state wise animal census data
M.all <-na.omit (M.all)
M.total <- sum(M.all)           ####total bovine population in the country
Mi <- Census_Data[match(row.names(SeroSurvData), row.names(Census_Data)), ]
diva.pos <- as.numeric(SeroSurvData[, 2])
##########Estimator
pi <- diva.pos / mi                ######state wise DIVA positivity rates (sucesses)
pi.perc <- round(pi * 100, 2)
qi <- 1 - pi                    ####rates of failure

stat_diva.tot <- round(Mi * pi)        #####state level numbers for DIVA tests

samp_var <- (pi * qi) / (mi - 1)         #######sample variance
Est_Var_pi <- (Mi - mi) / (Mi * (mi - 1)) * pi * qi  ###(Est value) Variance of the estimator

stand_err.pi <- abs(sqrt(Est_Var_pi))     ##### Standard error of the estimator
mar_err.pi <- 1.96 * stand_err.pi         #####margin error of the estimator
up_CI.pi <- pi + mar_err.pi
low_CI.pi <- pi - mar_err.pi

out1 <- cbind(pi.perc, signif(Est_Var_pi, 3), round(stand_err.pi, 3), round(mar_err.pi *100, 2), round(low_CI.pi, 2), round(up_CI.pi, 2), stat_diva.tot)
#row.names(out1) <-  row.names(SeroSurvData)
colnam <- c("DIVA_Positive(%)", "Var_DIVA_Prop.", "StError_Prop.", "MarginErr_Prop.(%)", "95%_CI_Low", "95%_CI_Upper", "Pred_Total")
StateEstimates <- data.frame(out1, row.names = row.names(SeroSurvData))
rm(out1)
colnames(StateEstimates) <- colnam
#print(out1)
#class(out1) <- "State Level Estimates"

####National level estimates
nat.diva.tot <- round ((N / n) * sum (stat_diva.tot))     #####est. value national diva total
nat.prop <- nat.diva.tot / M.total
Perc <- round(nat.prop *100, 2)

bet_var <- 1 / (n - 1) * sum((stat_diva.tot - sum(stat_diva.tot) / n) ^ 2)
within_var <- sum ((Mi * (Mi - mi) / mi) * samp_var)

var.nat.tot <- (N * (N - n) / n) * bet_var + (N / n) * within_var
SE.nat.tot <- round(abs(sqrt(var.nat.tot)), 2)

var.nat.prop <- var.nat.tot / M.total^2
SE.nat.prop <- abs(sqrt(var.nat.prop))
CV.nat.prop <- round(SE.nat.prop / nat.prop * 100, 2)
ME.nat.prop <- round(1.96 * SE.nat.prop * 100, 2)
up.CI.nat.prop <- round(nat.prop + 1.96 * SE.nat.prop , 2)
low.CI.nat.prop <- round(nat.prop - 1.96 * SE.nat.prop , 2)

NationalEstimates <- c(round(nat.prop, 3), Perc, round(SE.nat.prop, 3), CV.nat.prop, ME.nat.prop,  low.CI.nat.prop, up.CI.nat.prop, round(nat.diva.tot), SE.nat.tot)
names(NationalEstimates) <- c("DIVA_Pos_Prop", "Pos Rate(%)", "StError", "CV(%)", "MarError(%)", "Low_95%_CI", "Upper_95%_CI", "PredNatTotal", "SE_PredTotal")

finOut <- list(StateEstimates=StateEstimates, NationalEstimates=NationalEstimates)
return(finOut)
}

###########################Ends here ####################################

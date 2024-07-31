#' @title Uniform Error Index
#' @param actual Univariate Data Series
#' @param predicted Predicted Data Series from Models
#' @import FactoMineR factoextra Metrics
#' @return
#' \itemize{
#'\item ErrorSeries: Uniform Error Index Series
#'\item ErrorMetrics: Values of Different Error Measures i.e.,Relative Absolute Error (RAE),Mean Absolute Error(MAE),Median Absolute Error (MDAE),Mean Absolute Percent Error (MAPE),Root Mean Squared Error (RMSE),Mean Squared Error (MSE),Symmetric Mean Absolute Percentage Error(SMAPE), Sum of Squared Errors (SSE),Mean Uniform Error Index (MUEI).
#' }
#' @export
#' @examples
#' library("UEI")
#' actual<- as.ts(rnorm(50,100,50))
#' predicted<- as.ts(rnorm(50,110,60))
#' Result <- UEI(actual, predicted)
#' @references
#' \itemize{
#'\item Yeasin, M. and Paul, R.K., 2024. OptiSembleForecasting: optimization-based ensemble forecasting using MCS algorithm and PCA-based error index. The Journal of Supercomputing, 80(2), pp.1568-1597.
#' }
UEI<-function(actual,predicted){
  #loss measure
  l1 = (predicted - actual)^2
  l2 = (predicted^2 - actual^2)^2
  l3 = log(predicted^2) + actual^2 * predicted^(-2)
  l4 = (log(actual^2 * predicted^-2))^2
  l5 = abs(predicted - actual)
  l6 = abs(predicted^2 - actual^2)
  l7 = abs((actual-predicted))/((actual+predicted)/2)
  l8 = ((actual-predicted)^2)/sum((actual-mean(actual))^2)
  l9 = abs((actual-predicted))/(abs(actual))
  l10= abs(actual-predicted)/(abs(actual)+abs(predicted))
  l_bind<-cbind(l1,l2,l3, l4, l5, l6, l7, l8, l9, l10)
  l_bind<-as.data.frame(l_bind)
  colnames(l_bind)<-paste0("Err", seq(1:length(l_bind)))
  loss_all<-  l_bind
  #normalisation
  norm_loss<-NULL
  minMax <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  norm_loss_1<-lapply(loss_all, minMax)
  norm_loss<-as.data.frame(norm_loss_1)
  res_pca <- PCA(norm_loss,scale.unit= TRUE,graph = FALSE)
  eig_val <- get_eigenvalue(res_pca)
  no_pca<-sum(eig_val[,1]>= 1)
  var <- get_pca_var(res_pca)
  var_cont<-fviz_contrib(res_pca, choice = "var", axes = c(seq(1,no_pca,1)))
  Weight<-as.matrix(var_cont$data[,2]/100)
  final_loss<-as.matrix(norm_loss)%*%(as.matrix(Weight))
  error_series<-final_loss
  #Different error measure
  RAE<-rae(actual,predicted)
  MAE<-mae(actual,predicted)
  MDAE<-mdae(actual,predicted)
  MAPE<-mape(actual,predicted)
  MSE<-mse(actual,predicted)
  RMSE<-rmse(actual,predicted)
  SMAPE<-smape(actual,predicted)
  SSE<- sse(actual,predicted)
  MUEI<-mean(error_series)
  Error<-cbind(MAE,MDAE,MAPE,RMSE,MSE,SMAPE,MUEI)
  Result<-list(ErrorSeries=error_series,ErrorMetrics=Error)
  return(Result)
}

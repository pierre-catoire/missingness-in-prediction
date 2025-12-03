GeneratePerformanceMetricsTables = function(performanceMetrics,
                                            methodsKeys = c("PS","MARG","MARGMI","UI","SCI","SCIMI","MI","MIMI"),
                                            metrics = c("MSPE-OMU","MSPE-OMC","MSE"),
                                            path = "results/tables/"){
  for(model in names(performanceMetrics)){
    for(metric in c("MSPE-OMU","MSPE-OMC","MSE")){
      tab = data.frame(MISSPROP = performanceMetrics[[model]][["MISSPROP"]])
      for(methodKey in methodsKeys){
        tab[[methodKey]] = performanceMetrics[[model]][[metric]][[methodKey]]
      }
      write.csv(tab,
                file = paste("results/tables/pm_continuous_",model,"_",metric,".csv", sep = ""),
                row.names = FALSE)
    }
  }
}

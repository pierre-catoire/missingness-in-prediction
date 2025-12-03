ms = function(x,ref,na.rm=T){mean((x-ref)^2, na.rm = na.rm)}

GeneratePerformanceMetricsTables = function(performanceMetrics,
                                            methodsKeys = c("PS","CCS","MARG","UI","SCI","SCIMI","MI","MIMI"),
                                            metrics = c("MSPE_OMU","MSPE_OMC","MSE"),
                                            path = "results/tables/",
                                            writeTab = T,
                                            span = 1){
  result = list()
  for(model in names(performanceMetrics)){
    result[[model]] = list()
    for(metric in metrics){
      result[[model]][[metric]] = list()
      
      # --- create a table for MISSPROP and the performance metrics of each method
      tabPoints = data.frame(MISSPROP = performanceMetrics[[model]][["MISSPROP"]])
      tabLOESS = data.frame(MISSPROP = seq(0,0.7, by = 0.001))
      
      for(methodKey in c(methodsKeys,"PRAGMATIC_MU","PRAGMATIC_MC")){
        tabPoints[[methodKey]] = performanceMetrics[[model]][[metric]][[methodKey]]
        loessFit = loess(tabPoints[[methodKey]] ~ tabPoints[["MISSPROP"]],
                         control = loess.control(surface = "direct"),
                         span = span)
        tabLOESS[[methodKey]] = predict(loessFit,
                                        tabLOESS[["MISSPROP"]])
      }
      result[[model]][[metric]][["points"]] = tabPoints
      result[[model]][[metric]][["loess"]] = tabLOESS
      
      
      # --- save the tabs in separated csv files ---
      
      if(writeTab){
        write.csv(tabPoints,
                  file = paste("results/tables/pm_continuous_",model,"_",metric,"_points.csv", sep = ""),
                  row.names = FALSE)
        write.csv(tabLOESS,
                  file = paste("results/tables/pm_continuous_",model,"_",metric,"_LOESS.csv", sep = ""),
                  row.names = FALSE)
      }
    }
  }
  return(result)
}

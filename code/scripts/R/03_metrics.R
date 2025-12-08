ms = function(x,ref,na.rm=T){mean((x-ref)^2, na.rm = na.rm)}

computePerformanceMetrics = function(predictions,
                                     methodsKeys = c("PS","CCS","MARG","MARGMI","SCI","SCIMI","MI","MIMI"),
                                     references = c("PRAGMATIC_MU","PRAGMATIC_MC"),
                                     metrics = c("MSPE_OMU", "MSPE_OMC","MSE"),
                                     groups = c("ALL","COMPLETE","INCOMPLETE"),
                                     loessSeq = seq(0,0.7, by = 0.001),
                                     span = 1,
                                     writeTab = T){
  
  performanceMetrics = list("LOESSSEQ" = loessSeq)
  for(model in names(predictions)){
    performanceMetrics[[model]] = list("ALL" = list(),
                                       "COMPLETE" = list(),
                                       "INCOMPLETE" = list(),
                                       "MISSPROP" = numeric())
    for(group in groups){
      for(metric in metrics){
        performanceMetrics[[model]][[group]][[metric]] = list()
        for(methodKey in c(methodsKeys, references)){
          performanceMetrics[[model]][[group]][[metric]][[methodKey]] = list("POINTS" = numeric(),
                                                                             "LOESS" = numeric())
        }
      }
    }
  }
  
  for(model in names(predictions)){
    for(i in 1:length(predictions[[model]])){
      performanceMetrics[[model]][["MISSPROP"]][i] = mean(predictions[[model]][[i]][["M1"]])
      for(group in groups){
        if(group == "ALL"){
          id = rep(T,length(predictions[[model]][[i]][["M1"]]))
        }else if(group == "COMPLETE"){
          id = (predictions[[model]][[i]][["M1"]] == 0)
        }else if(group == "INCOMPLETE"){
          id = (predictions[[model]][[i]][["M1"]] == 1)
        }
        
        for(metric in c("MSPE_OMU", "MSPE_OMC","MSE")){
          if(metric == "MSPE_OMU"){
            refMetric = predictions[[model]][[i]][["ORACLE_MU"]][id]
          }else if(metric == "MSPE_OMC"){
            refMetric = predictions[[model]][[i]][["ORACLE_MC"]][id]
          }else if(metric == "MSE"){
            refMetric = predictions[[model]][[i]][["Y"]][id]
          }
          
          for(methodKey in c(methodsKeys, references)){
            message("Model: ", model, ", i:",i," Group: ", group, ", Metric: ", metric, ", MethodKey: ", methodKey)
            predMethod = predictions[[model]][[i]][[methodKey]][id]
            performanceMetrics[[model]][[group]][[metric]][[methodKey]][["POINTS"]][i] = ms(predMethod, refMetric)
          }
        }
      }
    }
    # Fit the LOESS regression
    for(group in groups){
      for(metric in metrics){
        if(writeTab){
          tabPoints = data.frame("MISSPROP" = performanceMetrics[[model]][["MISSPROP"]])
          tabLoess  = data.frame("MISSPROP" = loessSeq)
        }
        
        for(methodKey in c(methodsKeys, references)){
          loessFit = loess(performanceMetrics[[model]][[group]][[metric]][[methodKey]][["POINTS"]] ~
                             performanceMetrics[[model]][["MISSPROP"]],
                           control = loess.control(surface = "direct"),
                           span = span)
          loessPred = predict(loessFit, loessSeq)
          performanceMetrics[[model]][[group]][[metric]][[methodKey]][["LOESS"]] = loessPred
          if(writeTab){
            tabPoints[[methodKey]] = performanceMetrics[[model]][[group]][[metric]][[methodKey]][["POINTS"]]
            tabLoess[[methodKey]] = loessPred
          }
        }
        if(writeTab == T){
          write.csv(tabPoints, file = paste("outputs/continuous/tables/table",model,group,metric,"points.csv",sep = "_"))
          write.csv(tabLoess, paste("outputs/continuous/tables/table",model,group,metric,"loess.csv",sep = "_"))
        }
      }
    }
  }
  return(performanceMetrics)
}

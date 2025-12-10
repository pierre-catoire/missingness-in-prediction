ms = function(x,ref,na.rm=T){mean((x-ref)^2, na.rm = na.rm)}

computePerformanceMetrics = function(predictions,
                                     methodsKeys = c("PS","CCS","MARG","MARGMI","SCI","SCIMI","MI","MIMI"),
                                     references = c("PRAGMATIC_MU","PRAGMATIC_MC"),
                                     metrics = c("MSPE_OMU", "MSPE_OMC","MSE"),
                                     groups = c("ALL","COMPLETE","INCOMPLETE"),
                                     loessSeq = seq(0,0.7, by = 0.001),
                                     missInd = "M1",
                                     span = 1,
                                     writeTab = T,
                                     thresholdLOESS = 0.1){
  
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
      performanceMetrics[[model]][["MISSPROP"]][i] = 1-mean(as.numeric(predictions[[model]][[i]][[missInd]] == 0))
      for(group in groups){
        if(group == "ALL"){
          id = rep(T,length(predictions[[model]][[i]][[missInd]]))
        }else if(group == "COMPLETE"){
          id = (predictions[[model]][[i]][[missInd]] == 0)
        }else if(group == "INCOMPLETE"){
          id = (predictions[[model]][[i]][[missInd]] == 1)
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
          perfPoints = performanceMetrics[[model]][[group]][[metric]][[methodKey]][["POINTS"]]
          missProps = performanceMetrics[[model]][["MISSPROP"]]
          
          # --- don't fit LOESS on the lowest values for incomplete cases, as outliers produce 
          if(group == "INCOMPLETE"){
            slicer = missProps > thresholdLOESS
            missProps = missProps[slicer]
            perfPoints = perfPoints[slicer]
          }
          loessFit = loess(perfPoints ~ missProps,
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

#' Compute the performance metrics from a dataframe, as simulated_datasets
#' 
#' Used for datasets imported (notably for simulations with discrete predictors, imported from python scripts)

dataset2performanceMetrics = function(simulated_datasets,
                                      methodsKeys = c("PS","CCS","MARG","MARGMI","SCI","SCIMI","MI","MIMI"),
                                      identifier = "MISSCOEF",
                                      references = c("PRAGMATIC_MU","PRAGMATIC_MC"),
                                      metrics = c("MSPE_OMU", "MSPE_OMC","MSE"),
                                      groups = c("ALL","COMPLETE","INCOMPLETE"),
                                      loessSeq = seq(0,0.7, by = 0.001),
                                      missInd = "M1",
                                      span = 1,
                                      writeTab = F){
  
  performanceMetrics = list("LOESSSEQ" = loessSeq)
  for(model in unique(simulated_datasets[["MODEL"]])){
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
  
  for(model in unique(simulated_datasets[["MODEL"]])){
    i = 1
    separators = unique(simulated_datasets[simulated_datasets[["MODEL"]] == model,][[identifier]])
    
    for(sep in separators){
      dat = simulated_datasets[simulated_datasets[["MODEL"]] == model & simulated_datasets[[identifier]] == sep & simulated_datasets[["SET"]] == "TEST",]
      performanceMetrics[[model]][["MISSPROP"]][i] = 1-mean(as.numeric(dat[[missInd]] == 0))
      for(group in groups){
        if(group == "ALL"){
          id = rep(T,length(dat[[missInd]]))
        }else if(group == "COMPLETE"){
          id = (dat[[missInd]] == 0)
        }else if(group == "INCOMPLETE"){
          id = (dat[[missInd]] == 1)
        }
        
        for(metric in c("MSPE_OMU", "MSPE_OMC","MSE")){
          if(metric == "MSPE_OMU"){
            refMetric = dat[["ORACLE_MU"]][id]
          }else if(metric == "MSPE_OMC"){
            refMetric = dat[["ORACLE_MC"]][id]
          }else if(metric == "MSE"){
            refMetric = dat[["Y"]][id]
          }
          
          for(methodKey in c(methodsKeys, references)){
            predMethod = dat[[methodKey]][id]
            performanceMetrics[[model]][[group]][[metric]][[methodKey]][["POINTS"]][i] = ms(predMethod, refMetric)
          }
        }
      }
      i = i+1
    }
    # Fit the LOESS regression
    for(group in groups){
      for(metric in metrics){
        if(writeTab){
          tabPoints = data.frame("MISSPROP" = performanceMetrics[[model]][["MISSPROP"]])
          if(identifier == "MISSCOEF"){
            tabLoess  = data.frame("MISSPROP" = loessSeq)
          }
        }
        
        for(methodKey in c(methodsKeys, references)){
          if(identifier == "MISSCOEF"){
            perfPoints = performanceMetrics[[model]][[group]][[metric]][[methodKey]][["POINTS"]]
            missProps = performanceMetrics[[model]][["MISSPROP"]]
            
            # --- don't fit LOESS on the lowest values for incomplete cases, as outliers produce 
            if(group == "INCOMPLETE"){
              slicer = missProps > thresholdLOESS
              missProps = missProps[slicer]
              perfPoints = perfPoints[slicer]
            }
            loessFit = loess(perfPoints ~ missProps,
                             control = loess.control(surface = "direct"),
                             span = span)
            loessPred = predict(loessFit, loessSeq)
            performanceMetrics[[model]][[group]][[metric]][[methodKey]][["LOESS"]] = loessPred
            if(writeTab){
              tabPoints[[methodKey]] = performanceMetrics[[model]][[group]][[metric]][[methodKey]][["POINTS"]]
              tabLoess[[methodKey]] = loessPred
            }
          }
        }
        if(writeTab == T){
          write.csv(tabPoints, file = paste("outputs/continuous/tables/table",model,group,metric,"points.csv",sep = "_"))
          
          if(identifier == "MISSCOEF"){
            write.csv(tabLoess, paste("outputs/continuous/tables/table",model,group,metric,"loess.csv",sep = "_"))
          }
        }
      }
    }
  }
  return(performanceMetrics)
}

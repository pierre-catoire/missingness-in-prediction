load_package = function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

load_package("scales")
load_package("vioplot")

plotPerformanceMetrics = function(performanceMetrics,
                                  methodsKeys = c("SCI","MI","MARG","SCIMI","MIMI","UI","PS","CCS"),
                                  graphParameters = NULL,
                                  opacity = 1,
                                  ylims = list("MSPE_OMU" = list("ALL" = c(0,1.5),
                                                                 "SUBGROUPS" = c(0,3)),
                                               "MSPE_OMC" = list("ALL" = c(0,1.5),
                                                                 "SUBGROUPS" = c(0,3)),
                                               "MSE" = list("ALL" = c(0,5),
                                                            "SUBGROUPS" = c(0,10))),
                                  models = c("M1","M2","M3","M4","M5"),
                                  references = c("PRAGMATIC_MU","PRAGMATIC_MC"),
                                  metrics = c("MSPE_OMU","MSPE_OMC","MSE"),
                                  savePDF = F,
                                  storePath = "outputs/continuous/plots/plot"){
  
  if(is.null(graphParameters)){
    MC_colors = hue_pal(h = c(0, 90))(5)
    MU_colors = hue_pal(h = c(180, 270))(4)
    graphParameters = list("PS" = list("label" = "Pattern Submodels",
                                       "color" = MC_colors[1],
                                       "lty" = 1,
                                       "pch" = 16,
                                       "cex" = .5),
                           "CCS" = list("label" = "Complete Cases Submodels",
                                        "color" = MU_colors[1],
                                        "lty" = 1,
                                        "pch" = 16,
                                        "cex" = .5),
                           "MARG" = list("label" = "Marginalisation (without m. ind.)",
                                         "color" = MU_colors[2],
                                         "lty" = 1,
                                         "pch" = 16,
                                         "cex" = .5),
                           "MARGMI" = list("label" = "Marginalisation (with m. ind.)",
                                           "color" = MC_colors[2],
                                           "lty" = 1,
                                           "pch" = 16,
                                           "cex" = .5),
                           "UI" = list("label" = "Unconditional imputation",
                                       "color" = MC_colors[3],
                                       "lty" = 1,
                                       "pch" = 16,
                                       "cex" = .5),
                           "SCI" = list("label" = "Single conditional imputation (without m. ind.)",
                                        "color" = MU_colors[3],
                                        "lty" = 1,
                                        "pch" = 16,
                                        "cex" = .5),
                           "SCIMI" = list("label" = "Single conditional imputation (with m. ind.)",
                                          "color" = MC_colors[4],
                                          "lty" = 1,
                                          "pch" = 16,
                                          "cex" = .5),
                           "MI" = list("label" = "Multiple conditional imputation (without m. ind.)",
                                       "color" = MU_colors[4],
                                       "lty" = 1,
                                       "pch" = 16,
                                       "cex" = .5),
                           "MIMI" = list("label" = "Multiple conditional imputation (with m. ind.)",
                                         "color" = MC_colors[5],
                                         "lty" = 1,
                                         "pch" = 16,
                                         "cex" = .5),
                           "PRAGMATIC_MU" = list("label" = "MU probability (ref)",
                                                 "color" = "darkblue",
                                                 "lty" = 2,
                                                 "pch" = 16,
                                                 "cex" = 0),
                           "PRAGMATIC_MC" = list("label" = "MC probability (ref)",
                                                 "color" = "darkred",
                                                 "lty" = 2,
                                                 "pch" = 16,
                                                 "cex" = 0))
  }
  # Create vectors of methodLabels, colors, line types to pass to legend
  methodsLabels = c()
  cols = c()
  ltys = c()
  i = 1
  for(key in c(methodsKeys,references)){
    methodsLabels[i] = graphParameters[[key]][["label"]]
    names(methodsLabels)[i] = key
    cols[i] = graphParameters[[key]][["color"]]
    names(cols)[i] = key
    ltys[i] = graphParameters[[key]][["lty"]]
    names(ltys)[i] = key
    i=i+1
  }
  
  modelLabels = setNames(c("Model 1", "Model 2","Model 3", "Model 4", "Model 5"),
                         c("M1","M2","M3","M4","M5"))
  metricsLabels = setNames(c("Mean Squared Prediction Error (Oracle MU reference)",
                             "Mean Squared Predictions Error (Oracle MC reference)",
                             "Mean Squared Error"),
                           metrics)
  
  analysisGroupLabels = setNames(c("All cases","Complete cases","Incomplete cases"),
                                 c("ALL","COMPLETE","INCOMPLETE"))
  
  for(model in models){
    for(analysis in c("ALL","SUBGROUPS")){
      if(savePDF){
        pdf(paste(paste(storePath,model,analysis, sep = "_"),".pdf",sep = ""), width = 14, height = 7)
      }
      if(analysis == "ALL"){
        analysisGroups = "ALL"
        par(mfrow = c(1,length(metrics)))
      }
      else if(analysis == "SUBGROUPS"){
        par(mfrow = c(2,length(metrics)))
        analysisGroups = c("COMPLETE","INCOMPLETE")
      }
      for(analysisGroup in analysisGroups){
        for(metric in metrics){
          plot(1, type = "n",
               frame.plot = FALSE,
               xlab = "Missingness proportion",
               ylab = metricsLabels[metric],
               xlim = c(0, 0.7), 
               ylim = ylims[[metric]][[analysis]])
          title(main = paste(modelLabels[model], analysisGroupLabels[analysisGroup], sep = ", "))
          
          for(methodKey in c(methodsKeys,references)){
            points(performanceMetrics[[model]][["MISSPROP"]],
                   performanceMetrics[[model]][[analysisGroup]][[metric]][[methodKey]][["POINTS"]],
                   col = adjustcolor(graphParameters[[methodKey]][["color"]],
                                     alpha.f = opacity),
                   pch = graphParameters[[methodKey]][["pch"]],
                   cex = graphParameters[[methodKey]][["cex"]])
            lines(performanceMetrics[["LOESSSEQ"]],
                  performanceMetrics[[model]][[analysisGroup]][[metric]][[methodKey]][["LOESS"]],
                  col = graphParameters[[methodKey]][["color"]],
                  lty = graphParameters[[methodKey]][["lty"]])
          }
          legend("topleft",
                 bg="transparent",
                 bty = "n",
                 legend = methodsLabels[c(methodsKeys,references)],
                 lty = ltys,
                 col= cols)
        }
      }
      if(savePDF){dev.off()}
    }
  }
}

plotPerformanceMetricsM6M7 = function(performanceMetrics,
                                      methodsKeys = c("MARG","MARGMI","PS","MIA","CCS"),
                                      graphParameters = NULL,
                                      opacity = 1,
                                      ylims = list("MSPE_OMU" = list("ALL" = c(0,0.04),
                                                                     "SUBGROUPS" = c(0,0.08)),
                                                   "MSPE_OMC" = list("ALL" = c(0,0.04),
                                                                     "SUBGROUPS" = c(0,0.08)),
                                                   "MSE" = list("ALL" = c(0.15,0.3),
                                                                "SUBGROUPS" = c(0.15,0.3))),
                                      models = c("M6","M7"),
                                      references = c("PRAGMATIC_MU","PRAGMATIC_MC"),
                                      metrics = c("MSPE_OMU","MSPE_OMC","MSE"),
                                      savePDF = F,
                                      storePath = "outputs/discrete/plots/plot"){
  
  if(is.null(graphParameters)){
    MC_colors = hue_pal(h = c(0, 90))(6)
    MU_colors = hue_pal(h = c(180, 270))(4)
    graphParameters = list("PS" = list("label" = "Pattern Submodels",
                                       "color" = MC_colors[1],
                                       "lty" = 1,
                                       "pch" = 16,
                                       "cex" = .5),
                           "CCS" = list("label" = "Complete Cases Submodels",
                                        "color" = MU_colors[1],
                                        "lty" = 1,
                                        "pch" = 16,
                                        "cex" = .5),
                           "MARG" = list("label" = "Marginalisation (without m. ind.)",
                                         "color" = MU_colors[2],
                                         "lty" = 1,
                                         "pch" = 16,
                                         "cex" = .5),
                           "MARGMI" = list("label" = "Marginalisation (with m. ind.)",
                                           "color" = MC_colors[2],
                                           "lty" = 1,
                                           "pch" = 16,
                                           "cex" = .5),
                           "UI" = list("label" = "Unconditional imputation",
                                       "color" = MC_colors[3],
                                       "lty" = 1,
                                       "pch" = 16,
                                       "cex" = .5),
                           "SCI" = list("label" = "Single conditional imputation (without m. ind.)",
                                        "color" = MU_colors[3],
                                        "lty" = 1,
                                        "pch" = 16,
                                        "cex" = .5),
                           "SCIMI" = list("label" = "Single conditional imputation (with m. ind.)",
                                          "color" = MC_colors[4],
                                          "lty" = 1,
                                          "pch" = 16,
                                          "cex" = .5),
                           "MI" = list("label" = "Multiple conditional imputation (without m. ind.)",
                                       "color" = MU_colors[4],
                                       "lty" = 1,
                                       "pch" = 16,
                                       "cex" = .5),
                           "MIMI" = list("label" = "Multiple conditional imputation (with m. ind.)",
                                         "color" = MC_colors[5],
                                         "lty" = 1,
                                         "pch" = 16,
                                         "cex" = .5),
                           "MIA" = list("label" = "Missingness Incorporated as Attribute (MIA)",
                                        "color" = MC_colors[6],
                                        "lty" = 1,
                                        "pch" = 16,
                                        "cex" = .5),
                           "PRAGMATIC_MU" = list("label" = "MU probability (reference)",
                                                 "color" = "darkblue",
                                                 "lty" = 2,
                                                 "pch" = 16,
                                                 "cex" = 0),
                           "PRAGMATIC_MC" = list("label" = "MC probability (reference)",
                                                 "color" = "darkred",
                                                 "lty" = 2,
                                                 "pch" = 16,
                                                 "cex" = 0))
  }
  
  methodsLabels = c()
  cols = c()
  i = 1
  for(key in c(methodsKeys,references)){
    methodsLabels[i] = graphParameters[[key]][["label"]]
    names(methodsLabels)[i] = key
    cols[i] = adjustcolor(graphParameters[[key]][["color"]], alpha.f = opacity)
    names(cols)[i] = key
    i=i+1
  }
  
  modelLabels = setNames(c("Model 6", "Model 7"),
                         c("M6","M7"))
  metricsLabels = setNames(c("MSPE (Oracle MU reference)",
                             "MSPE (Oracle MC reference)",
                             "Mean Squared Error"),
                           metrics)
  
  analysisGroupLabels = setNames(c("All cases","Complete cases","Incomplete cases"),
                                 c("ALL","COMPLETE","INCOMPLETE"))
  
  for(model in models){
    for(analysis in c("ALL","SUBGROUPS")){
      if(savePDF){
        pdf(paste(paste(storePath,model,analysis, sep = "_"),".pdf",sep = ""), width = 14, height = 7)
      }
      if(analysis == "ALL"){
        analysisGroups = "ALL"
        par(mfrow = c(1,length(metrics)))
      }
      else if(analysis == "SUBGROUPS"){
        par(mfrow = c(2,length(metrics)))
        analysisGroups = c("COMPLETE","INCOMPLETE")
      }
      for(analysisGroup in analysisGroups){
        for(metric in metrics){
          dat = data.frame(matrix(NA,
                                  nrow = length(performanceMetrics[[model]][["MISSPROP"]]),
                                  ncol = length(methodsKeys)))
          colnames(dat) = methodsKeys
          for(methodKey in methodsKeys){
            dat[[methodKey]] = performanceMetrics[[model]][[analysisGroup]][[metric]][[methodKey]][["POINTS"]]
          }
          vioplot(dat,
                  frame.plot = FALSE,
                  col = cols,
                  ylim = ylims[[metric]][[analysis]])
          for(reference in references){
            abline(h = mean(performanceMetrics[[model]][[analysisGroup]][[metric]][[reference]][["POINTS"]]),
                   col = graphParameters[[reference]][["color"]],
                   lty = graphParameters[[reference]][["lty"]])
          }
          title(main = paste(paste(modelLabels[model], analysisGroupLabels[analysisGroup],sep = ", "),
                             "\n", metricsLabels[metric], sep = ""))
          legend("topleft",
                 bg="transparent",
                 bty = "n",
                 lty = c(2,2),
                 lwd = 2,
                 legend = methodsLabels[references],
                 col = cols[references])
        }
      }
      if(savePDF){dev.off()}
    }
  }
}




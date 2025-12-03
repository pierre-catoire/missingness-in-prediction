library(scales)

plotPerformanceMetrics = function(performanceMetricsTableList,
                                  methodsKeys = c("SCI","MI","MARG","SCIMI","MIMI","UI","PS","CCS"),
                                  graphParameters = NULL,
                                  opacity = .2){
  
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
                                                 "cex" = .5),
                           "PRAGMATIC_MC" = list("label" = "MC probability (ref)",
                                                 "color" = "darkred",
                                                 "lty" = 2,
                                                 "pch" = 16,
                                                 "cex" = .5))
  }
  
  
  
  references = c("PRAGMATIC_MU","PRAGMATIC_MC")
  ylims = list("MSPE_OMU" = c(0,1.5),
               "MSPE_OMC" = c(0,1.5),
               "MSE" = c(0,5))
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
  
  metricsLabels = setNames(c("Mean Squared Prediction Error (Oracle MU reference)",
                             "Mean Squared Predictions Error (Oracle MC reference)",
                             "Mean Squared Error"),
                           c("MSPE_OMU","MSPE_OMC","MSE"))
  
  for(model in names(performanceMetricsTableList)){
    par(mfrow = c(1,length(performanceMetricsTableList[[model]])))
    for(metric in names(performanceMetricsTableList[[model]])){
      plot(1, type = "n",
           frame.plot = FALSE,
           xlab = "Missingness proportion",
           ylab = metricsLabels[metric],
           xlim = c(0, 0.7), 
           ylim = ylims[[metric]],
           main = paste(model,", ", metric, sep = ""))
      for(methodKey in c(methodsKeys,references)){
        points(performanceMetricsTableList[[model]][[metric]][["points"]][["MISSPROP"]],
               performanceMetricsTableList[[model]][[metric]][["points"]][[methodKey]],
               col = adjustcolor(graphParameters[[methodKey]][["color"]],
                                 alpha.f = opacity),
               pch = graphParameters[[methodKey]][["pch"]],
               cex = graphParameters[[methodKey]][["cex"]])
        lines(performanceMetricsTableList[[model]][[metric]][["loess"]][["MISSPROP"]],
              performanceMetricsTableList[[model]][[metric]][["loess"]][[methodKey]],
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
}

# Issues:
# - CCS gives same as PS
# - MI and SCI give always MU ???
# - MIMI and SCIMI biased for MNAR
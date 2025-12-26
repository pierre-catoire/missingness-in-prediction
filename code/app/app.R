library(shiny)
library(shinyBS)
library(scales)
library(rsconnect)
#load("code/app/www/data/performanceMetrics.RData")
load("www/data/performanceMetrics.RData")
methods_available <- list(
  
  continuous = list(
    MU = c(
      "Marginalisation with EM (without miss. ind.)" = "MARG",
      "Multiple Conditional Imputation (without miss. ind.)" = "MI",
      "Single Conditional Imputation (without miss. ind.)" = "SCI",
      "Pragmatic MU reference" = "PRAGMATIC_MU"
    ),
    MC = c(
      "Pattern Submodels" = "PS",
      "Marginalisation with EM (with miss. ind.)" = "MARGMI",
      "Multiple Conditional Imputation (with miss. ind.)" = "MIMI",
      "Single Conditional Imputation (with miss. ind.)" = "SCIMI",
      "Pragmatic MC reference" = "PRAGMATIC_MC"
    ),
    OTHER = c("Complete Cases Submodels" = "CCS")
  ),
  
  discrete = list(
    MU = c(
      "Marginalisation with EM (without miss. ind.)",
      "Pragmatic MU reference"
    ),
    MC = c(
      "Pattern Submodels",
      "MIA",
      "Marginalisation with EM (with miss. ind.)",
      "Pragmatic MC reference"
    )
  )
)

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

ylims = list("MSPE_OMU" = list("ALL" = c(0,1.5),
                               "SUBGROUPS" = c(0,3)),
             "MSPE_OMC" = list("ALL" = c(0,1.5),
                               "SUBGROUPS" = c(0,3)),
             "MSE" = list("ALL" = c(0,5),
                          "SUBGROUPS" = c(0,10)))

metricsLabels = c("MSPE_OMU" = "Mean Squared Prediction Error (MU reference)",
                  "MSPE_OMC" = "Mean Squared Prediction Error (MC reference)",
                  "MSE" = "Mean Squared Error / Brier")
# --- UI ---

ui <- fluidPage(
  
  titlePanel("Prediction with Missing Information"),
  
  sidebarLayout(
    sidebarPanel(
      
      ## Level 1
      selectInput(
        "model",
        "Model:",
        choices = c("Model 1 (MCAR)" = "M1",
                    "Model 2 (MARX-YM)" = "M2",
                    "Model 3 (MNARX-YO, IMO, NICO)" = "M3",
                    "Model 4 (MNARX-YO, IMO, ICO)" = "M4",
                    "Model 5 (MARX-YO, IMO, ICO)" = "M5"),
        selected = 1
      ),
      
      tags$div(
        tags$strong("Predictive performance metric"),
        tags$a(
          id = "metric",
          icon("info-circle"),
          href = "#",
          style = "margin-left:6px;"
        ),
        bsPopover(
          id = "metric",
          title = "Predictive performance metric",
          content = "MSPE is Mean Squared difference between prediction and reference probability. MU reference is P(Y|X). MC reference is P(Y|X,M). Mean Squared Error/Brier is the mean squared difference between the prediction and the observed Y.",
          placement = "right",
          trigger = "click"
        )
      ),
      
      radioButtons(
        "metric",
        label = NULL,
        choices = c("Mean Squared Prediction Error (MU reference)" = "MSPE_OMU",
                    "Mean Squared Prediction Error (MC reference)"   = "MSPE_OMC",
                    "Mean Squared Error/Brier" = "MSE"),
        selected = NULL
      ),
      
      radioButtons(
        "var_type",
        "Variable type:",
        choices = c("Continuous" = "continuous",
                    "Discrete"   = "discrete"),
        selected = "continuous"
      ),
      
      radioButtons(
        "plot_elements",
        "Plot elements:",
        choices = c("Simulation result points" = "POINTS",
                    "LOESS regression line"   = "LINES",
                    "Both" = "BOTH"),
        selected = "BOTH"
      ),
      
      ## Level 2 — MU
      tags$div(
        tags$strong("Missingness-unconditional methods"),
        tags$a(
          id = "info_MU",
          icon("info-circle"),
          href = "#",
          style = "margin-left:6px;"
        ),
        bsPopover(
          id = "info_MU",
          title = "Missingness-unconditional methods",
          content = "Methods targeting Missingness-unconditional prediction, i.e. expected to be close to the MU reference. Pragmatic reference is P(Y,Xo(m)).",
          placement = "right",
          trigger = "click"
        )
      ),
      
      checkboxGroupInput(
        "methods_MU",
        label = NULL,
        choices = NULL
      ),
      
      ## Level 2 — MC
      tags$div(
        tags$strong("Missingness-conditional methods"),
        tags$a(
          id = "info_MC",
          icon("info-circle"),
          href = "#",
          style = "margin-left:6px;"
        ),
        bsPopover(
          id = "info_MC",
          title = "Missingness-conditional methods",
          content = "Methods targeting Missingness-conditional prediction, i.e. expected to be close to the MC reference. Pragmatic reference is P(Y,Xo(m),M=m).",
          placement = "right",
          trigger = "click"
        )
      ),
      
      checkboxGroupInput(
        "methods_MC",
        label = NULL,
        choices = NULL
      ),
      
      tags$div(
        tags$strong("Other methods"),
        tags$a(
          id = "info_OTHER",
          icon("info-circle"),
          href = "#",
          style = "margin-left:6px;"
        ),
        bsPopover(
          id = "info_OTHER",
          title = "Other methods",
          content = "Methods not targeting either Missingness-unconditional or missingness-conditional prediction.",
          placement = "right",
          trigger = "click"
        )
      ),
      
      checkboxGroupInput(
        "methods_OTHER",
        label = NULL,
        choices = NULL
      ),
      
    ),
    
    mainPanel(
      plotOutput("plot", height = "80vh"),
      
      hr(),
      
      fluidRow(
        column(6, uiOutput("plot_explanation")),
        column(6, uiOutput("graph_panel"))
      )
    )
  )
)


server <- function(input, output, session) {
  
  observeEvent(input$var_type, {
    
    updateCheckboxGroupInput(
      session,
      "methods_MU",
      choices = methods_available[[input$var_type]]$MU,
      selected = "Pragmatic MU reference"
    )
    
    updateCheckboxGroupInput(
      session,
      "methods_MC",
      choices = methods_available[[input$var_type]]$MC,
      selected = "Pragmatic MC reference"
    )
    
    updateCheckboxGroupInput(
      session,
      "methods_OTHER",
      choices = methods_available[[input$var_type]]$OTHER,
      selected = NULL
    )
    
  })
  
  output$plot <- renderPlot({
    plot(1, type = "n",
         frame.plot = FALSE,
         xlab = "Missingness proportion",
         ylab = metricsLabels[input$metric],
         xlim = c(0, 0.7), 
         ylim = ylims[[input$metric]][["ALL"]])
    if(input$plot_elements %in% c("POINTS","BOTH")){
      for(methodKey in c(input$methods_MU,input$methods_MC, input$methods_OTHER)){
        points(performanceMetrics[[input$model]][["MISSPROP"]],
               performanceMetrics[[input$model]][["ALL"]][[input$metric]][[methodKey]][["POINTS"]],
               pch = graphParameters[[methodKey]][["pch"]],
               cex = graphParameters[[methodKey]][["cex"]],
               col = adjustcolor(graphParameters[[methodKey]][["color"]],
                                 alpha.f = .25))
      }
    }
    if(input$plot_elements %in% c("LINES","BOTH")){
      for(methodKey in c(input$methods_MU,input$methods_MC, input$methods_OTHER)){
        lines(performanceMetrics[["LOESSSEQ"]],
              performanceMetrics[[input$model]][["ALL"]][[input$metric]][[methodKey]][["LOESS"]],
              lty = graphParameters[[methodKey]][["lty"]],
              col = graphParameters[[methodKey]][["color"]])
      }
    }
  })
  
  output$plot_explanation <- renderUI({
    
    withMathJax(
      bsCollapse(
        id = "plot_explanation",
        open = NULL,
        
        bsCollapsePanel(
          title = "Explanation of the graph",
          
          HTML("
        <p>
        This figure reports predictive performance, measured by
        $$\\text{MSPE}_{\\text{OMU}} = \\frac{1}{N} \\underset{i = 1}{\\overset{N}{\\sum}}\\left[(\\Pr(Y_i \\mid \\mathbf{X}_i) - \\hat{\\Pr}(Y_i))^2\\right].$$
        </p>
        
        <p>
        <strong>Properties of the model</strong> (see directed acyclic graph on the right):
        </p>

        <ul>
          <li>
            <strong>MCAR: </strong> the observation of the predictors is independent from the outcome and both observed and missing predictors. 
          </li>
          <li>
            <strong>MARX-YM: </strong> the observation of the predictors is independent from the outcome and the missing predictors, given the observed predictors. 
          </li>
          <li>
            <strong>MARX-YO: </strong> the observation of the predictors is independent fromthe missing predictors, given the outcome and the observed predictors. 
          </li>
          <li>
            <strong>NIMO: </strong> the outcome is independent from the observation of the predictors, given the value of the observed predictors
          </li>
          <li>
            <strong>NICO: </strong> the outcome is independent from the observation of the predictors, given the value of the observed predictors, when all predictors are observed.
          </li>
        </ul>
        <p>
        <strong>Displayed elements:</strong>
        </p>
        <p>
        <strong>Evaluated methods:</strong>
        <ul>
          <li>
            <i>Missingness-unconditioned methods: </i>
          </li>
          <ul>
          <li>
          Test
          </li>
          </ul>
        </ul>
        </p>
        ")
        )
      )
    )
  })
  
  output$graph_panel <- renderUI({
    
    bsCollapse(
      id = "graph_panel",
      open = NULL,
      
      bsCollapsePanel(
        title = "Graphical model",
        
        tags$div(
          style = "
      display: flex;
      justify-content: center;
    ",
          tags$embed(
            src = paste("img/",input$model,".png",sep = ""),
            style = "width:75%; max-width:75%;"
          )
        )
      )
    )
  })
  
}

shinyApp(ui, server)
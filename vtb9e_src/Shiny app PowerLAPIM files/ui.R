ui <- shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Power analysis for the longitudinal APIM to select the number of dyads in intensive longitudinal studies"),

  # Sidebar with controls to select the outputs to compute power
  sidebarPanel(
      shinyjs::useShinyjs(debug = TRUE),
      id = "side-panel",
      withMathJax(),
      width = 9,

    tags$style(
      "li a {
        font-size: 20px;
      }
    ",
    "p, div { 
      font-size: 18px;
    }",
    ),

  # Input: Selector for choosing model ----
      selectInput(inputId = "Model",
                  label = "Choose L-APIM (more information in panel About the Method):",
                  choices = c("Model 1: For distinguishable partners"=1,
                              "Model 2: For indistinguishable partners"=2, 
                              "Model 3: Group differences in actor and partner effects for distinguishable partners"=3, 
                              "Model 4: Group differences in actor and partner effects for indistinguishable partners"=4,
                              "Model 5: Continuous time-varying moderator C for distinguishable partners"=5, 
                              "Model 6: Continuous time-varying moderator C for indistinguishable partners"=6,  
                              "Model 7: Dichotomous time-varying moderator D for distinguishable partners"=7, 
                              "Model 8: Dichotomous time-varying moderator D for indistinguishable partners"=8, 
                              "Model 9: Quadratic actor and partner effects for distinguishable partners"=9, 
                              "Model 10: Quadratic actor and partner effects for indistinguishable partners"=10,
                              "Model 11: Group differences in quadratic actor and partner effects for distinguishable partners"=11,
                              "Model 12: Group differences in quadratic actor and partner effects for indistinguishable partners"=12,
                              "Model 13: Continuous time-varying moderator C with quadratic effects for distinguishable partners"=13,
                              "Model 14: Continuous time-varying moderator C with quadratic effects for indistinguishable partners"=14,
                              "Model 15: Dichotomous time-varying moderator D with quadratic effects for distinguishable partners"=15,
                              "Model 16: Dichotomous time-varying moderator D with quadratic effects for indistinguishable partners"=16,
                              "Model 17: AR(1) L-APIM. For distinguishable partners"=17,
                              "Model 18: AR(1) L-APIM. For indistinguishable partners"=18, 
                              "Model 19: AR(1) L-APIM. Group differences in actor and partner effects for distinguishable partners"=19, 
                              "Model 20: AR(1) L-APIM. Group differences in actor and partner effects for indistinguishable partners"=20,
                              "Model 21: AR(1) L-APIM. Continuous time-varying moderator C for distinguishable partners"=21, 
                              "Model 22: AR(1) L-APIM. Continuous time-varying moderator C for indistinguishable partners"=22,  
                              "Model 23: AR(1) L-APIM. Dichotomous time-varying moderator D for distinguishable partners"=23, 
                              "Model 24: AR(1) L-APIM. Dichotomous time-varying moderator D for indistinguishable partners"=24, 
                              "Model 25: AR(1) L-APIM. Quadratic actor and partner effects for distinguishable partners"=25, 
                              "Model 26: AR(1) L-APIM. Quadratic actor and partner effects for indistinguishable partners"=26,
                              "Model 27: AR(1) L-APIM. Group differences in quadratic actor and partner effects for distinguishable partners"=27,
                              "Model 28: AR(1) L-APIM. Group differences in quadratic actor and partner effects for indistinguishable partners"=28,
                              "Model 29: AR(1) L-APIM. Continuous time-varying moderator C with quadratic effects for distinguishable partners"=29,
                              "Model 30: AR(1) L-APIM. Continuous time-varying moderator C with quadratic effects for indistinguishable partners"=30,
                              "Model 31: AR(1) L-APIM. Dichotomous time-varying moderator D with quadratic effects for distinguishable partners"=31,
                              "Model 32: AR(1) L-APIM. Dichotomous time-varying moderator D with quadratic effects for indistinguishable partners"=32)),

# Add formulas

conditionalPanel(
condition = "input.Model == '3' || input.Model == '4' || input.Model == '11' || input.Model == '12' || 
         input.Model == '19' || input.Model == '20' || input.Model == '27' || input.Model == '28'", 
         helpText("Number of dyads: there are two possibilities either introduce an increasing sequence of comma-separated positive integers (e.g., 60, 70, 80, 90), where the length of the sequence must be the same in the two groups. The second option allows computing power for only one value for the number of dyads for each group."),
         textInput("N0.dyad","Number of dyads in Group 0 (reference group)", NULL),
         textInput("N1.dyad","Number of dyads in Group 1", NULL)
),

conditionalPanel(
condition = "input.Model == '1' || input.Model == '2' || input.Model == '5' || input.Model == '6' || input.Model == '7' || input.Model == '8' || 
input.Model == '9' || input.Model == '10' || input.Model == '13' || input.Model == '14' || input.Model == '15' || input.Model == '16' || 
input.Model == '17' || input.Model == '18' || input.Model == '21' || input.Model == '22' || input.Model == '23' || input.Model == '24' || 
input.Model == '25' || input.Model == '26' || input.Model == '29' || input.Model == '30' || input.Model == '31' || input.Model == '32'", 
         helpText("Number of dyads: there are two possibilities either introduce an increasing sequence of comma-separated positive integers (e.g., 60, 70, 80, 90). The second option allows computing power for only one value for the number of dyads."),
         textInput("N.dyad","Number of dyads", NULL)
),

    numericInput("T.obs", "Number of time points", NULL), 

conditionalPanel(
condition = "input.Model == '1' || input.Model == '5' || input.Model == '7' || input.Model == '9' || input.Model == '13' || input.Model == '15' || 
input.Model == '17' || input.Model == '21' || input.Model == '23' || input.Model == '25' || input.Model == '29' || input.Model == '31'", 
    numericInput("c.F", "Fixed intercept for partner A: \\( c_{A} \\)", NULL),
    numericInput("c.M", "Fixed intercept for partner B: \\( c_{B} \\)", NULL),
    numericInput("a.FF", "Fixed actor effect for partner A: \\( a_{AA} \\)", NULL),
    numericInput("p.MF", "Fixed partner effect for partner A: \\( p_{BA} \\)", NULL),
    numericInput("a.MM", "Fixed actor effect for partner B: \\( a_{BB} \\)", NULL),
    numericInput("p.FM", "Fixed partner effect for partner B: \\( p_{AB} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '2' || input.Model == '6' || input.Model == '8' || input.Model == '10' || input.Model == '14' || 
input.Model == '16' || input.Model == '18' || input.Model == '22' || input.Model == '24' || input.Model == '26' || 
input.Model == '30' || input.Model == '32'", 
    numericInput("c", "Fixed intercept: \\( c \\)", NULL),
    numericInput("a", "Fixed actor effect: \\( a \\)", NULL),
    numericInput("p", "Fixed partner effect: \\( p \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '3' || input.Model == '11' || input.Model == '19' || input.Model == '27'", 
    numericInput("c.F0", "Fixed intercept for partner A in Group 0: \\( c_{A0} \\)", NULL),
    numericInput("c.F1", "Difference in the fixed intercept between Group 0 and 1 for partner A: \\( c_{A1} \\)", NULL),
    numericInput("c.M0", "Fixed intercept for partner B in Group 0: \\( c_{B0} \\)", NULL),
    numericInput("c.M1", "Difference in the fixed intercept between Group 0 and 1 for partner B: \\( c_{B1} \\)", NULL),
    numericInput("a.FF0", "Fixed actor effect for partner A in Group 0: \\( a_{AA0} \\)", NULL),
    numericInput("a.FF1", "Difference in the fixed actor effect between Group 0 and 1 for partner A: \\( a_{AA1} \\)", NULL),
    numericInput("p.MF0", "Fixed partner effect for partner A in Group 0: \\( p_{BA0} \\)", NULL),
    numericInput("p.MF1", "Difference in the fixed partner effect between Group 0 and 1 for partner A: \\( p_{BA1} \\)", NULL),
    numericInput("a.MM0", "Fixed actor effect for partner B in Group 0: \\( a_{BB0} \\)", NULL),
    numericInput("a.MM1", "Difference in the fixed actor effect between Group 0 and 1 for partner B: \\( a_{BB1} \\)", NULL),
    numericInput("p.FM0", "Fixed partner effect for partner B in Group 0: \\( p_{AM0} \\)", NULL),
    numericInput("p.FM1", "Difference in the fixed partner effect between Group 0 and 1 for partner B: \\( p_{AB1} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '4' || input.Model == '12' || input.Model == '20' || input.Model == '28'", 
    numericInput("c0", "Fixed intercept in Group 0: \\( c_{0} \\)", NULL),
    numericInput("c1", "Difference in the fixed intercept between Group 0 and 1: \\( c_{1} \\)", NULL),
    numericInput("a0", "Fixed actor effect in Group 0: \\( a_{0} \\)", NULL),
    numericInput("a1", "Difference in the fixed actor effect between Group 0 and 1: \\( a_{1} \\)", NULL),
    numericInput("p0", "Fixed partner effect in Group 0: \\( p_{0} \\)", NULL),
    numericInput("p1", "Difference in the fixed partner effect between Group 0 and 1: \\( p_{1} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '5' || input.Model == '13' || input.Model == '21' || input.Model == '29'", 
    numericInput("b.F", "Effect of the time-varying continuous variable C for partner A: \\( b_{A} \\)", NULL),
    numericInput("b.FF", "Moderation effect of the time-varying continuous variable C in the actor effect for partner A: \\( b_{AA} \\)", NULL),
    numericInput("b.MF", "Moderation effect of the time-varying continuous variable C in the partner effect for partner A: \\( b_{BA} \\)", NULL),
    numericInput("b.M", "Effect of the time-varying continuous variable C for partner B: \\( b_{B} \\)", NULL),
    numericInput("b.MM", "Moderation effect of the time-varying continuous variable C in the actor effect for partner B: \\( b_{BB} \\)", NULL),
    numericInput("b.FM", "Moderation effect of the time-varying continuous variable C in the partner effect for partner B: \\( b_{AB} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '7' || input.Model == '15' || input.Model == '23' || input.Model == '31'", 
    numericInput("d.F", "Effect of the time-varying dichotomous variable D for partner A: \\( d_{A} \\)", NULL),
    numericInput("d.FF", "Moderation effect of the time-varying dichotomous variable D in the actor effect for partner A: \\( d_{AA} \\)", NULL),
    numericInput("d.MF", "Moderation effect of the time-varying dichotomous variable D in the partner effect for partner A: \\( d_{BA} \\)", NULL),
    numericInput("d.M", "Effect of the time-varying dichotomous variable D for partner B: \\( d_{B} \\)", NULL),
    numericInput("d.MM", "Moderation effect of the time-varying dichotomous variable D in the actor effect for partner B: \\( d_{BB} \\)", NULL),
    numericInput("d.FM", "Moderation effect of the time-varying dichotomous variable D in the partner effect for partner B: \\( d_{AB} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '6' || input.Model == '14' || input.Model == '22' || input.Model == '30'", 
    numericInput("b", "Effect of the time-varying continuous variable: \\( b \\)", NULL),
    numericInput("b.a", "Moderation effect of the time-varying continuous variable C in the actor effect: \\( b_{a} \\)", NULL),
    numericInput("b.p", "Moderation effect of the time-varying continuous variable C in the partner effect: \\( b_{p} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '8' || input.Model == '16' || input.Model == '24' || input.Model == '32'", 
    numericInput("d", "Effect of the time-varying continuous variable: \\( d \\)", NULL),
    numericInput("d.a", "Moderation effect of the time-varying dichotomous variable D in the actor effect: \\( d_{a} \\)", NULL),
    numericInput("d.p", "Moderation effect of the time-varying dichotomous variable D in the partner effect: \\( d_{p} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '9' || input.Model == '13' || input.Model == '15' || input.Model == '25' || input.Model == '29' || input.Model == '31'", 
    numericInput("a.FF2", "Fixed quadratic actor effect for partner A: \\( a_{AA2} \\)", NULL),
    numericInput("p.MF2", "Fixed quadratic partner effect for partner A: \\( p_{BA2} \\)", NULL),
    numericInput("a.MM2", "Fixed quadratic actor effect for partner B: \\( a_{BB2} \\)", NULL),
    numericInput("p.FM2", "Fixed quadratic partner effect for partner B: \\( p_{AB2} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '10' || input.Model == '14' || input.Model == '16' || input.Model == '26' || input.Model == '30' || input.Model == '32'", 
    numericInput("a.2", "Fixed quadratic actor effect: \\( a_2 \\)", NULL),
    numericInput("p.2", "Fixed quadratic partner effect: \\( p_2 \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '11' || input.Model == '27'", 
    numericInput("a.FF02", "Fixed quadratic actor effect for partner A in Group 0: \\( a_{AA02} \\)", NULL),
    numericInput("a.FF12", "Difference in the fixed quadratic actor effect between Group 0 and 1 for partner A: \\( a_{AA12} \\)", NULL),
    numericInput("p.MF02", "Fixed quadratic partner effect for partner A in Group 0: \\( p_{BA02} \\)", NULL),
    numericInput("p.MF12", "Difference in the fixed quadratic partner effect between Group 0 and 1 for partner A: \\( p_{BA12} \\)", NULL),
    numericInput("a.MM02", "Fixed quadratic actor effect for partner B in Group 0: \\( a_{BB02} \\)", NULL),
    numericInput("a.MM12", "Difference in the fixed quadratic actor effect between Group 0 and 1 for partner B: \\( a_{BB12} \\)", NULL),
    numericInput("p.FM02", "Fixed quadratic partner effect for partner B in Group 0: \\( p_{AB02} \\)", NULL),
    numericInput("p.FM12", "Difference in the fixed quadratic partner effect between Group 0 and 1 for partner B: \\( p_{AB12} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '12' || input.Model == '28'", 
    numericInput("a02", "Fixed quadratic actor effect in Group 0: \\( a_{02} \\)", NULL),
    numericInput("a12", "Difference in the fixed quadratic actor effect between Group 0 and 1: \\( a_{12} \\)", NULL),
    numericInput("p02", "Fixed quadratic partner effect in Group 0: \\( p_{02} \\)", NULL),
    numericInput("p12", "Difference in the fixed quadratic partner effect between Group 0 and 1: \\( p_{12} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '13' || input.Model == '29'", 
    numericInput("b.FF2", "Moderation effect of the time-varying continuous variable C in the quadratic actor effect for partner A: \\( b_{AA2} \\)", NULL),
    numericInput("b.MF2", "Moderation effect of the time-varying continuous variable C in the quadratic partner effect for partner A: \\( b_{BA2} \\)", NULL),
    numericInput("b.MM2", "Moderation effect of the time-varying continuous variable C in the quadratic actor effect for partner B: \\( b_{BB2} \\)", NULL),
    numericInput("b.FM2", "Moderation effect of the time-varying continuous variable C in the quadratic partner effect for partner B: \\( b_{AB2} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '15' || input.Model == '31'", 
    numericInput("d.FF2", "Moderation effect of the time-varying dichotomous variable D in the quadratic actor effect for partner A: \\( d_{AA2} \\)", NULL),
    numericInput("d.MF2", "Moderation effect of the time-varying dichotomous variable D in the quadratic partner effect for partner A: \\( d_{BA2} \\)", NULL),
    numericInput("d.MM2", "Moderation effect of the time-varying dichotomous variable D in the quadratic actor effect for partner B: \\( d_{BB2} \\)", NULL),
    numericInput("d.FM2", "Moderation effect of the time-varying dichotomous variable D in the quadratic partner effect for partner B: \\( d_{AB2} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '14' || input.Model == '30'", 
    numericInput("b.a2", "Moderation effect of the time-varying continuous variable C in the quadratic actor effect: \\( b_{a2} \\)", NULL),
    numericInput("b.p2", "Moderation effect of the time-varying continuous variable C in the quadratic partner effect: \\( b_{p2} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '16' || input.Model == '32'", 
    numericInput("d.a2", "Moderation effect of the time-varying dichotomous variable D in the quadratic actor effect: \\( d_{a2} \\)", NULL),
    numericInput("d.p2", "Moderation effect of the time-varying dichotomous variable D in the quadratic partner effect: \\( d_{p2} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '17' || input.Model == '21' || input.Model == '23' || input.Model == '25' || input.Model == '29' || input.Model == '31'", 
    numericInput("rho.YF", "Autoregressive effect for partner A: \\( \\rho_{YA} \\)", NULL),
    numericInput("rho.YM", "Autoregressive effect for partner B: \\( \\rho_{YB} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '18' || input.Model == '22' || input.Model == '24' || input.Model == '26' || input.Model == '30' || input.Model == '32'", 
    numericInput("rho.Y", "Autoregressive effect: \\( \\rho_{Y} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '19' || input.Model == '27'", 
    numericInput("rho.YF0", "Autoregressive effect for partner A in Group 0: \\( \\rho_{YA0} \\)", NULL),
    numericInput("rho.YF1", "Difference in the autoregressive effect for partner A between Group 0 and 1: \\( \\rho_{YA1} \\)", NULL),
    numericInput("rho.YM0", "Autoregressive effect for partner B in Group 0: \\( \\rho_{YB0} \\)", NULL),
    numericInput("rho.YM1", "Difference in the autoregressive effect for partner B between Group 0 and 1: \\( \\rho_{YB1} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '20' || input.Model == '28'", 
    numericInput("rho.Y0", "Autoregressive effect in Group 0: \\( \\rho_{Y0} \\)", NULL),
    numericInput("rho.Y1", "Difference in the autoregressive effect between Group 0 and 1: \\( \\rho_{Y1} \\)", NULL)
),

    numericInput("sigma.eps.F", "Standard deviation of Level 1 errors for partner A: \\( \\sigma_{\\epsilon_A} \\)", NULL), 
    numericInput("sigma.eps.M", "Standard deviation of Level 1 errors for partner B: \\( \\sigma_{\\epsilon_B} \\)", NULL), 
    numericInput("rho.eps.FM", "Correlation between the Level 1 errors for partners A and B: \\( \\rho_{\\epsilon_{AB}} \\)", NULL), 

conditionalPanel(
condition = "input.Model == '1' || input.Model == '3' || input.Model == '5' || input.Model == '7' || input.Model == '9' || 
input.Model == '11' || input.Model == '13' || input.Model == '15' || input.Model == '17' || input.Model == '19' || 
input.Model == '21' || input.Model == '23' || input.Model == '25' || input.Model == '27' || input.Model == '29' || input.Model == '31'",
    numericInput("sigma.nu.F", "Standard deviation of the random intercept for partner A: \\( \\sigma_{\\nu_A} \\)", NULL),
    numericInput("sigma.nu.M", "Standard deviation of the random intercept for partner B: \\( \\sigma_{\\nu_B} \\)", NULL), 
    numericInput("rho.nu.F.M", "Correlation between the random intercepts for partners A and B: \\( \\rho_{\\nu_{AB}} \\)", NULL)
),

conditionalPanel(
condition = "input.Model == '2' || input.Model == '4' || input.Model == '6' || input.Model == '8' || input.Model == '10' || 
input.Model == '12' || input.Model == '14' || input.Model == '16' || input.Model == '18' || input.Model == '20' || 
input.Model == '22' || input.Model == '24' || input.Model == '26' || input.Model == '28' || input.Model == '30' || input.Model == '32'",
    numericInput("sigma.nu", "Standard deviation of the random intercept: \\( \\sigma_{\\nu} \\)", NULL)),

conditionalPanel(
condition = "input.Model == '1' || input.Model == '2' || input.Model == '5' || input.Model == '6' || input.Model == '7' || 
input.Model == '8' || input.Model == '9' || input.Model == '10' || input.Model == '13' || input.Model == '14' || 
input.Model == '15' || input.Model == '16' || input.Model == '17' || input.Model == '18' || input.Model == '21' || 
input.Model == '22' || input.Model == '23' || input.Model == '24' || input.Model == '25' || input.Model == '26' || 
input.Model == '29' || input.Model == '30' || input.Model == '31' || input.Model == '32'", 
    numericInput("mu.XF", "Mean of time-varying variable X for partner A:", NULL), 
    numericInput("sigma.XF", "Standard deviation of time-varying variable X for partner A:", NULL),
    numericInput("mu.XM", "Mean of time-varying variable X for partner B:", NULL), 
    numericInput("sigma.XM", "Standard deviation of time-varying variable X for partner B:", NULL),
    numericInput("rho.X", "Correlation between time-varying variables X for partners A and B:", NULL) 
),

conditionalPanel(
condition = "input.Model == '3' || input.Model == '4' || input.Model == '11' || input.Model == '12' || 
input.Model == '19' || input.Model == '20' || input.Model == '27' || input.Model == '28'", 
    numericInput("mu.XF0", "Mean of time-varying variable X in Group 0 for partner A:", NULL), 
    numericInput("sigma.XF0", "Standard deviation of time-varying variable X in Group 0 for partner A:", NULL),
    numericInput("mu.XM0", "Mean of time-varying variable X in Group 0 for partner B:", NULL), 
    numericInput("sigma.XM0", "Standard deviation of time-varying variable X in Group 0 for partner B:", NULL),
    numericInput("rho.X0", "Correlation between time-varying variables X for partners A and B in Group 0:", NULL), 
    numericInput("mu.XF1", "Mean of time-varying variable X in Group 1 for partner A:", NULL), 
    numericInput("sigma.XF1", "Standard deviation of time-varying variable X in Group 1 for partner A:", NULL),
    numericInput("mu.XM1", "Mean of time-varying variable X in Group 1 for partner B:", NULL), 
    numericInput("sigma.XM1", "Standard deviation of time-varying variable X in Group 1 for partner B:", NULL),
    numericInput("rho.X1", "Correlation between time-varying variables X for partners A and B in Group 1:", NULL)  
),

checkboxInput(inputId = "is.center.X", label = strong("Person-mean center the time-varying predictor X"), value = TRUE),

conditionalPanel(
condition = "input.Model == '5' || input.Model == '6' || input.Model == '13' || input.Model == '14' || 
input.Model == '21' || input.Model == '22' || input.Model == '29' || input.Model == '30'", 
    numericInput("mu.W", "Mean of the time-varying continuous variable C:", NULL), 
    numericInput("sigma.W", "Standard deviation of the time-varying continuous variable C:", NULL), 
    checkboxInput(inputId = "is.center.W", label = strong("Center the time-varying continuous variable C"), value = TRUE)
),

conditionalPanel(
condition = "input.Model == '7' || input.Model == '8' || input.Model == '15' || input.Model == '16' || 
input.Model == '23' || input.Model == '24' || input.Model == '31' || input.Model == '32'", 
    numericInput("prob.D", "Probability that the time-varying dichotomous variable D is 1:", NULL)
),

selectInput(inputId = "is.REML",
                  label = "Choose the method to fit linear mixed-effects model",
                  choices = c("Maximizing the log-likelihood"=FALSE, 
                              "Maximizing the restricted log-likelihood"=TRUE)),

    numericInput("alpha", "Type I error: \\( \\alpha \\)", 0.05), 
    numericInput("R", "Monte Carlo Replicates", 1000),
    helpText("We recommend using at least 1000 Monte Carlo replicates"),

    actionButton(inputId = "input_action", label = "Compute Power"),

    actionButton("reset_button", "Reset Page"),

    helpText("Note:",
             "To switch models and set new parameters click the Reset Page button."),

    helpText("Contact: ginette.lafit@kuleuven.be"),
    helpText("Lafit, G., Sels, L., Adolf, J., Loeys, T., & Ceulemans, E. (2021, June 29). PowerLAPIM: An Application to Conduct Power Analysis for Longitudinal Actor-Partner Interdependence Models that Include Quadratic Effects. https://doi.org/10.31234/osf.io/mnce4")  

),


  mainPanel(
    tabsetPanel(
        
      tabPanel("Power Analysis",
      tags$h4("Simulation Progress"),        
      verbatimTextOutput("text"),
      tags$h4("Monte Carlo replicates that converged"),        
      verbatimTextOutput("replicates"),
      tags$h4("Power Analysis"),
      tags$h5("The power curve is shown when the vector of the number of dyads includes at least two elements."), 
      plotOutput("powerplot", height = "2800px")), 

      tabPanel("Summary Fixed Effects", 
      tags$h4("Summary Fixed Effects"),
      tags$h5("Note:"), 
      tags$h5("Mean is the average of the estimated parameter over the Monte Carlo replicates"),
      tags$h5("Std.error is the standard error of the estimated parameter over the Monte Carlo replicates"),
      tags$h5("Bias is the average of the difference between the estimated parameter and true parameter over the Monte Carlo replicates"),
      tags$h5("Confidence interval width is the average width of the (1-alpha)% confidence intervals over the Monte Carlo replicates"),
      tags$h5("(1-alpha)% Coverage is the average of the (1-alpha)% confidence intervals that include the true parameter over the Monte Carlo replicates"),
      tags$h5("Empirical power is the number of times the null hypothesis is rejected over the Monte Carlo replicates"),
      tags$h5("Std.error is the standard error of the empirical power over the Monte Carlo replicates"),

      formattableOutput("power", width ='1200px')
      ),

      tabPanel("Summary Variance Components", 
      tags$h4("Summary Variance Components"),
      tags$h5("Note:"), 
      tags$h5("Mean is the average of the estimated parameter over the Monte Carlo replicates"),
      tags$h5("Std.error is the standard error of the estimated parameter over the Monte Carlo replicates"),
      tags$h5("Bias is the average of the difference between the estimated parameter and true parameter over the Monte Carlo replicates"),
      formattableOutput("covariance", width = '1100px')
      ),
    
      tabPanel("About the Method",

      tags$h4("Simulation Approach to Estimate Power in Multilevel Linear Models"),
      tags$p("Algorithm:"),
      tags$p("1. Given a model based on the hypothesized theory, set up the population parameters. The parameter values can be decided from previous studies or a pilot study."),
      tags$p("2. Set the sample size and generate a data set based on the model and its population parameters."),
      tags$p("3. Test the significance of the null hypothesis (i.e. the hypothesized effect is zero) using the generated data with a Wald test."),
      tags$p("4. Repeat steps 2 and 3 for R times, where R is the number of Monte Carlo replications."),
      tags$p("5. Compute the power which is the number of times the hypothesized effect is significant over the Monte Carlo replications."),

      tags$h4("Tutorial"),
      tags$p(a("Lafit et al. (2021)",href="https://psyarxiv.com/mnce4/"),"includes a hands-on tutorial for conducting simulation-based power for the L-APIM."),

      tags$h4("Population models of interest"),
      tags$p("The", a("supplementary material",href="https://osf.io/vtb9e/"),"includes a description of the population models included in the application."))

    ), width = 10)


))
# Load packages
library(nlme)
library(MASS)
library(tidyverse)
library(future.apply)
library(gridExtra)
library(formattable)
library(htmltools)
library(shiny)
library(DT)
library(ggplot2)
library(gridExtra)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(shinyjs)

server <- shinyServer(function(input, output){

  power.sim = eventReactive(input$input_action, {

     withProgress(message = 'Running... (check estimated computational time)',
                   detail = 'Go get a cup of tea...',style = getShinyOption("progress.style",
     default = "old"),{
                     
     withCallingHandlers({
     shinyjs::html("text", "")

    Sim.model.Dyad.IL(input$Model,input$N.dyad,input$N0.dyad,input$N1.dyad,input$T.obs,  
    input$c.F,input$c.M,input$a.FF,input$p.MF,input$a.FF2,input$p.MF2,input$a.MM,input$p.FM,input$a.MM2,input$p.FM2,
    input$c,input$a,input$a.2,input$p,input$p.2,
    input$c.F0,input$c.F1,input$c.M0,input$c.M1,input$a.FF0,input$a.FF1,input$a.FF02,input$a.FF12,input$p.MF0,input$p.MF1,input$p.MF02,input$p.MF12,
    input$a.MM0,input$a.MM1,input$a.MM02,input$a.MM12,input$p.FM0,input$p.FM1,input$p.FM02,input$p.FM12,
    input$c0,input$c1,input$a0,input$a1,input$a02,input$a12,input$p0,input$p1,input$p02,input$p12,
    input$b.F,input$b.M,input$b.FF,input$b.MF,input$b.MM,input$b.FM,input$b.FF2,input$b.MF2,input$b.MM2,input$b.FM2,
    input$d.F,input$d.M,input$d.FF,input$d.MF,input$d.MM,input$d.FM,input$d.FF2,input$d.MF2,input$d.MM2,input$d.FM2,
    input$b,input$b.a,input$b.a2,input$b.p,input$b.p2,
    input$d,input$d.a,input$d.a2,input$d.p,input$d.p2,
    input$rho.YF,input$rho.YM,input$rho.Y,input$rho.YF0,input$rho.YF1,input$rho.YM0,input$rho.YM1,
    input$sigma.eps.F,input$sigma.eps.M,input$rho.eps.FM,
    input$sigma.nu.F,input$sigma.nu.M,input$rho.nu.F.M,input$sigma.nu,
    input$mu.XF,input$sigma.XF,input$mu.XM,input$sigma.XM,input$rho.X,
    input$mu.XF0,input$mu.XF1,input$sigma.XF0,input$sigma.XF1,input$mu.XM0,input$mu.XM1,input$sigma.XM0,input$sigma.XM1,input$rho.X0,input$rho.X1,
    input$mu.W,input$sigma.W,input$prob.D,
    input$is.center.X,input$is.center.W,input$R,input$alpha,input$is.REML) 
      },
        message = function(m) {
          shinyjs::html(id = "text", html = m$message, add = TRUE)
      })
    })
  })
  
  observeEvent(input$reset_button, {
      shinyjs::reset("side-panel")
    })

   output$powerplot <-renderPlot({
   
     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==1){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.a.FF,Power.p.MF,Power.a.MM,Power.p.FM,ncol=1,heights=c(4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==2){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.a,Power.p,ncol=1,heights=c(4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==3){
     power.IL.c.F0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F0.se')])
     data.IL.c.F0 = data.frame(Participants=c(paste(power.IL.c.F0$N.Dyad.Group.0.,power.IL.c.F0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F0$c.F0,se=power.IL.c.F0$c.F0.se)
     Power.c.F0 = ggplot(data = data.IL.c.F0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A in Group 0", x = "Number of dyads")

     power.IL.c.F1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F1','c.F1.se')])
     data.IL.c.F1 = data.frame(Participants=c(paste(power.IL.c.F1$N.Dyad.Group.0.,power.IL.c.F1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F1$c.F1,se=power.IL.c.F1$c.F1.se)
     Power.c.F1 = ggplot(data = data.IL.c.F1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.c.M0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M0','c.M0.se')])
     data.IL.c.M0 = data.frame(Participants=c(paste(power.IL.c.M0$N.Dyad.Group.0.,power.IL.c.M0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M0$c.M0,se=power.IL.c.M0$c.M0.se)
     Power.c.M0 = ggplot(data = data.IL.c.M0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B in Group 0", x = "Number of dyads")

     power.IL.c.M1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M1','c.M1.se')])
     data.IL.c.M1 = data.frame(Participants=c(paste(power.IL.c.M1$N.Dyad.Group.0.,power.IL.c.M1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M1$c.M1,se=power.IL.c.M1$c.M1.se)
     Power.c.M1 = ggplot(data = data.IL.c.M1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.a.FF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF0','a.FF0.se')])
     data.IL.a.FF0 = data.frame(Participants=c(paste(power.IL.a.FF0$N.Dyad.Group.0.,power.IL.a.FF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF0$a.FF0,se=power.IL.a.FF0$a.FF0.se)
     Power.a.FF0 = ggplot(data = data.IL.a.FF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A in Group 0", x = "Number of dyads")

     power.IL.a.FF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF1','a.FF1.se')])
     data.IL.a.FF1 = data.frame(Participants=c(paste(power.IL.a.FF1$N.Dyad.Group.0.,power.IL.a.FF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF1$a.FF1,se=power.IL.a.FF1$a.FF1.se)
     Power.a.FF1 = ggplot(data = data.IL.a.FF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed actor effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.p.MF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF0','p.MF0.se')])
     data.IL.p.MF0 = data.frame(Participants=c(paste(power.IL.p.MF0$N.Dyad.Group.0.,power.IL.p.MF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF0$p.MF0,se=power.IL.p.MF0$p.MF0.se)
     Power.p.MF0 = ggplot(data = data.IL.p.MF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A in Group 0", x = "Number of dyads")

     power.IL.p.MF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF1','p.MF1.se')])
     data.IL.p.MF1 = data.frame(Participants=c(paste(power.IL.p.MF1$N.Dyad.Group.0.,power.IL.p.MF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF1$p.MF1,se=power.IL.p.MF1$p.MF1.se)
     Power.p.MF1 = ggplot(data = data.IL.p.MF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed partner effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.a.MM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM0','a.MM0.se')])
     data.IL.a.MM0 = data.frame(Participants=c(paste(power.IL.a.MM0$N.Dyad.Group.0.,power.IL.a.MM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM0$a.MM0,se=power.IL.a.MM0$a.MM0.se)
     Power.a.MM0 = ggplot(data = data.IL.a.MM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B in Group 0", x = "Number of dyads")

     power.IL.a.MM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM1','a.MM1.se')])
     data.IL.a.MM1 = data.frame(Participants=c(paste(power.IL.a.MM1$N.Dyad.Group.0.,power.IL.a.MM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM1$a.MM1,se=power.IL.a.MM1$a.MM1.se)
     Power.a.MM1 = ggplot(data = data.IL.a.MM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed actor effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.p.FM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM0','p.FM0.se')])
     data.IL.p.FM0 = data.frame(Participants=c(paste(power.IL.p.FM0$N.Dyad.Group.0.,power.IL.p.FM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM0$p.FM0,se=power.IL.p.FM0$p.FM0.se)
     Power.p.FM0 = ggplot(data = data.IL.p.FM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B in Group 0", x = "Number of dyads")

     power.IL.p.FM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM1','p.FM1.se')])
     data.IL.p.FM1 = data.frame(Participants=c(paste(power.IL.p.FM1$N.Dyad.Group.0.,power.IL.p.FM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM1$p.FM1,se=power.IL.p.FM1$p.FM1.se)
     Power.p.FM1 = ggplot(data = data.IL.p.FM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed partner effect of partner B between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F0,Power.c.F1,Power.c.M0,Power.c.M1,
     Power.a.FF0,Power.a.FF1,Power.p.MF0,Power.p.MF1,
     Power.a.MM0,Power.a.MM1,Power.p.FM0,Power.p.MF1,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==4){
     power.IL.c0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c0.se')])
     data.IL.c0 = data.frame(Participants=c(paste(power.IL.c0$N.Dyad.Group.0.,power.IL.c0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c0$c0,se=power.IL.c0$c0.se)
     Power.c0 = ggplot(data = data.IL.c0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept in Group 0", x = "Number of dyads")

     power.IL.c1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c1','c1.se')])
     data.IL.c1 = data.frame(Participants=c(paste(power.IL.c1$N.Dyad.Group.0.,power.IL.c1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c1$c1,se=power.IL.c1$c1.se)
     Power.c1 = ggplot(data = data.IL.c1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept between Group 0 and 1", x = "Number of dyads")

     power.IL.a0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a0','a0.se')])
     data.IL.a0 = data.frame(Participants=c(paste(power.IL.a0$N.Dyad.Group.0.,power.IL.a0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a0$a0,se=power.IL.a0$a0.se)
     Power.a0 = ggplot(data = data.IL.a0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect F in Group 0", x = "Number of dyads")

     power.IL.a1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a1','a1.se')])
     data.IL.a1 = data.frame(Participants=c(paste(power.IL.a1$N.Dyad.Group.0.,power.IL.a1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a1$a1,se=power.IL.a1$a1.se)
     Power.a1 = ggplot(data = data.IL.a1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed actor effect between Group 0 and 1", x = "Number of dyads")

     power.IL.p0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p0','p0.se')])
     data.IL.p0 = data.frame(Participants=c(paste(power.IL.p0$N.Dyad.Group.0.,power.IL.p0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p0$p0,se=power.IL.p0$p0.se)
     Power.p0 = ggplot(data = data.IL.p0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect in Group 0", x = "Number of dyads")

     power.IL.p1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p1','p1.se')])
     data.IL.p1 = data.frame(Participants=c(paste(power.IL.p1$N.Dyad.Group.0.,power.IL.p1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p1$p1,se=power.IL.p1$p1.se)
     Power.p1 = ggplot(data = data.IL.p1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed partner effect between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c0,Power.c1,Power.a0,Power.a1,Power.p0,Power.p1,ncol=1,heights=c(4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==5){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A", x = "Number of dyads")

     power.IL.b.F = data.frame(power.sim()$power[,c('N dyad','b.F','b.F.se')])
     data.IL.b.F = data.frame(Participants=c(paste(power.IL.b.F$N.dyad,sep=";")),
     Power=power.IL.b.F$b.F,se=power.IL.b.F$b.F.se)
     Power.b.F = ggplot(data = data.IL.b.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner A", x = "Number of dyads")

     power.IL.b.FF = data.frame(power.sim()$power[,c('N dyad','b.FF','b.FF.se')])
     data.IL.b.FF = data.frame(Participants=c(paste(power.IL.b.FF$N.dyad,sep=";")),
     Power=power.IL.b.FF$b.FF,se=power.IL.b.FF$b.FF.se)
     Power.b.FF = ggplot(data = data.IL.b.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the actor effect of partner A", x = "Number of dyads")

     power.IL.b.MF = data.frame(power.sim()$power[,c('N dyad','b.MF','b.MF.se')])
     data.IL.b.MF = data.frame(Participants=c(paste(power.IL.b.MF$N.dyad,sep=";")),
     Power=power.IL.b.MF$b.MF,se=power.IL.b.MF$b.MF.se)
     Power.b.MF = ggplot(data = data.IL.b.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B", x = "Number of dyads")

     power.IL.b.M = data.frame(power.sim()$power[,c('N dyad','b.M','b.M.se')])
     data.IL.b.M = data.frame(Participants=c(paste(power.IL.b.M$N.dyad,sep=";")),
     Power=power.IL.b.M$b.M,se=power.IL.b.M$b.M.se)
     Power.b.M = ggplot(data = data.IL.b.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner B", x = "Number of dyads")

     power.IL.b.MM = data.frame(power.sim()$power[,c('N dyad','b.MM','b.MM.se')])
     data.IL.b.MM = data.frame(Participants=c(paste(power.IL.b.MM$N.dyad,sep=";")),
     Power=power.IL.b.MM$b.MM,se=power.IL.b.MM$b.MM.se)
     Power.b.MM = ggplot(data = data.IL.b.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the actor effect of partner B", x = "Number of dyads")

     power.IL.b.FM = data.frame(power.sim()$power[,c('N dyad','b.FM','b.FM.se')])
     data.IL.b.FM = data.frame(Participants=c(paste(power.IL.b.FM$N.dyad,sep=";")),
     Power=power.IL.b.FM$b.FM,se=power.IL.b.FM$b.FM.se)
     Power.b.FM = ggplot(data = data.IL.b.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.a.FF,Power.p.MF,
     Power.b.F,Power.b.FF,Power.b.MF,
     Power.a.MM,Power.p.FM,Power.b.M,Power.b.MM,Power.b.FM,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==6){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect", x = "Number of dyads")

     power.IL.b = data.frame(power.sim()$power[,c('N dyad','b','b.se')])
     data.IL.b = data.frame(Participants=c(paste(power.IL.b$N.dyad,sep=";")),
     Power=power.IL.b$b,se=power.IL.b$b.se)
     Power.b = ggplot(data = data.IL.b, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C", x = "Number of dyads")

     power.IL.b.a = data.frame(power.sim()$power[,c('N dyad','b.a','b.a.se')])
     data.IL.b.a = data.frame(Participants=c(paste(power.IL.b.a$N.dyad,sep=";")),
     Power=power.IL.b.a$b.a,se=power.IL.b.a$b.a.se)
     Power.b.a = ggplot(data = data.IL.b.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the actor effect", x = "Number of dyads")

     power.IL.b.p = data.frame(power.sim()$power[,c('N dyad','b.p','b.p.se')])
     data.IL.b.p = data.frame(Participants=c(paste(power.IL.b.p$N.dyad,sep=";")),
     Power=power.IL.b.p$b.p,se=power.IL.b.p$b.p.se)
     Power.b.p = ggplot(data = data.IL.b.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.a,Power.p,
     Power.b,Power.b.a,Power.b.p,ncol=1,heights=c(4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==7){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A", x = "Number of dyads")

     power.IL.d.F = data.frame(power.sim()$power[,c('N dyad','d.F','d.F.se')])
     data.IL.d.F = data.frame(Participants=c(paste(power.IL.d.F$N.dyad,sep=";")),
     Power=power.IL.d.F$d.F,se=power.IL.d.F$d.F.se)
     Power.d.F = ggplot(data = data.IL.d.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner A", x = "Number of dyads")

     power.IL.d.FF = data.frame(power.sim()$power[,c('N dyad','d.FF','d.FF.se')])
     data.IL.d.FF = data.frame(Participants=c(paste(power.IL.d.FF$N.dyad,sep=";")),
     Power=power.IL.d.FF$d.FF,se=power.IL.d.FF$d.FF.se)
     Power.d.FF = ggplot(data = data.IL.d.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the actor effect of partner A", x = "Number of dyads")

     power.IL.d.MF = data.frame(power.sim()$power[,c('N dyad','d.MF','d.MF.se')])
     data.IL.d.MF = data.frame(Participants=c(paste(power.IL.d.MF$N.dyad,sep=";")),
     Power=power.IL.d.MF$d.MF,se=power.IL.d.MF$d.MF.se)
     Power.d.MF = ggplot(data = data.IL.d.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B", x = "Number of dyads")

     power.IL.d.M = data.frame(power.sim()$power[,c('N dyad','d.M','d.M.se')])
     data.IL.d.M = data.frame(Participants=c(paste(power.IL.d.M$N.dyad,sep=";")),
     Power=power.IL.d.M$d.M,se=power.IL.d.M$d.M.se)
     Power.d.M = ggplot(data = data.IL.d.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner B", x = "Number of dyads")

     power.IL.d.MM = data.frame(power.sim()$power[,c('N dyad','d.MM','d.MM.se')])
     data.IL.d.MM = data.frame(Participants=c(paste(power.IL.d.MM$N.dyad,sep=";")),
     Power=power.IL.d.MM$d.MM,se=power.IL.d.MM$d.MM.se)
     Power.d.MM = ggplot(data = data.IL.d.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the actor effect of partner B", x = "Number of dyads")

     power.IL.d.FM = data.frame(power.sim()$power[,c('N dyad','d.FM','d.FM.se')])
     data.IL.d.FM = data.frame(Participants=c(paste(power.IL.d.FM$N.dyad,sep=";")),
     Power=power.IL.d.FM$d.FM,se=power.IL.d.FM$d.FM.se)
     Power.d.FM = ggplot(data = data.IL.d.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.a.FF,Power.p.MF,
     Power.d.F,Power.d.FF,Power.d.MF,
     Power.a.MM,Power.p.FM,Power.d.M,Power.d.MM,Power.d.FM,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==8){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect", x = "Number of dyads")

     power.IL.d = data.frame(power.sim()$power[,c('N dyad','d','d.se')])
     data.IL.d = data.frame(Participants=c(paste(power.IL.d$N.dyad,sep=";")),
     Power=power.IL.d$d,se=power.IL.d$d.se)
     Power.d = ggplot(data = data.IL.d, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D", x = "Number of dyads")

     power.IL.d.a = data.frame(power.sim()$power[,c('N dyad','d.a','d.a.se')])
     data.IL.d.a = data.frame(Participants=c(paste(power.IL.d.a$N.dyad,sep=";")),
     Power=power.IL.d.a$d.a,se=power.IL.d.a$d.a.se)
     Power.d.a = ggplot(data = data.IL.d.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the actor effect", x = "Number of dyads")

     power.IL.d.p = data.frame(power.sim()$power[,c('N dyad','d.p','d.p.se')])
     data.IL.d.p = data.frame(Participants=c(paste(power.IL.d.p$N.dyad,sep=";")),
     Power=power.IL.d.p$d.p,se=power.IL.d.p$d.p.se)
     Power.d.p = ggplot(data = data.IL.d.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.a,Power.p,
     Power.d,Power.d.a,Power.d.p,ncol=1,heights=c(4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==9){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A", x = "Number of dyads")

     power.IL.a.FF2 = data.frame(power.sim()$power[,c('N dyad','a.FF2','a.FF2.se')])
     data.IL.a.FF2 = data.frame(Participants=c(paste(power.IL.a.FF2$N.dyad,sep=";")),
     Power=power.IL.a.FF2$a.FF2,se=power.IL.a.FF2$a.FF2.se)
     Power.a.FF2 = ggplot(data = data.IL.a.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A", x = "Number of dyads")

     power.IL.p.MF2 = data.frame(power.sim()$power[,c('N dyad','p.MF2','p.MF2.se')])
     data.IL.p.MF2 = data.frame(Participants=c(paste(power.IL.p.MF2$N.dyad,sep=";")),
     Power=power.IL.p.MF2$p.MF2,se=power.IL.p.MF2$p.MF2.se)
     Power.p.MF2 = ggplot(data = data.IL.p.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B", x = "Number of dyads")

     power.IL.a.MM2 = data.frame(power.sim()$power[,c('N dyad','a.MM2','a.MM2.se')])
     data.IL.a.MM2 = data.frame(Participants=c(paste(power.IL.a.MM2$N.dyad,sep=";")),
     Power=power.IL.a.MM2$a.MM2,se=power.IL.a.MM2$a.MM2.se)
     Power.a.MM2 = ggplot(data = data.IL.a.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B", x = "Number of dyads")

     power.IL.p.FM2 = data.frame(power.sim()$power[,c('N dyad','p.FM2','p.FM2.se')])
     data.IL.p.FM2 = data.frame(Participants=c(paste(power.IL.p.FM2$N.dyad,sep=";")),
     Power=power.IL.p.FM2$p.FM2,se=power.IL.p.FM2$p.FM2.se)
     Power.p.FM2 = ggplot(data = data.IL.p.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.a.FF,Power.a.FF2,
     Power.p.MF,Power.p.MF2,Power.a.MM,Power.a.MM2,Power.p.FM,Power.p.FM2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==10){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect", x = "Number of dyads")

     power.IL.a.2 = data.frame(power.sim()$power[,c('N dyad','a.2','a.2.se')])
     data.IL.a.2 = data.frame(Participants=c(paste(power.IL.a.2$N.dyad,sep=";")),
     Power=power.IL.a.2$a.2,se=power.IL.a.2$a.2.se)
     Power.a.2 = ggplot(data = data.IL.a.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect", x = "Number of dyads")

     power.IL.p.2 = data.frame(power.sim()$power[,c('N dyad','p.2','p.2.se')])
     data.IL.p.2 = data.frame(Participants=c(paste(power.IL.p.2$N.dyad,sep=";")),
     Power=power.IL.p.2$p.2,se=power.IL.p.2$p.2.se)
     Power.p.2 = ggplot(data = data.IL.p.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.a,Power.a.2,Power.p,Power.p.2,ncol=1,heights=c(4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==11){
     power.IL.c.F0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F0.se')])
     data.IL.c.F0 = data.frame(Participants=c(paste(power.IL.c.F0$N.Dyad.Group.0.,power.IL.c.F0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F0$c.F0,se=power.IL.c.F0$c.F0.se)
     Power.c.F0 = ggplot(data = data.IL.c.F0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A in Group 0", x = "Number of dyads")

     power.IL.c.F1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F1','c.F1.se')])
     data.IL.c.F1 = data.frame(Participants=c(paste(power.IL.c.F1$N.Dyad.Group.0.,power.IL.c.F1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F1$c.F1,se=power.IL.c.F1$c.F1.se)
     Power.c.F1 = ggplot(data = data.IL.c.F1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.c.M0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M0','c.M0.se')])
     data.IL.c.M0 = data.frame(Participants=c(paste(power.IL.c.M0$N.Dyad.Group.0.,power.IL.c.M0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M0$c.M0,se=power.IL.c.M0$c.M0.se)
     Power.c.M0 = ggplot(data = data.IL.c.M0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed intercept of partner B in Group 0", x = "Number of dyads")

     power.IL.c.M1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M1','c.M1.se')])
     data.IL.c.M1 = data.frame(Participants=c(paste(power.IL.c.M1$N.Dyad.Group.0.,power.IL.c.M1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M1$c.M1,se=power.IL.c.M1$c.M1.se)
     Power.c.M1 = ggplot(data = data.IL.c.M1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the linear fixed intercept of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.a.FF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF0','a.FF0.se')])
     data.IL.a.FF0 = data.frame(Participants=c(paste(power.IL.a.FF0$N.Dyad.Group.0.,power.IL.a.FF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF0$a.FF0,se=power.IL.a.FF0$a.FF0.se)
     Power.a.FF0 = ggplot(data = data.IL.a.FF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A in Group 0", x = "Number of dyads")

     power.IL.a.FF02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF02','a.FF02.se')])
     data.IL.a.FF02 = data.frame(Participants=c(paste(power.IL.a.FF02$N.Dyad.Group.0.,power.IL.a.FF02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF02$a.FF02,se=power.IL.a.FF02$a.FF02.se)
     Power.a.FF02 = ggplot(data = data.IL.a.FF02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A in Group 0", x = "Number of dyads")

     power.IL.a.FF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF1','a.FF1.se')])
     data.IL.a.FF1 = data.frame(Participants=c(paste(power.IL.a.FF1$N.Dyad.Group.0.,power.IL.a.FF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF1$a.FF1,se=power.IL.a.FF1$a.FF1.se)
     Power.a.FF1 = ggplot(data = data.IL.a.FF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed actor effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.a.FF12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF12','a.FF12.se')])
     data.IL.a.FF12 = data.frame(Participants=c(paste(power.IL.a.FF12$N.Dyad.Group.0.,power.IL.a.FF12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF12$a.FF12,se=power.IL.a.FF12$a.FF12.se)
     Power.a.FF12 = ggplot(data = data.IL.a.FF12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed actor effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.p.MF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF0','p.MF0.se')])
     data.IL.p.MF0 = data.frame(Participants=c(paste(power.IL.p.MF0$N.Dyad.Group.0.,power.IL.p.MF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF0$p.MF0,se=power.IL.p.MF0$p.MF0.se)
     Power.p.MF0 = ggplot(data = data.IL.p.MF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A in Group 0", x = "Number of dyads")

     power.IL.p.MF02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF02','p.MF02.se')])
     data.IL.p.MF02 = data.frame(Participants=c(paste(power.IL.p.MF02$N.Dyad.Group.0.,power.IL.p.MF02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF02$p.MF02,se=power.IL.p.MF02$p.MF02.se)
     Power.p.MF02 = ggplot(data = data.IL.p.MF02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A in Group 0", x = "Number of dyads")

     power.IL.p.MF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF1','p.MF1.se')])
     data.IL.p.MF1 = data.frame(Participants=c(paste(power.IL.p.MF1$N.Dyad.Group.0.,power.IL.p.MF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF1$p.MF1,se=power.IL.p.MF1$p.MF1.se)
     Power.p.MF1 = ggplot(data = data.IL.p.MF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.p.MF12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF12','p.MF12.se')])
     data.IL.p.MF12 = data.frame(Participants=c(paste(power.IL.p.MF12$N.Dyad.Group.0.,power.IL.p.MF12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF12$p.MF12,se=power.IL.p.MF12$p.MF12.se)
     Power.p.MF12 = ggplot(data = data.IL.p.MF12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed partner effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.a.MM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM0','a.MM0.se')])
     data.IL.a.MM0 = data.frame(Participants=c(paste(power.IL.a.MM0$N.Dyad.Group.0.,power.IL.a.MM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM0$a.MM0,se=power.IL.a.MM0$a.MM0.se)
     Power.a.MM0 = ggplot(data = data.IL.a.MM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B in Group 0", x = "Number of dyads")

     power.IL.a.MM02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM02','a.MM02.se')])
     data.IL.a.MM02 = data.frame(Participants=c(paste(power.IL.a.MM02$N.Dyad.Group.0.,power.IL.a.MM02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM02$a.MM02,se=power.IL.a.MM02$a.MM02.se)
     Power.a.MM02 = ggplot(data = data.IL.a.MM02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B in Group 0", x = "Number of dyads")

     power.IL.a.MM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM1','a.MM1.se')])
     data.IL.a.MM1 = data.frame(Participants=c(paste(power.IL.a.MM1$N.Dyad.Group.0.,power.IL.a.MM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM1$a.MM1,se=power.IL.a.MM1$a.MM1.se)
     Power.a.MM1 = ggplot(data = data.IL.a.MM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed actor effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.a.MM12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM12','a.MM12.se')])
     data.IL.a.MM12 = data.frame(Participants=c(paste(power.IL.a.MM12$N.Dyad.Group.0.,power.IL.a.MM12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM12$a.MM12,se=power.IL.a.MM12$a.MM12.se)
     Power.a.MM12 = ggplot(data = data.IL.a.MM12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed actor effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.p.FM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM0','p.FM0.se')])
     data.IL.p.FM0 = data.frame(Participants=c(paste(power.IL.p.FM0$N.Dyad.Group.0.,power.IL.p.FM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM0$p.FM0,se=power.IL.p.FM0$p.FM0.se)
     Power.p.FM0 = ggplot(data = data.IL.p.FM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B in Group 0", x = "Number of dyads")

     power.IL.p.FM02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM02','p.FM02.se')])
     data.IL.p.FM02 = data.frame(Participants=c(paste(power.IL.p.FM02$N.Dyad.Group.0.,power.IL.p.FM02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM02$p.FM02,se=power.IL.p.FM02$p.FM02.se)
     Power.p.FM02 = ggplot(data = data.IL.p.FM02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B in Group 0", x = "Number of dyads")

     power.IL.p.FM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM1','p.FM1.se')])
     data.IL.p.FM1 = data.frame(Participants=c(paste(power.IL.p.FM1$N.Dyad.Group.0.,power.IL.p.FM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM1$p.FM1,se=power.IL.p.FM1$p.FM1.se)
     Power.p.FM1 = ggplot(data = data.IL.p.FM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.p.FM12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM12','p.FM12.se')])
     data.IL.p.FM12 = data.frame(Participants=c(paste(power.IL.p.FM12$N.Dyad.Group.0.,power.IL.p.FM12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM12$p.FM12,se=power.IL.p.FM12$p.FM12.se)
     Power.p.FM12 = ggplot(data = data.IL.p.FM12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed partner effect of partner B between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F0,Power.c.F1,Power.c.M0,Power.c.M1,
     Power.a.FF0,Power.a.FF02,Power.a.FF1,Power.a.FF12,
     Power.p.MF0,Power.p.MF02,Power.p.MF1,Power.p.MF12,
     Power.a.MM0,Power.a.MM02,Power.a.MM1,Power.a.MM12,
     Power.p.FM0,Power.p.FM02,Power.p.MF1,Power.p.MF12,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==12){
     power.IL.c0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c0.se')])
     data.IL.c0 = data.frame(Participants=c(paste(power.IL.c0$N.Dyad.Group.0.,power.IL.c0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c0$c0,se=power.IL.c0$c0.se)
     Power.c0 = ggplot(data = data.IL.c0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept in Group 0", x = "Number of dyads")

     power.IL.c1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c1','c1.se')])
     data.IL.c1 = data.frame(Participants=c(paste(power.IL.c1$N.Dyad.Group.0.,power.IL.c1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c1$c1,se=power.IL.c1$c1.se)
     Power.c1 = ggplot(data = data.IL.c1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept between Group 0 and 1", x = "Number of dyads")

     power.IL.a0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a0','a0.se')])
     data.IL.a0 = data.frame(Participants=c(paste(power.IL.a0$N.Dyad.Group.0.,power.IL.a0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a0$a0,se=power.IL.a0$a0.se)
     Power.a0 = ggplot(data = data.IL.a0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect F in Group 0", x = "Number of dyads")

     power.IL.a02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a02','a02.se')])
     data.IL.a02 = data.frame(Participants=c(paste(power.IL.a02$N.Dyad.Group.0.,power.IL.a02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a02$a02,se=power.IL.a02$a02.se)
     Power.a02 = ggplot(data = data.IL.a02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect F in Group 0", x = "Number of dyads")

     power.IL.a1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a1','a1.se')])
     data.IL.a1 = data.frame(Participants=c(paste(power.IL.a1$N.Dyad.Group.0.,power.IL.a1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a1$a1,se=power.IL.a1$a1.se)
     Power.a1 = ggplot(data = data.IL.a1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed actor effect between Group 0 and 1", x = "Number of dyads")

     power.IL.a12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a12','a12.se')])
     data.IL.a12 = data.frame(Participants=c(paste(power.IL.a12$N.Dyad.Group.0.,power.IL.a12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a12$a12,se=power.IL.a12$a12.se)
     Power.a12 = ggplot(data = data.IL.a12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed actor effect between Group 0 and 1", x = "Number of dyads")

     power.IL.p0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p0','p0.se')])
     data.IL.p0 = data.frame(Participants=c(paste(power.IL.p0$N.Dyad.Group.0.,power.IL.p0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p0$p0,se=power.IL.p0$p0.se)
     Power.p0 = ggplot(data = data.IL.p0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect in Group 0", x = "Number of dyads")

     power.IL.p02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p02','p02.se')])
     data.IL.p02 = data.frame(Participants=c(paste(power.IL.p02$N.Dyad.Group.0.,power.IL.p02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p02$p02,se=power.IL.p02$p02.se)
     Power.p02 = ggplot(data = data.IL.p02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect in Group 0", x = "Number of dyads")

     power.IL.p1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p1','p1.se')])
     data.IL.p1 = data.frame(Participants=c(paste(power.IL.p1$N.Dyad.Group.0.,power.IL.p1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p1$p1,se=power.IL.p1$p1.se)
     Power.p1 = ggplot(data = data.IL.p1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect between Group 0 and 1", x = "Number of dyads")

     power.IL.p12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p12','p12.se')])
     data.IL.p12 = data.frame(Participants=c(paste(power.IL.p12$N.Dyad.Group.0.,power.IL.p12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p12$p12,se=power.IL.p12$p12.se)
     Power.p12 = ggplot(data = data.IL.p12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c0,Power.c1,Power.a0,Power.a02,Power.a1,Power.a12,
     Power.p0,Power.p02,Power.p1,Power.p12,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==13){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A", x = "Number of dyads")

     power.IL.a.FF2 = data.frame(power.sim()$power[,c('N dyad','a.FF2','a.FF2.se')])
     data.IL.a.FF2 = data.frame(Participants=c(paste(power.IL.a.FF2$N.dyad,sep=";")),
     Power=power.IL.a.FF2$a.FF2,se=power.IL.a.FF2$a.FF2.se)
     Power.a.FF2 = ggplot(data = data.IL.a.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A", x = "Number of dyads")

     power.IL.p.MF2 = data.frame(power.sim()$power[,c('N dyad','p.MF2','p.MF2.se')])
     data.IL.p.MF2 = data.frame(Participants=c(paste(power.IL.p.MF2$N.dyad,sep=";")),
     Power=power.IL.p.MF2$p.MF2,se=power.IL.p.MF2$p.MF2.se)
     Power.p.MF2 = ggplot(data = data.IL.p.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A", x = "Number of dyads")

     power.IL.b.F = data.frame(power.sim()$power[,c('N dyad','b.F','b.F.se')])
     data.IL.b.F = data.frame(Participants=c(paste(power.IL.b.F$N.dyad,sep=";")),
     Power=power.IL.b.F$b.F,se=power.IL.b.F$b.F.se)
     Power.b.F = ggplot(data = data.IL.b.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner A", x = "Number of dyads")

     power.IL.b.FF = data.frame(power.sim()$power[,c('N dyad','b.FF','b.FF.se')])
     data.IL.b.FF = data.frame(Participants=c(paste(power.IL.b.FF$N.dyad,sep=";")),
     Power=power.IL.b.FF$b.FF,se=power.IL.b.FF$b.FF.se)
     Power.b.FF = ggplot(data = data.IL.b.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the actor effect of partner A", x = "Number of dyads")

     power.IL.b.FF2 = data.frame(power.sim()$power[,c('N dyad','b.FF2','b.FF2.se')])
     data.IL.b.FF2 = data.frame(Participants=c(paste(power.IL.b.FF2$N.dyad,sep=";")),
     Power=power.IL.b.FF2$b.FF2,se=power.IL.b.FF2$b.FF2.se)
     Power.b.FF2 = ggplot(data = data.IL.b.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the actor effect of partner A", x = "Number of dyads")

     power.IL.b.MF = data.frame(power.sim()$power[,c('N dyad','b.MF','b.MF.se')])
     data.IL.b.MF = data.frame(Participants=c(paste(power.IL.b.MF$N.dyad,sep=";")),
     Power=power.IL.b.MF$b.MF,se=power.IL.b.MF$b.MF.se)
     Power.b.MF = ggplot(data = data.IL.b.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the partner effect of partner A", x = "Number of dyads")

     power.IL.b.MF2 = data.frame(power.sim()$power[,c('N dyad','b.MF2','b.MF2.se')])
     data.IL.b.MF2 = data.frame(Participants=c(paste(power.IL.b.MF2$N.dyad,sep=";")),
     Power=power.IL.b.MF2$b.MF2,se=power.IL.b.MF2$b.MF2.se)
     Power.b.MF2 = ggplot(data = data.IL.b.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B", x = "Number of dyads")

     power.IL.a.MM2 = data.frame(power.sim()$power[,c('N dyad','a.MM2','a.MM2.se')])
     data.IL.a.MM2 = data.frame(Participants=c(paste(power.IL.a.MM2$N.dyad,sep=";")),
     Power=power.IL.a.MM2$a.MM2,se=power.IL.a.MM2$a.MM2.se)
     Power.a.MM2 = ggplot(data = data.IL.a.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B", x = "Number of dyads")

     power.IL.p.FM2 = data.frame(power.sim()$power[,c('N dyad','p.FM2','p.FM2.se')])
     data.IL.p.FM2 = data.frame(Participants=c(paste(power.IL.p.FM2$N.dyad,sep=";")),
     Power=power.IL.p.FM2$p.FM2,se=power.IL.p.FM2$p.FM2.se)
     Power.p.FM2 = ggplot(data = data.IL.p.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B", x = "Number of dyads")

     power.IL.b.M = data.frame(power.sim()$power[,c('N dyad','b.M','b.M.se')])
     data.IL.b.M = data.frame(Participants=c(paste(power.IL.b.M$N.dyad,sep=";")),
     Power=power.IL.b.M$b.M,se=power.IL.b.M$b.M.se)
     Power.b.M = ggplot(data = data.IL.b.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner B", x = "Number of dyads")

     power.IL.b.MM = data.frame(power.sim()$power[,c('N dyad','b.MM','b.MM.se')])
     data.IL.b.MM = data.frame(Participants=c(paste(power.IL.b.MM$N.dyad,sep=";")),
     Power=power.IL.b.MM$b.MM,se=power.IL.b.MM$b.MM.se)
     Power.b.MM = ggplot(data = data.IL.b.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the actor effect of partner B", x = "Number of dyads")

     power.IL.b.MM2 = data.frame(power.sim()$power[,c('N dyad','b.MM2','b.MM2.se')])
     data.IL.b.MM2 = data.frame(Participants=c(paste(power.IL.b.MM2$N.dyad,sep=";")),
     Power=power.IL.b.MM2$b.MM2,se=power.IL.b.MM2$b.MM2.se)
     Power.b.MM2 = ggplot(data = data.IL.b.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the actor effect of partner B", x = "Number of dyads")

     power.IL.b.FM = data.frame(power.sim()$power[,c('N dyad','b.FM','b.FM.se')])
     data.IL.b.FM = data.frame(Participants=c(paste(power.IL.b.FM$N.dyad,sep=";")),
     Power=power.IL.b.FM$b.FM,se=power.IL.b.FM$b.FM.se)
     Power.b.FM = ggplot(data = data.IL.b.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the partner effect of partner B", x = "Number of dyads")

     power.IL.b.FM2 = data.frame(power.sim()$power[,c('N dyad','b.FM2','b.FM2.se')])
     data.IL.b.FM2 = data.frame(Participants=c(paste(power.IL.b.FM2$N.dyad,sep=";")),
     Power=power.IL.b.FM2$b.FM2,se=power.IL.b.FM2$b.FM2.se)
     Power.b.FM2 = ggplot(data = data.IL.b.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.a.FF,Power.a.FF2,
     Power.p.MF,Power.p.MF2,Power.b.F,
     Power.b.FF,Power.b.FF2,Power.b.MF,Power.b.MF2,
     Power.a.MM,Power.a.MM2,Power.p.FM,Power.p.FM2,Power.b.M,
     Power.b.MM,Power.b.MM2,Power.b.FM,Power.b.FM2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==14){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect", x = "Number of dyads")

     power.IL.a.2 = data.frame(power.sim()$power[,c('N dyad','a.2','a.2.se')])
     data.IL.a.2 = data.frame(Participants=c(paste(power.IL.a.2$N.dyad,sep=";")),
     Power=power.IL.a.2$a.2,se=power.IL.a.2$a.2.se)
     Power.a.2 = ggplot(data = data.IL.a.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect", x = "Number of dyads")

     power.IL.p.2 = data.frame(power.sim()$power[,c('N dyad','p.2','p.2.se')])
     data.IL.p.2 = data.frame(Participants=c(paste(power.IL.p.2$N.dyad,sep=";")),
     Power=power.IL.p.2$p.2,se=power.IL.p.2$p.2.se)
     Power.p.2 = ggplot(data = data.IL.p.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect", x = "Number of dyads")

     power.IL.b = data.frame(power.sim()$power[,c('N dyad','b','b.se')])
     data.IL.b = data.frame(Participants=c(paste(power.IL.b$N.dyad,sep=";")),
     Power=power.IL.b$b,se=power.IL.b$b.se)
     Power.b = ggplot(data = data.IL.b, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C", x = "Number of dyads")

     power.IL.b.a = data.frame(power.sim()$power[,c('N dyad','b.a','b.a.se')])
     data.IL.b.a = data.frame(Participants=c(paste(power.IL.b.a$N.dyad,sep=";")),
     Power=power.IL.b.a$b.a,se=power.IL.b.a$b.a.se)
     Power.b.a = ggplot(data = data.IL.b.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the actor effect", x = "Number of dyads")

     power.IL.b.a2 = data.frame(power.sim()$power[,c('N dyad','b.a2','b.a2.se')])
     data.IL.b.a2 = data.frame(Participants=c(paste(power.IL.b.a2$N.dyad,sep=";")),
     Power=power.IL.b.a2$b.a2,se=power.IL.b.a2$b.a2.se)
     Power.b.a2 = ggplot(data = data.IL.b.a2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.a2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the actor effect", x = "Number of dyads")

     power.IL.b.p = data.frame(power.sim()$power[,c('N dyad','b.p','b.p.se')])
     data.IL.b.p = data.frame(Participants=c(paste(power.IL.b.p$N.dyad,sep=";")),
     Power=power.IL.b.p$b.p,se=power.IL.b.p$b.p.se)
     Power.b.p = ggplot(data = data.IL.b.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the partner effect", x = "Number of dyads")

     power.IL.b.p2 = data.frame(power.sim()$power[,c('N dyad','b.p2','b.p2.se')])
     data.IL.b.p2 = data.frame(Participants=c(paste(power.IL.b.p2$N.dyad,sep=";")),
     Power=power.IL.b.p2$b.p2,se=power.IL.b.p2$b.p2.se)
     Power.b.p2 = ggplot(data = data.IL.b.p2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.p2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.a,Power.a.2,Power.p,Power.p.2,
     Power.b,Power.b.a,Power.b.a2,Power.b.p,Power.b.p2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==15){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A", x = "Number of dyads")

     power.IL.a.FF2 = data.frame(power.sim()$power[,c('N dyad','a.FF2','a.FF2.se')])
     data.IL.a.FF2 = data.frame(Participants=c(paste(power.IL.a.FF2$N.dyad,sep=";")),
     Power=power.IL.a.FF2$a.FF2,se=power.IL.a.FF2$a.FF2.se)
     Power.a.FF2 = ggplot(data = data.IL.a.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A", x = "Number of dyads")

     power.IL.p.MF2 = data.frame(power.sim()$power[,c('N dyad','p.MF2','p.MF2.se')])
     data.IL.p.MF2 = data.frame(Participants=c(paste(power.IL.p.MF2$N.dyad,sep=";")),
     Power=power.IL.p.MF2$p.MF2,se=power.IL.p.MF2$p.MF2.se)
     Power.p.MF2 = ggplot(data = data.IL.p.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A", x = "Number of dyads")

     power.IL.d.F = data.frame(power.sim()$power[,c('N dyad','d.F','d.F.se')])
     data.IL.d.F = data.frame(Participants=c(paste(power.IL.d.F$N.dyad,sep=";")),
     Power=power.IL.d.F$d.F,se=power.IL.d.F$d.F.se)
     Power.d.F = ggplot(data = data.IL.d.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner A", x = "Number of dyads")

     power.IL.d.FF = data.frame(power.sim()$power[,c('N dyad','d.FF','d.FF.se')])
     data.IL.d.FF = data.frame(Participants=c(paste(power.IL.d.FF$N.dyad,sep=";")),
     Power=power.IL.d.FF$d.FF,se=power.IL.d.FF$d.FF.se)
     Power.d.FF = ggplot(data = data.IL.d.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the actor effect of partner A", x = "Number of dyads")

     power.IL.d.FF2 = data.frame(power.sim()$power[,c('N dyad','d.FF2','d.FF2.se')])
     data.IL.d.FF2 = data.frame(Participants=c(paste(power.IL.d.FF2$N.dyad,sep=";")),
     Power=power.IL.d.FF2$d.FF2,se=power.IL.d.FF2$d.FF2.se)
     Power.d.FF2 = ggplot(data = data.IL.d.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the actor effect of partner A", x = "Number of dyads")

     power.IL.d.MF = data.frame(power.sim()$power[,c('N dyad','d.MF','d.MF.se')])
     data.IL.d.MF = data.frame(Participants=c(paste(power.IL.d.MF$N.dyad,sep=";")),
     Power=power.IL.d.MF$d.MF,se=power.IL.d.MF$d.MF.se)
     Power.d.MF = ggplot(data = data.IL.d.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the partner effect of partner A", x = "Number of dyads")

     power.IL.d.MF2 = data.frame(power.sim()$power[,c('N dyad','d.MF2','d.MF2.se')])
     data.IL.d.MF2 = data.frame(Participants=c(paste(power.IL.d.MF2$N.dyad,sep=";")),
     Power=power.IL.d.MF2$d.MF2,se=power.IL.d.MF2$d.MF2.se)
     Power.d.MF2 = ggplot(data = data.IL.d.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B", x = "Number of dyads")

     power.IL.a.MM2 = data.frame(power.sim()$power[,c('N dyad','a.MM2','a.MM2.se')])
     data.IL.a.MM2 = data.frame(Participants=c(paste(power.IL.a.MM2$N.dyad,sep=";")),
     Power=power.IL.a.MM2$a.MM2,se=power.IL.a.MM2$a.MM2.se)
     Power.a.MM2 = ggplot(data = data.IL.a.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B", x = "Number of dyads")

     power.IL.p.FM2 = data.frame(power.sim()$power[,c('N dyad','p.FM2','p.FM2.se')])
     data.IL.p.FM2 = data.frame(Participants=c(paste(power.IL.p.FM2$N.dyad,sep=";")),
     Power=power.IL.p.FM2$p.FM2,se=power.IL.p.FM2$p.FM2.se)
     Power.p.FM2 = ggplot(data = data.IL.p.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B", x = "Number of dyads")

     power.IL.d.M = data.frame(power.sim()$power[,c('N dyad','d.M','d.M.se')])
     data.IL.d.M = data.frame(Participants=c(paste(power.IL.d.M$N.dyad,sep=";")),
     Power=power.IL.d.M$d.M,se=power.IL.d.M$d.M.se)
     Power.d.M = ggplot(data = data.IL.d.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner B", x = "Number of dyads")

     power.IL.d.MM = data.frame(power.sim()$power[,c('N dyad','d.MM','d.MM.se')])
     data.IL.d.MM = data.frame(Participants=c(paste(power.IL.d.MM$N.dyad,sep=";")),
     Power=power.IL.d.MM$d.MM,se=power.IL.d.MM$d.MM.se)
     Power.d.MM = ggplot(data = data.IL.d.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the actor effect of partner B", x = "Number of dyads")

     power.IL.d.MM2 = data.frame(power.sim()$power[,c('N dyad','d.MM2','d.MM2.se')])
     data.IL.d.MM2 = data.frame(Participants=c(paste(power.IL.d.MM2$N.dyad,sep=";")),
     Power=power.IL.d.MM2$d.MM2,se=power.IL.d.MM2$d.MM2.se)
     Power.d.MM2 = ggplot(data = data.IL.d.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the actor effect of partner B", x = "Number of dyads")

     power.IL.d.FM = data.frame(power.sim()$power[,c('N dyad','d.FM','d.FM.se')])
     data.IL.d.FM = data.frame(Participants=c(paste(power.IL.d.FM$N.dyad,sep=";")),
     Power=power.IL.d.FM$d.FM,se=power.IL.d.FM$d.FM.se)
     Power.d.FM = ggplot(data = data.IL.d.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the partner effect of partner B", x = "Number of dyads")

     power.IL.d.FM2 = data.frame(power.sim()$power[,c('N dyad','d.FM2','d.FM2.se')])
     data.IL.d.FM2 = data.frame(Participants=c(paste(power.IL.d.FM2$N.dyad,sep=";")),
     Power=power.IL.d.FM2$d.FM2,se=power.IL.d.FM2$d.FM2.se)
     Power.d.FM2 = ggplot(data = data.IL.d.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.a.FF,Power.a.FF2,
     Power.p.MF,Power.p.MF2,Power.d.F,
     Power.d.FF,Power.d.FF2,Power.d.MF,Power.d.MF2,
     Power.a.MM,Power.a.MM2,Power.p.FM,Power.p.FM2,
     Power.d.M,Power.d.MM,Power.d.MM2,Power.d.FM,Power.d.FM2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==16){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect", x = "Number of dyads")

     power.IL.a.2 = data.frame(power.sim()$power[,c('N dyad','a.2','a.2.se')])
     data.IL.a.2 = data.frame(Participants=c(paste(power.IL.a.2$N.dyad,sep=";")),
     Power=power.IL.a.2$a.2,se=power.IL.a.2$a.2.se)
     Power.a.2 = ggplot(data = data.IL.a.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect", x = "Number of dyads")

     power.IL.p.2 = data.frame(power.sim()$power[,c('N dyad','p.2','p.2.se')])
     data.IL.p.2 = data.frame(Participants=c(paste(power.IL.p.2$N.dyad,sep=";")),
     Power=power.IL.p.2$p.2,se=power.IL.p.2$p.2.se)
     Power.p.2 = ggplot(data = data.IL.p.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect", x = "Number of dyads")

     power.IL.d = data.frame(power.sim()$power[,c('N dyad','d','d.se')])
     data.IL.d = data.frame(Participants=c(paste(power.IL.d$N.dyad,sep=";")),
     Power=power.IL.d$d,se=power.IL.d$d.se)
     Power.d = ggplot(data = data.IL.d, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D", x = "Number of dyads")

     power.IL.d.a = data.frame(power.sim()$power[,c('N dyad','d.a','d.a.se')])
     data.IL.d.a = data.frame(Participants=c(paste(power.IL.d.a$N.dyad,sep=";")),
     Power=power.IL.d.a$d.a,se=power.IL.d.a$d.a.se)
     Power.d.a = ggplot(data = data.IL.d.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the actor effect", x = "Number of dyads")

     power.IL.d.a2 = data.frame(power.sim()$power[,c('N dyad','d.a2','d.a2.se')])
     data.IL.d.a2 = data.frame(Participants=c(paste(power.IL.d.a2$N.dyad,sep=";")),
     Power=power.IL.d.a2$d.a2,se=power.IL.d.a2$d.a2.se)
     Power.d.a2 = ggplot(data = data.IL.d.a2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.a2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the actor effect", x = "Number of dyads")

     power.IL.d.p = data.frame(power.sim()$power[,c('N dyad','d.p','d.p.se')])
     data.IL.d.p = data.frame(Participants=c(paste(power.IL.d.p$N.dyad,sep=";")),
     Power=power.IL.d.p$d.p,se=power.IL.d.p$d.p.se)
     Power.d.p = ggplot(data = data.IL.d.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the partner effect", x = "Number of dyads")

     power.IL.d.p2 = data.frame(power.sim()$power[,c('N dyad','d.p2','d.p2.se')])
     data.IL.d.p2 = data.frame(Participants=c(paste(power.IL.d.p2$N.dyad,sep=";")),
     Power=power.IL.d.p2$d.p2,se=power.IL.d.p2$d.p2.se)
     Power.d.p2 = ggplot(data = data.IL.d.p2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.p2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.a,Power.a.2,Power.p,Power.p.2,
     Power.d,Power.d.a,Power.d.a2,Power.d.p,Power.d.p2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==17){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.rho.YF = data.frame(power.sim()$power[,c('N dyad','rho.YF','rho.YF.se')])
     data.IL.rho.YF = data.frame(Participants=c(paste(power.IL.rho.YF$N.dyad,sep=";")),
     Power=power.IL.rho.YF$rho.YF,se=power.IL.rho.YF$rho.YF.se)
     Power.rho.YF = ggplot(data = data.IL.rho.YF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A", x = "Number of dyads")

     power.IL.rho.YM = data.frame(power.sim()$power[,c('N dyad','rho.YM','rho.YM.se')])
     data.IL.rho.YM = data.frame(Participants=c(paste(power.IL.rho.YM$N.dyad,sep=";")),
     Power=power.IL.rho.YM$rho.YM,se=power.IL.rho.YM$rho.YM.se)
     Power.rho.YM = ggplot(data = data.IL.rho.YM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.rho.YF,Power.rho.YM,Power.a.FF,Power.p.MF,Power.a.MM,Power.p.FM,ncol=1,heights=c(4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==18){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.rho.Y = data.frame(power.sim()$power[,c('N dyad','rho.Y','rho.Y.se')])
     data.IL.rho.Y = data.frame(Participants=c(paste(power.IL.rho.Y$N.dyad,sep=";")),
     Power=power.IL.rho.Y$rho.Y,se=power.IL.rho.Y$rho.Y.se)
     Power.rho.Y = ggplot(data = data.IL.rho.Y, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.rho.Y,Power.a,Power.p,ncol=1,heights=c(4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==19){
     power.IL.c.F0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F0.se')])
     data.IL.c.F0 = data.frame(Participants=c(paste(power.IL.c.F0$N.Dyad.Group.0.,power.IL.c.F0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F0$c.F0,se=power.IL.c.F0$c.F0.se)
     Power.c.F0 = ggplot(data = data.IL.c.F0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A in Group 0", x = "Number of dyads")

     power.IL.c.F1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F1','c.F1.se')])
     data.IL.c.F1 = data.frame(Participants=c(paste(power.IL.c.F1$N.Dyad.Group.0.,power.IL.c.F1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F1$c.F1,se=power.IL.c.F1$c.F1.se)
     Power.c.F1 = ggplot(data = data.IL.c.F1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.c.M0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M0','c.M0.se')])
     data.IL.c.M0 = data.frame(Participants=c(paste(power.IL.c.M0$N.Dyad.Group.0.,power.IL.c.M0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M0$c.M0,se=power.IL.c.M0$c.M0.se)
     Power.c.M0 = ggplot(data = data.IL.c.M0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B in Group 0", x = "Number of dyads")

     power.IL.c.M1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M1','c.M1.se')])
     data.IL.c.M1 = data.frame(Participants=c(paste(power.IL.c.M1$N.Dyad.Group.0.,power.IL.c.M1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M1$c.M1,se=power.IL.c.M1$c.M1.se)
     Power.c.M1 = ggplot(data = data.IL.c.M1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.rho.YF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YF0','rho.YF0.se')])
     data.IL.rho.YF0 = data.frame(Participants=c(paste(power.IL.rho.YF0$N.Dyad.Group.0.,power.IL.rho.YF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YF0$rho.YF0,se=power.IL.rho.YF0$rho.YF0.se)
     Power.rho.YF0 = ggplot(data = data.IL.rho.YF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A in Group 0", x = "Number of dyads")

     power.IL.rho.YF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YF1','rho.YF1.se')])
     data.IL.rho.YF1 = data.frame(Participants=c(paste(power.IL.rho.YF1$N.Dyad.Group.0.,power.IL.rho.YF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YF1$rho.YF1,se=power.IL.rho.YF1$rho.YF1.se)
     Power.rho.YF1 = ggplot(data = data.IL.rho.YF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed autoregressive effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.rho.YM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YM0','rho.YM0.se')])
     data.IL.rho.YM0 = data.frame(Participants=c(paste(power.IL.rho.YM0$N.Dyad.Group.0.,power.IL.rho.YM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YM0$rho.YM0,se=power.IL.rho.YM0$rho.YM0.se)
     Power.rho.YM0 = ggplot(data = data.IL.rho.YM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B in Group 0", x = "Number of dyads")

     power.IL.rho.YM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YM1','rho.YM1.se')])
     data.IL.rho.YM1 = data.frame(Participants=c(paste(power.IL.rho.YM1$N.Dyad.Group.0.,power.IL.rho.YM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YM1$rho.YM1,se=power.IL.rho.YM1$rho.YM1.se)
     Power.rho.YM1 = ggplot(data = data.IL.rho.YM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed autoregressive effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.a.FF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF0','a.FF0.se')])
     data.IL.a.FF0 = data.frame(Participants=c(paste(power.IL.a.FF0$N.Dyad.Group.0.,power.IL.a.FF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF0$a.FF0,se=power.IL.a.FF0$a.FF0.se)
     Power.a.FF0 = ggplot(data = data.IL.a.FF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A in Group 0", x = "Number of dyads")

     power.IL.a.FF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF1','a.FF1.se')])
     data.IL.a.FF1 = data.frame(Participants=c(paste(power.IL.a.FF1$N.Dyad.Group.0.,power.IL.a.FF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF1$a.FF1,se=power.IL.a.FF1$a.FF1.se)
     Power.a.FF1 = ggplot(data = data.IL.a.FF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed actor effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.p.MF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF0','p.MF0.se')])
     data.IL.p.MF0 = data.frame(Participants=c(paste(power.IL.p.MF0$N.Dyad.Group.0.,power.IL.p.MF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF0$p.MF0,se=power.IL.p.MF0$p.MF0.se)
     Power.p.MF0 = ggplot(data = data.IL.p.MF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A in Group 0", x = "Number of dyads")

     power.IL.p.MF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF1','p.MF1.se')])
     data.IL.p.MF1 = data.frame(Participants=c(paste(power.IL.p.MF1$N.Dyad.Group.0.,power.IL.p.MF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF1$p.MF1,se=power.IL.p.MF1$p.MF1.se)
     Power.p.MF1 = ggplot(data = data.IL.p.MF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed partner effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.a.MM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM0','a.MM0.se')])
     data.IL.a.MM0 = data.frame(Participants=c(paste(power.IL.a.MM0$N.Dyad.Group.0.,power.IL.a.MM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM0$a.MM0,se=power.IL.a.MM0$a.MM0.se)
     Power.a.MM0 = ggplot(data = data.IL.a.MM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B in Group 0", x = "Number of dyads")

     power.IL.a.MM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM1','a.MM1.se')])
     data.IL.a.MM1 = data.frame(Participants=c(paste(power.IL.a.MM1$N.Dyad.Group.0.,power.IL.a.MM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM1$a.MM1,se=power.IL.a.MM1$a.MM1.se)
     Power.a.MM1 = ggplot(data = data.IL.a.MM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed actor effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.p.FM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM0','p.FM0.se')])
     data.IL.p.FM0 = data.frame(Participants=c(paste(power.IL.p.FM0$N.Dyad.Group.0.,power.IL.p.FM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM0$p.FM0,se=power.IL.p.FM0$p.FM0.se)
     Power.p.FM0 = ggplot(data = data.IL.p.FM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B in Group 0", x = "Number of dyads")

     power.IL.p.FM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM1','p.FM1.se')])
     data.IL.p.FM1 = data.frame(Participants=c(paste(power.IL.p.FM1$N.Dyad.Group.0.,power.IL.p.FM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM1$p.FM1,se=power.IL.p.FM1$p.FM1.se)
     Power.p.FM1 = ggplot(data = data.IL.p.FM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed partner effect of partner B between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F0,Power.c.F1,Power.c.M0,Power.c.M1,
     Power.rho.YF0,Power.rho.YF1,Power.rho.YM0,Power.rho.YM1,
     Power.a.FF0,Power.a.FF1,Power.p.MF0,Power.p.MF1,
     Power.a.MM0,Power.a.MM1,Power.p.FM0,Power.p.MF1,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==20){
     power.IL.c0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c0.se')])
     data.IL.c0 = data.frame(Participants=c(paste(power.IL.c0$N.Dyad.Group.0.,power.IL.c0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c0$c0,se=power.IL.c0$c0.se)
     Power.c0 = ggplot(data = data.IL.c0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept in Group 0", x = "Number of dyads")

     power.IL.c1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c1','c1.se')])
     data.IL.c1 = data.frame(Participants=c(paste(power.IL.c1$N.Dyad.Group.0.,power.IL.c1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c1$c1,se=power.IL.c1$c1.se)
     Power.c1 = ggplot(data = data.IL.c1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept between Group 0 and 1", x = "Number of dyads")

     power.IL.rho.Y0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.Y0','rho.Y0.se')])
     data.IL.rho.Y0 = data.frame(Participants=c(paste(power.IL.rho.Y0$N.Dyad.Group.0.,power.IL.rho.YF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.Y0$rho.Y0,se=power.IL.rho.Y0$rho.Y0.se)
     Power.rho.Y0 = ggplot(data = data.IL.rho.Y0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect in Group 0", x = "Number of dyads")

     power.IL.rho.Y1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.Y1','rho.Y1.se')])
     data.IL.rho.Y1 = data.frame(Participants=c(paste(power.IL.rho.Y1$N.Dyad.Group.0.,power.IL.rho.Y1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.Y1$rho.Y1,se=power.IL.rho.Y1$rho.Y1.se)
     Power.rho.Y1 = ggplot(data = data.IL.rho.Y1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed autoregressive effect between Group 0 and 1", x = "Number of dyads")

     power.IL.a0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a0','a0.se')])
     data.IL.a0 = data.frame(Participants=c(paste(power.IL.a0$N.Dyad.Group.0.,power.IL.a0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a0$a0,se=power.IL.a0$a0.se)
     Power.a0 = ggplot(data = data.IL.a0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect F in Group 0", x = "Number of dyads")

     power.IL.a1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a1','a1.se')])
     data.IL.a1 = data.frame(Participants=c(paste(power.IL.a1$N.Dyad.Group.0.,power.IL.a1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a1$a1,se=power.IL.a1$a1.se)
     Power.a1 = ggplot(data = data.IL.a1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed actor effect between Group 0 and 1", x = "Number of dyads")

     power.IL.p0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p0','p0.se')])
     data.IL.p0 = data.frame(Participants=c(paste(power.IL.p0$N.Dyad.Group.0.,power.IL.p0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p0$p0,se=power.IL.p0$p0.se)
     Power.p0 = ggplot(data = data.IL.p0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect in Group 0", x = "Number of dyads")

     power.IL.p1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p1','p1.se')])
     data.IL.p1 = data.frame(Participants=c(paste(power.IL.p1$N.Dyad.Group.0.,power.IL.p1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p1$p1,se=power.IL.p1$p1.se)
     Power.p1 = ggplot(data = data.IL.p1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the fixed partner effect between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c0,Power.c1,Power.rho.Y0,Power.rho.Y1,Power.a0,Power.a1,Power.p0,Power.p1,ncol=1,heights=c(4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==21){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.rho.YF = data.frame(power.sim()$power[,c('N dyad','rho.YF','rho.YF.se')])
     data.IL.rho.YF = data.frame(Participants=c(paste(power.IL.rho.YF$N.dyad,sep=";")),
     Power=power.IL.rho.YF$rho.YF,se=power.IL.rho.YF$rho.YF.se)
     Power.rho.YF = ggplot(data = data.IL.rho.YF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A", x = "Number of dyads")

     power.IL.rho.YM = data.frame(power.sim()$power[,c('N dyad','rho.YM','rho.YM.se')])
     data.IL.rho.YM = data.frame(Participants=c(paste(power.IL.rho.YM$N.dyad,sep=";")),
     Power=power.IL.rho.YM$rho.YM,se=power.IL.rho.YM$rho.YM.se)
     Power.rho.YM = ggplot(data = data.IL.rho.YM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A", x = "Number of dyads")

     power.IL.b.F = data.frame(power.sim()$power[,c('N dyad','b.F','b.F.se')])
     data.IL.b.F = data.frame(Participants=c(paste(power.IL.b.F$N.dyad,sep=";")),
     Power=power.IL.b.F$b.F,se=power.IL.b.F$b.F.se)
     Power.b.F = ggplot(data = data.IL.b.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner A", x = "Number of dyads")

     power.IL.b.FF = data.frame(power.sim()$power[,c('N dyad','b.FF','b.FF.se')])
     data.IL.b.FF = data.frame(Participants=c(paste(power.IL.b.FF$N.dyad,sep=";")),
     Power=power.IL.b.FF$b.FF,se=power.IL.b.FF$b.FF.se)
     Power.b.FF = ggplot(data = data.IL.b.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the actor effect of partner A", x = "Number of dyads")

     power.IL.b.MF = data.frame(power.sim()$power[,c('N dyad','b.MF','b.MF.se')])
     data.IL.b.MF = data.frame(Participants=c(paste(power.IL.b.MF$N.dyad,sep=";")),
     Power=power.IL.b.MF$b.MF,se=power.IL.b.MF$b.MF.se)
     Power.b.MF = ggplot(data = data.IL.b.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B", x = "Number of dyads")

     power.IL.b.M = data.frame(power.sim()$power[,c('N dyad','b.M','b.M.se')])
     data.IL.b.M = data.frame(Participants=c(paste(power.IL.b.M$N.dyad,sep=";")),
     Power=power.IL.b.M$b.M,se=power.IL.b.M$b.M.se)
     Power.b.M = ggplot(data = data.IL.b.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner B", x = "Number of dyads")

     power.IL.b.MM = data.frame(power.sim()$power[,c('N dyad','b.MM','b.MM.se')])
     data.IL.b.MM = data.frame(Participants=c(paste(power.IL.b.MM$N.dyad,sep=";")),
     Power=power.IL.b.MM$b.MM,se=power.IL.b.MM$b.MM.se)
     Power.b.MM = ggplot(data = data.IL.b.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the actor effect of partner B", x = "Number of dyads")

     power.IL.b.FM = data.frame(power.sim()$power[,c('N dyad','b.FM','b.FM.se')])
     data.IL.b.FM = data.frame(Participants=c(paste(power.IL.b.FM$N.dyad,sep=";")),
     Power=power.IL.b.FM$b.FM,se=power.IL.b.FM$b.FM.se)
     Power.b.FM = ggplot(data = data.IL.b.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.rho.YF,Power.rho.YM,Power.a.FF,Power.p.MF,
     Power.b.F,Power.b.FF,Power.b.MF,
     Power.a.MM,Power.p.FM,Power.b.M,Power.b.MM,Power.b.FM,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==22){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.rho.Y = data.frame(power.sim()$power[,c('N dyad','rho.Y','rho.Y.se')])
     data.IL.rho.Y = data.frame(Participants=c(paste(power.IL.rho.Y$N.dyad,sep=";")),
     Power=power.IL.rho.Y$rho.Y,se=power.IL.rho.Y$rho.Y.se)
     Power.rho.Y = ggplot(data = data.IL.rho.Y, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect", x = "Number of dyads")

     power.IL.b = data.frame(power.sim()$power[,c('N dyad','b','b.se')])
     data.IL.b = data.frame(Participants=c(paste(power.IL.b$N.dyad,sep=";")),
     Power=power.IL.b$b,se=power.IL.b$b.se)
     Power.b = ggplot(data = data.IL.b, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C", x = "Number of dyads")

     power.IL.b.a = data.frame(power.sim()$power[,c('N dyad','b.a','b.a.se')])
     data.IL.b.a = data.frame(Participants=c(paste(power.IL.b.a$N.dyad,sep=";")),
     Power=power.IL.b.a$b.a,se=power.IL.b.a$b.a.se)
     Power.b.a = ggplot(data = data.IL.b.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the actor effect", x = "Number of dyads")

     power.IL.b.p = data.frame(power.sim()$power[,c('N dyad','b.p','b.p.se')])
     data.IL.b.p = data.frame(Participants=c(paste(power.IL.b.p$N.dyad,sep=";")),
     Power=power.IL.b.p$b.p,se=power.IL.b.p$b.p.se)
     Power.b.p = ggplot(data = data.IL.b.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of continuous predictor C on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.rho.Y,Power.a,Power.p,
     Power.b,Power.b.a,Power.b.p,ncol=1,heights=c(4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==23){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.rho.YF = data.frame(power.sim()$power[,c('N dyad','rho.YF','rho.YF.se')])
     data.IL.rho.YF = data.frame(Participants=c(paste(power.IL.rho.YF$N.dyad,sep=";")),
     Power=power.IL.rho.YF$rho.YF,se=power.IL.rho.YF$rho.YF.se)
     Power.rho.YF = ggplot(data = data.IL.rho.YF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A", x = "Number of dyads")

     power.IL.rho.YM = data.frame(power.sim()$power[,c('N dyad','rho.YM','rho.YM.se')])
     data.IL.rho.YM = data.frame(Participants=c(paste(power.IL.rho.YM$N.dyad,sep=";")),
     Power=power.IL.rho.YM$rho.YM,se=power.IL.rho.YM$rho.YM.se)
     Power.rho.YM = ggplot(data = data.IL.rho.YM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner A", x = "Number of dyads")

     power.IL.d.F = data.frame(power.sim()$power[,c('N dyad','d.F','d.F.se')])
     data.IL.d.F = data.frame(Participants=c(paste(power.IL.d.F$N.dyad,sep=";")),
     Power=power.IL.d.F$d.F,se=power.IL.d.F$d.F.se)
     Power.d.F = ggplot(data = data.IL.d.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner A", x = "Number of dyads")

     power.IL.d.FF = data.frame(power.sim()$power[,c('N dyad','d.FF','d.FF.se')])
     data.IL.d.FF = data.frame(Participants=c(paste(power.IL.d.FF$N.dyad,sep=";")),
     Power=power.IL.d.FF$d.FF,se=power.IL.d.FF$d.FF.se)
     Power.d.FF = ggplot(data = data.IL.d.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the actor effect of partner A", x = "Number of dyads")

     power.IL.d.MF = data.frame(power.sim()$power[,c('N dyad','d.MF','d.MF.se')])
     data.IL.d.MF = data.frame(Participants=c(paste(power.IL.d.MF$N.dyad,sep=";")),
     Power=power.IL.d.MF$d.MF,se=power.IL.d.MF$d.MF.se)
     Power.d.MF = ggplot(data = data.IL.d.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect of partner B", x = "Number of dyads")

     power.IL.d.M = data.frame(power.sim()$power[,c('N dyad','d.M','d.M.se')])
     data.IL.d.M = data.frame(Participants=c(paste(power.IL.d.M$N.dyad,sep=";")),
     Power=power.IL.d.M$d.M,se=power.IL.d.M$d.M.se)
     Power.d.M = ggplot(data = data.IL.d.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner B", x = "Number of dyads")

     power.IL.d.MM = data.frame(power.sim()$power[,c('N dyad','d.MM','d.MM.se')])
     data.IL.d.MM = data.frame(Participants=c(paste(power.IL.d.MM$N.dyad,sep=";")),
     Power=power.IL.d.MM$d.MM,se=power.IL.d.MM$d.MM.se)
     Power.d.MM = ggplot(data = data.IL.d.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the actor effect of partner B", x = "Number of dyads")

     power.IL.d.FM = data.frame(power.sim()$power[,c('N dyad','d.FM','d.FM.se')])
     data.IL.d.FM = data.frame(Participants=c(paste(power.IL.d.FM$N.dyad,sep=";")),
     Power=power.IL.d.FM$d.FM,se=power.IL.d.FM$d.FM.se)
     Power.d.FM = ggplot(data = data.IL.d.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.rho.YF,Power.rho.YM,Power.a.FF,Power.p.MF,
     Power.d.F,Power.d.FF,Power.d.MF,
     Power.a.MM,Power.p.FM,Power.d.M,Power.d.MM,Power.d.FM,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==24){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.rho.Y = data.frame(power.sim()$power[,c('N dyad','rho.Y','rho.Y.se')])
     data.IL.rho.Y = data.frame(Participants=c(paste(power.IL.rho.Y$N.dyad,sep=";")),
     Power=power.IL.rho.Y$rho.Y,se=power.IL.rho.Y$rho.Y.se)
     Power.rho.Y = ggplot(data = data.IL.rho.Y, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed partner effect", x = "Number of dyads")

     power.IL.d = data.frame(power.sim()$power[,c('N dyad','d','d.se')])
     data.IL.d = data.frame(Participants=c(paste(power.IL.d$N.dyad,sep=";")),
     Power=power.IL.d$d,se=power.IL.d$d.se)
     Power.d = ggplot(data = data.IL.d, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D", x = "Number of dyads")

     power.IL.d.a = data.frame(power.sim()$power[,c('N dyad','d.a','d.a.se')])
     data.IL.d.a = data.frame(Participants=c(paste(power.IL.d.a$N.dyad,sep=";")),
     Power=power.IL.d.a$d.a,se=power.IL.d.a$d.a.se)
     Power.d.a = ggplot(data = data.IL.d.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the actor effect", x = "Number of dyads")

     power.IL.d.p = data.frame(power.sim()$power[,c('N dyad','d.p','d.p.se')])
     data.IL.d.p = data.frame(Participants=c(paste(power.IL.d.p$N.dyad,sep=";")),
     Power=power.IL.d.p$d.p,se=power.IL.d.p$d.p.se)
     Power.d.p = ggplot(data = data.IL.d.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed moderation effect of dichotomous predictor D on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.rho.Y,Power.a,Power.p,
     Power.d,Power.d.a,Power.d.p,ncol=1,heights=c(4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==25){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.rho.YF = data.frame(power.sim()$power[,c('N dyad','rho.YF','rho.YF.se')])
     data.IL.rho.YF = data.frame(Participants=c(paste(power.IL.rho.YF$N.dyad,sep=";")),
     Power=power.IL.rho.YF$rho.YF,se=power.IL.rho.YF$rho.YF.se)
     Power.rho.YF = ggplot(data = data.IL.rho.YF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A", x = "Number of dyads")

     power.IL.rho.YM = data.frame(power.sim()$power[,c('N dyad','rho.YM','rho.YM.se')])
     data.IL.rho.YM = data.frame(Participants=c(paste(power.IL.rho.YM$N.dyad,sep=";")),
     Power=power.IL.rho.YM$rho.YM,se=power.IL.rho.YM$rho.YM.se)
     Power.rho.YM = ggplot(data = data.IL.rho.YM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A", x = "Number of dyads")

     power.IL.a.FF2 = data.frame(power.sim()$power[,c('N dyad','a.FF2','a.FF2.se')])
     data.IL.a.FF2 = data.frame(Participants=c(paste(power.IL.a.FF2$N.dyad,sep=";")),
     Power=power.IL.a.FF2$a.FF2,se=power.IL.a.FF2$a.FF2.se)
     Power.a.FF2 = ggplot(data = data.IL.a.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A", x = "Number of dyads")

     power.IL.p.MF2 = data.frame(power.sim()$power[,c('N dyad','p.MF2','p.MF2.se')])
     data.IL.p.MF2 = data.frame(Participants=c(paste(power.IL.p.MF2$N.dyad,sep=";")),
     Power=power.IL.p.MF2$p.MF2,se=power.IL.p.MF2$p.MF2.se)
     Power.p.MF2 = ggplot(data = data.IL.p.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B", x = "Number of dyads")

     power.IL.a.MM2 = data.frame(power.sim()$power[,c('N dyad','a.MM2','a.MM2.se')])
     data.IL.a.MM2 = data.frame(Participants=c(paste(power.IL.a.MM2$N.dyad,sep=";")),
     Power=power.IL.a.MM2$a.MM2,se=power.IL.a.MM2$a.MM2.se)
     Power.a.MM2 = ggplot(data = data.IL.a.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B", x = "Number of dyads")

     power.IL.p.FM2 = data.frame(power.sim()$power[,c('N dyad','p.FM2','p.FM2.se')])
     data.IL.p.FM2 = data.frame(Participants=c(paste(power.IL.p.FM2$N.dyad,sep=";")),
     Power=power.IL.p.FM2$p.FM2,se=power.IL.p.FM2$p.FM2.se)
     Power.p.FM2 = ggplot(data = data.IL.p.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.rho.YF,Power.rho.YM,Power.a.FF,Power.a.FF2,
     Power.p.MF,Power.p.MF2,Power.a.MM,Power.a.MM2,Power.p.FM,Power.p.FM2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==26){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.rho.Y = data.frame(power.sim()$power[,c('N dyad','rho.Y','rho.Y.se')])
     data.IL.rho.Y = data.frame(Participants=c(paste(power.IL.rho.Y$N.dyad,sep=";")),
     Power=power.IL.rho.Y$rho.Y,se=power.IL.rho.Y$rho.Y.se)
     Power.rho.Y = ggplot(data = data.IL.rho.Y, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect", x = "Number of dyads")

     power.IL.a.2 = data.frame(power.sim()$power[,c('N dyad','a.2','a.2.se')])
     data.IL.a.2 = data.frame(Participants=c(paste(power.IL.a.2$N.dyad,sep=";")),
     Power=power.IL.a.2$a.2,se=power.IL.a.2$a.2.se)
     Power.a.2 = ggplot(data = data.IL.a.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect", x = "Number of dyads")

     power.IL.p.2 = data.frame(power.sim()$power[,c('N dyad','p.2','p.2.se')])
     data.IL.p.2 = data.frame(Participants=c(paste(power.IL.p.2$N.dyad,sep=";")),
     Power=power.IL.p.2$p.2,se=power.IL.p.2$p.2.se)
     Power.p.2 = ggplot(data = data.IL.p.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.rho.Y,Power.a,Power.a.2,Power.p,Power.p.2,ncol=1,heights=c(4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==27){
     power.IL.c.F0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F0','c.F0.se')])
     data.IL.c.F0 = data.frame(Participants=c(paste(power.IL.c.F0$N.Dyad.Group.0.,power.IL.c.F0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F0$c.F0,se=power.IL.c.F0$c.F0.se)
     Power.c.F0 = ggplot(data = data.IL.c.F0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A in Group 0", x = "Number of dyads")

     power.IL.c.F1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.F1','c.F1.se')])
     data.IL.c.F1 = data.frame(Participants=c(paste(power.IL.c.F1$N.Dyad.Group.0.,power.IL.c.F1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.F1$c.F1,se=power.IL.c.F1$c.F1.se)
     Power.c.F1 = ggplot(data = data.IL.c.F1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.c.M0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M0','c.M0.se')])
     data.IL.c.M0 = data.frame(Participants=c(paste(power.IL.c.M0$N.Dyad.Group.0.,power.IL.c.M0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M0$c.M0,se=power.IL.c.M0$c.M0.se)
     Power.c.M0 = ggplot(data = data.IL.c.M0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed intercept of partner B in Group 0", x = "Number of dyads")

     power.IL.c.M1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c.M1','c.M1.se')])
     data.IL.c.M1 = data.frame(Participants=c(paste(power.IL.c.M1$N.Dyad.Group.0.,power.IL.c.M1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c.M1$c.M1,se=power.IL.c.M1$c.M1.se)
     Power.c.M1 = ggplot(data = data.IL.c.M1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the linear fixed intercept of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.rho.YF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YF0','rho.YF0.se')])
     data.IL.rho.YF0 = data.frame(Participants=c(paste(power.IL.rho.YF0$N.Dyad.Group.0.,power.IL.rho.YF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YF0$rho.YF0,se=power.IL.rho.YF0$rho.YF0.se)
     Power.rho.YF0 = ggplot(data = data.IL.rho.YF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A in Group 0", x = "Number of dyads")

     power.IL.rho.YF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YF1','rho.YF1.se')])
     data.IL.rho.YF1 = data.frame(Participants=c(paste(power.IL.rho.YF1$N.Dyad.Group.0.,power.IL.rho.YF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YF1$rho.YF1,se=power.IL.rho.YF1$rho.YF1.se)
     Power.rho.YF1 = ggplot(data = data.IL.rho.YF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed autoregressive effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.rho.YM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YM0','rho.YM0.se')])
     data.IL.rho.YM0 = data.frame(Participants=c(paste(power.IL.rho.YM0$N.Dyad.Group.0.,power.IL.rho.YM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YM0$rho.YM0,se=power.IL.rho.YM0$rho.YM0.se)
     Power.rho.YM0 = ggplot(data = data.IL.rho.YM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B in Group 0", x = "Number of dyads")

     power.IL.rho.YM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.YM1','rho.YM1.se')])
     data.IL.rho.YM1 = data.frame(Participants=c(paste(power.IL.rho.YM1$N.Dyad.Group.0.,power.IL.rho.YM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.YM1$rho.YM1,se=power.IL.rho.YM1$rho.YM1.se)
     Power.rho.YM1 = ggplot(data = data.IL.rho.YM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed autoregressive effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.a.FF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF0','a.FF0.se')])
     data.IL.a.FF0 = data.frame(Participants=c(paste(power.IL.a.FF0$N.Dyad.Group.0.,power.IL.a.FF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF0$a.FF0,se=power.IL.a.FF0$a.FF0.se)
     Power.a.FF0 = ggplot(data = data.IL.a.FF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A in Group 0", x = "Number of dyads")

     power.IL.a.FF02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF02','a.FF02.se')])
     data.IL.a.FF02 = data.frame(Participants=c(paste(power.IL.a.FF02$N.Dyad.Group.0.,power.IL.a.FF02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF02$a.FF02,se=power.IL.a.FF02$a.FF02.se)
     Power.a.FF02 = ggplot(data = data.IL.a.FF02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A in Group 0", x = "Number of dyads")

     power.IL.a.FF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF1','a.FF1.se')])
     data.IL.a.FF1 = data.frame(Participants=c(paste(power.IL.a.FF1$N.Dyad.Group.0.,power.IL.a.FF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF1$a.FF1,se=power.IL.a.FF1$a.FF1.se)
     Power.a.FF1 = ggplot(data = data.IL.a.FF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed actor effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.a.FF12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.FF12','a.FF12.se')])
     data.IL.a.FF12 = data.frame(Participants=c(paste(power.IL.a.FF12$N.Dyad.Group.0.,power.IL.a.FF12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.FF12$a.FF12,se=power.IL.a.FF12$a.FF12.se)
     Power.a.FF12 = ggplot(data = data.IL.a.FF12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed actor effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.p.MF0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF0','p.MF0.se')])
     data.IL.p.MF0 = data.frame(Participants=c(paste(power.IL.p.MF0$N.Dyad.Group.0.,power.IL.p.MF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF0$p.MF0,se=power.IL.p.MF0$p.MF0.se)
     Power.p.MF0 = ggplot(data = data.IL.p.MF0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A in Group 0", x = "Number of dyads")

     power.IL.p.MF02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF02','p.MF02.se')])
     data.IL.p.MF02 = data.frame(Participants=c(paste(power.IL.p.MF02$N.Dyad.Group.0.,power.IL.p.MF02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF02$p.MF02,se=power.IL.p.MF02$p.MF02.se)
     Power.p.MF02 = ggplot(data = data.IL.p.MF02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A in Group 0", x = "Number of dyads")

     power.IL.p.MF1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF1','p.MF1.se')])
     data.IL.p.MF1 = data.frame(Participants=c(paste(power.IL.p.MF1$N.Dyad.Group.0.,power.IL.p.MF1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF1$p.MF1,se=power.IL.p.MF1$p.MF1.se)
     Power.p.MF1 = ggplot(data = data.IL.p.MF1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.p.MF12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.MF12','p.MF12.se')])
     data.IL.p.MF12 = data.frame(Participants=c(paste(power.IL.p.MF12$N.Dyad.Group.0.,power.IL.p.MF12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.MF12$p.MF12,se=power.IL.p.MF12$p.MF12.se)
     Power.p.MF12 = ggplot(data = data.IL.p.MF12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed partner effect of partner A between Group 0 and 1", x = "Number of dyads")

     power.IL.a.MM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM0','a.MM0.se')])
     data.IL.a.MM0 = data.frame(Participants=c(paste(power.IL.a.MM0$N.Dyad.Group.0.,power.IL.a.MM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM0$a.MM0,se=power.IL.a.MM0$a.MM0.se)
     Power.a.MM0 = ggplot(data = data.IL.a.MM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B in Group 0", x = "Number of dyads")

     power.IL.a.MM02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM02','a.MM02.se')])
     data.IL.a.MM02 = data.frame(Participants=c(paste(power.IL.a.MM02$N.Dyad.Group.0.,power.IL.a.MM02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM02$a.MM02,se=power.IL.a.MM02$a.MM02.se)
     Power.a.MM02 = ggplot(data = data.IL.a.MM02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B in Group 0", x = "Number of dyads")

     power.IL.a.MM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM1','a.MM1.se')])
     data.IL.a.MM1 = data.frame(Participants=c(paste(power.IL.a.MM1$N.Dyad.Group.0.,power.IL.a.MM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM1$a.MM1,se=power.IL.a.MM1$a.MM1.se)
     Power.a.MM1 = ggplot(data = data.IL.a.MM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed actor effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.a.MM12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a.MM12','a.MM12.se')])
     data.IL.a.MM12 = data.frame(Participants=c(paste(power.IL.a.MM12$N.Dyad.Group.0.,power.IL.a.MM12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a.MM12$a.MM12,se=power.IL.a.MM12$a.MM12.se)
     Power.a.MM12 = ggplot(data = data.IL.a.MM12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed actor effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.p.FM0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM0','p.FM0.se')])
     data.IL.p.FM0 = data.frame(Participants=c(paste(power.IL.p.FM0$N.Dyad.Group.0.,power.IL.p.FM0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM0$p.FM0,se=power.IL.p.FM0$p.FM0.se)
     Power.p.FM0 = ggplot(data = data.IL.p.FM0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B in Group 0", x = "Number of dyads")

     power.IL.p.FM02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM02','p.FM02.se')])
     data.IL.p.FM02 = data.frame(Participants=c(paste(power.IL.p.FM02$N.Dyad.Group.0.,power.IL.p.FM02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM02$p.FM02,se=power.IL.p.FM02$p.FM02.se)
     Power.p.FM02 = ggplot(data = data.IL.p.FM02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B in Group 0", x = "Number of dyads")

     power.IL.p.FM1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM1','p.FM1.se')])
     data.IL.p.FM1 = data.frame(Participants=c(paste(power.IL.p.FM1$N.Dyad.Group.0.,power.IL.p.FM1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM1$p.FM1,se=power.IL.p.FM1$p.FM1.se)
     Power.p.FM1 = ggplot(data = data.IL.p.FM1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect of partner B between Group 0 and 1", x = "Number of dyads")

     power.IL.p.FM12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p.FM12','p.FM12.se')])
     data.IL.p.FM12 = data.frame(Participants=c(paste(power.IL.p.FM12$N.Dyad.Group.0.,power.IL.p.FM12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p.FM12$p.FM12,se=power.IL.p.FM12$p.FM12.se)
     Power.p.FM12 = ggplot(data = data.IL.p.FM12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed partner effect of partner B between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F0,Power.c.F1,Power.c.M0,Power.c.M1,
     Power.rho.YF0,Power.rho.YF1,Power.rho.YM0,Power.rho.YM1,
     Power.a.FF0,Power.a.FF02,Power.a.FF1,Power.a.FF12,
     Power.p.MF0,Power.p.MF02,Power.p.MF1,Power.p.MF12,
     Power.a.MM0,Power.a.MM02,Power.a.MM1,Power.a.MM12,
     Power.p.FM0,Power.p.FM02,Power.p.MF1,Power.p.MF12,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N0.dyad,","))))!=1){
     if (input$Model==28){
     power.IL.c0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c0','c0.se')])
     data.IL.c0 = data.frame(Participants=c(paste(power.IL.c0$N.Dyad.Group.0.,power.IL.c0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c0$c0,se=power.IL.c0$c0.se)
     Power.c0 = ggplot(data = data.IL.c0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept in Group 0", x = "Number of dyads")

     power.IL.c1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','c1','c1.se')])
     data.IL.c1 = data.frame(Participants=c(paste(power.IL.c1$N.Dyad.Group.0.,power.IL.c1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.c1$c1,se=power.IL.c1$c1.se)
     Power.c1 = ggplot(data = data.IL.c1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed intercept between Group 0 and 1", x = "Number of dyads")

     power.IL.rho.Y0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.Y0','rho.Y0.se')])
     data.IL.rho.Y0 = data.frame(Participants=c(paste(power.IL.rho.Y0$N.Dyad.Group.0.,power.IL.rho.YF0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.Y0$rho.Y0,se=power.IL.rho.Y0$rho.Y0.se)
     Power.rho.Y0 = ggplot(data = data.IL.rho.Y0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect in Group 0", x = "Number of dyads")

     power.IL.rho.Y1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','rho.Y1','rho.Y1.se')])
     data.IL.rho.Y1 = data.frame(Participants=c(paste(power.IL.rho.Y1$N.Dyad.Group.0.,power.IL.rho.Y1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.rho.Y1$rho.Y1,se=power.IL.rho.Y1$rho.Y1.se)
     Power.rho.Y1 = ggplot(data = data.IL.rho.Y1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the difference in the fixed autoregressive effect between Group 0 and 1", x = "Number of dyads")

     power.IL.a0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a0','a0.se')])
     data.IL.a0 = data.frame(Participants=c(paste(power.IL.a0$N.Dyad.Group.0.,power.IL.a0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a0$a0,se=power.IL.a0$a0.se)
     Power.a0 = ggplot(data = data.IL.a0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect F in Group 0", x = "Number of dyads")

     power.IL.a02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a02','a02.se')])
     data.IL.a02 = data.frame(Participants=c(paste(power.IL.a02$N.Dyad.Group.0.,power.IL.a02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a02$a02,se=power.IL.a02$a02.se)
     Power.a02 = ggplot(data = data.IL.a02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect F in Group 0", x = "Number of dyads")

     power.IL.a1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a1','a1.se')])
     data.IL.a1 = data.frame(Participants=c(paste(power.IL.a1$N.Dyad.Group.0.,power.IL.a1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a1$a1,se=power.IL.a1$a1.se)
     Power.a1 = ggplot(data = data.IL.a1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed actor effect between Group 0 and 1", x = "Number of dyads")

     power.IL.a12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','a12','a12.se')])
     data.IL.a12 = data.frame(Participants=c(paste(power.IL.a12$N.Dyad.Group.0.,power.IL.a12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.a12$a12,se=power.IL.a12$a12.se)
     Power.a12 = ggplot(data = data.IL.a12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the quadratic fixed actor effect between Group 0 and 1", x = "Number of dyads")

     power.IL.p0 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p0','p0.se')])
     data.IL.p0 = data.frame(Participants=c(paste(power.IL.p0$N.Dyad.Group.0.,power.IL.p0$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p0$p0,se=power.IL.p0$p0.se)
     Power.p0 = ggplot(data = data.IL.p0, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p0$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect in Group 0", x = "Number of dyads")

     power.IL.p02 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p02','p02.se')])
     data.IL.p02 = data.frame(Participants=c(paste(power.IL.p02$N.Dyad.Group.0.,power.IL.p02$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p02$p02,se=power.IL.p02$p02.se)
     Power.p02 = ggplot(data = data.IL.p02, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p02$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect in Group 0", x = "Number of dyads")

     power.IL.p1 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p1','p1.se')])
     data.IL.p1 = data.frame(Participants=c(paste(power.IL.p1$N.Dyad.Group.0.,power.IL.p1$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p1$p1,se=power.IL.p1$p1.se)
     Power.p1 = ggplot(data = data.IL.p1, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p1$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect between Group 0 and 1", x = "Number of dyads")

     power.IL.p12 = data.frame(power.sim()$power[,c('N Dyad(Group=0)','N Dyad(Group=1)','p12','p12.se')])
     data.IL.p12 = data.frame(Participants=c(paste(power.IL.p12$N.Dyad.Group.0.,power.IL.p12$N.Dyad.Group.1.,sep=";")),
     Power=power.IL.p12$p12,se=power.IL.p12$p12.se)
     Power.p12 = ggplot(data = data.IL.p12, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p12$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the differences in the linear fixed partner effect between Group 0 and 1", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c0,Power.c1,Power.rho.Y0,Power.rho.Y1,
     Power.a0,Power.a02,Power.a1,Power.a12,
     Power.p0,Power.p02,Power.p1,Power.p12,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==29){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B", x = "Number of dyads")

     power.IL.rho.YF = data.frame(power.sim()$power[,c('N dyad','rho.YF','rho.YF.se')])
     data.IL.rho.YF = data.frame(Participants=c(paste(power.IL.rho.YF$N.dyad,sep=";")),
     Power=power.IL.rho.YF$rho.YF,se=power.IL.rho.YF$rho.YF.se)
     Power.rho.YF = ggplot(data = data.IL.rho.YF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A", x = "Number of dyads")

     power.IL.rho.YM = data.frame(power.sim()$power[,c('N dyad','rho.YM','rho.YM.se')])
     data.IL.rho.YM = data.frame(Participants=c(paste(power.IL.rho.YM$N.dyad,sep=";")),
     Power=power.IL.rho.YM$rho.YM,se=power.IL.rho.YM$rho.YM.se)
     Power.rho.YM = ggplot(data = data.IL.rho.YM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A", x = "Number of dyads")

     power.IL.a.FF2 = data.frame(power.sim()$power[,c('N dyad','a.FF2','a.FF2.se')])
     data.IL.a.FF2 = data.frame(Participants=c(paste(power.IL.a.FF2$N.dyad,sep=";")),
     Power=power.IL.a.FF2$a.FF2,se=power.IL.a.FF2$a.FF2.se)
     Power.a.FF2 = ggplot(data = data.IL.a.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A", x = "Number of dyads")

     power.IL.p.MF2 = data.frame(power.sim()$power[,c('N dyad','p.MF2','p.MF2.se')])
     data.IL.p.MF2 = data.frame(Participants=c(paste(power.IL.p.MF2$N.dyad,sep=";")),
     Power=power.IL.p.MF2$p.MF2,se=power.IL.p.MF2$p.MF2.se)
     Power.p.MF2 = ggplot(data = data.IL.p.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A", x = "Number of dyads")

     power.IL.b.F = data.frame(power.sim()$power[,c('N dyad','b.F','b.F.se')])
     data.IL.b.F = data.frame(Participants=c(paste(power.IL.b.F$N.dyad,sep=";")),
     Power=power.IL.b.F$b.F,se=power.IL.b.F$b.F.se)
     Power.b.F = ggplot(data = data.IL.b.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner A", x = "Number of dyads")

     power.IL.b.FF = data.frame(power.sim()$power[,c('N dyad','b.FF','b.FF.se')])
     data.IL.b.FF = data.frame(Participants=c(paste(power.IL.b.FF$N.dyad,sep=";")),
     Power=power.IL.b.FF$b.FF,se=power.IL.b.FF$b.FF.se)
     Power.b.FF = ggplot(data = data.IL.b.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the actor effect of partner A", x = "Number of dyads")

     power.IL.b.FF2 = data.frame(power.sim()$power[,c('N dyad','b.FF2','b.FF2.se')])
     data.IL.b.FF2 = data.frame(Participants=c(paste(power.IL.b.FF2$N.dyad,sep=";")),
     Power=power.IL.b.FF2$b.FF2,se=power.IL.b.FF2$b.FF2.se)
     Power.b.FF2 = ggplot(data = data.IL.b.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the actor effect of partner A", x = "Number of dyads")

     power.IL.b.MF = data.frame(power.sim()$power[,c('N dyad','b.MF','b.MF.se')])
     data.IL.b.MF = data.frame(Participants=c(paste(power.IL.b.MF$N.dyad,sep=";")),
     Power=power.IL.b.MF$b.MF,se=power.IL.b.MF$b.MF.se)
     Power.b.MF = ggplot(data = data.IL.b.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the partner effect of partner A", x = "Number of dyads")

     power.IL.b.MF2 = data.frame(power.sim()$power[,c('N dyad','b.MF2','b.MF2.se')])
     data.IL.b.MF2 = data.frame(Participants=c(paste(power.IL.b.MF2$N.dyad,sep=";")),
     Power=power.IL.b.MF2$b.MF2,se=power.IL.b.MF2$b.MF2.se)
     Power.b.MF2 = ggplot(data = data.IL.b.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B", x = "Number of dyads")

     power.IL.a.MM2 = data.frame(power.sim()$power[,c('N dyad','a.MM2','a.MM2.se')])
     data.IL.a.MM2 = data.frame(Participants=c(paste(power.IL.a.MM2$N.dyad,sep=";")),
     Power=power.IL.a.MM2$a.MM2,se=power.IL.a.MM2$a.MM2.se)
     Power.a.MM2 = ggplot(data = data.IL.a.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B", x = "Number of dyads")

     power.IL.p.FM2 = data.frame(power.sim()$power[,c('N dyad','p.FM2','p.FM2.se')])
     data.IL.p.FM2 = data.frame(Participants=c(paste(power.IL.p.FM2$N.dyad,sep=";")),
     Power=power.IL.p.FM2$p.FM2,se=power.IL.p.FM2$p.FM2.se)
     Power.p.FM2 = ggplot(data = data.IL.p.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B", x = "Number of dyads")

     power.IL.b.M = data.frame(power.sim()$power[,c('N dyad','b.M','b.M.se')])
     data.IL.b.M = data.frame(Participants=c(paste(power.IL.b.M$N.dyad,sep=";")),
     Power=power.IL.b.M$b.M,se=power.IL.b.M$b.M.se)
     Power.b.M = ggplot(data = data.IL.b.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C of partner B", x = "Number of dyads")

     power.IL.b.MM = data.frame(power.sim()$power[,c('N dyad','b.MM','b.MM.se')])
     data.IL.b.MM = data.frame(Participants=c(paste(power.IL.b.MM$N.dyad,sep=";")),
     Power=power.IL.b.MM$b.MM,se=power.IL.b.MM$b.MM.se)
     Power.b.MM = ggplot(data = data.IL.b.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the actor effect of partner B", x = "Number of dyads")

     power.IL.b.MM2 = data.frame(power.sim()$power[,c('N dyad','b.MM2','b.MM2.se')])
     data.IL.b.MM2 = data.frame(Participants=c(paste(power.IL.b.MM2$N.dyad,sep=";")),
     Power=power.IL.b.MM2$b.MM2,se=power.IL.b.MM2$b.MM2.se)
     Power.b.MM2 = ggplot(data = data.IL.b.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the actor effect of partner B", x = "Number of dyads")

     power.IL.b.FM = data.frame(power.sim()$power[,c('N dyad','b.FM','b.FM.se')])
     data.IL.b.FM = data.frame(Participants=c(paste(power.IL.b.FM$N.dyad,sep=";")),
     Power=power.IL.b.FM$b.FM,se=power.IL.b.FM$b.FM.se)
     Power.b.FM = ggplot(data = data.IL.b.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the partner effect of partner B", x = "Number of dyads")

     power.IL.b.FM2 = data.frame(power.sim()$power[,c('N dyad','b.FM2','b.FM2.se')])
     data.IL.b.FM2 = data.frame(Participants=c(paste(power.IL.b.FM2$N.dyad,sep=";")),
     Power=power.IL.b.FM2$b.FM2,se=power.IL.b.FM2$b.FM2.se)
     Power.b.FM2 = ggplot(data = data.IL.b.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.rho.YF,Power.rho.YM,Power.a.FF,Power.a.FF2,
     Power.p.MF,Power.p.MF2,Power.b.F,
     Power.b.FF,Power.b.FF2,Power.b.MF,Power.b.MF2,
     Power.a.MM,Power.a.MM2,Power.p.FM,Power.p.FM2,Power.b.M,
     Power.b.MM,Power.b.MM2,Power.b.FM,Power.b.FM2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==30){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.rho.Y = data.frame(power.sim()$power[,c('N dyad','rho.Y','rho.Y.se')])
     data.IL.rho.Y = data.frame(Participants=c(paste(power.IL.rho.Y$N.dyad,sep=";")),
     Power=power.IL.rho.Y$rho.Y,se=power.IL.rho.Y$rho.Y.se)
     Power.rho.Y = ggplot(data = data.IL.rho.Y, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect", x = "Number of dyads")

     power.IL.a.2 = data.frame(power.sim()$power[,c('N dyad','a.2','a.2.se')])
     data.IL.a.2 = data.frame(Participants=c(paste(power.IL.a.2$N.dyad,sep=";")),
     Power=power.IL.a.2$a.2,se=power.IL.a.2$a.2.se)
     Power.a.2 = ggplot(data = data.IL.a.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect", x = "Number of dyads")

     power.IL.p.2 = data.frame(power.sim()$power[,c('N dyad','p.2','p.2.se')])
     data.IL.p.2 = data.frame(Participants=c(paste(power.IL.p.2$N.dyad,sep=";")),
     Power=power.IL.p.2$p.2,se=power.IL.p.2$p.2.se)
     Power.p.2 = ggplot(data = data.IL.p.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect", x = "Number of dyads")

     power.IL.b = data.frame(power.sim()$power[,c('N dyad','b','b.se')])
     data.IL.b = data.frame(Participants=c(paste(power.IL.b$N.dyad,sep=";")),
     Power=power.IL.b$b,se=power.IL.b$b.se)
     Power.b = ggplot(data = data.IL.b, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor C", x = "Number of dyads")

     power.IL.b.a = data.frame(power.sim()$power[,c('N dyad','b.a','b.a.se')])
     data.IL.b.a = data.frame(Participants=c(paste(power.IL.b.a$N.dyad,sep=";")),
     Power=power.IL.b.a$b.a,se=power.IL.b.a$b.a.se)
     Power.b.a = ggplot(data = data.IL.b.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the actor effect", x = "Number of dyads")

     power.IL.b.a2 = data.frame(power.sim()$power[,c('N dyad','b.a2','b.a2.se')])
     data.IL.b.a2 = data.frame(Participants=c(paste(power.IL.b.a2$N.dyad,sep=";")),
     Power=power.IL.b.a2$b.a2,se=power.IL.b.a2$b.a2.se)
     Power.b.a2 = ggplot(data = data.IL.b.a2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.a2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the actor effect", x = "Number of dyads")

     power.IL.b.p = data.frame(power.sim()$power[,c('N dyad','b.p','b.p.se')])
     data.IL.b.p = data.frame(Participants=c(paste(power.IL.b.p$N.dyad,sep=";")),
     Power=power.IL.b.p$b.p,se=power.IL.b.p$b.p.se)
     Power.b.p = ggplot(data = data.IL.b.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of continuous predictor C on the partner effect", x = "Number of dyads")

     power.IL.b.p2 = data.frame(power.sim()$power[,c('N dyad','b.p2','b.p2.se')])
     data.IL.b.p2 = data.frame(Participants=c(paste(power.IL.b.p2$N.dyad,sep=";")),
     Power=power.IL.b.p2$b.p2,se=power.IL.b.p2$b.p2.se)
     Power.b.p2 = ggplot(data = data.IL.b.p2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.b.p2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of continuous predictor C on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.rho.Y,Power.a,Power.a.2,Power.p,Power.p.2,
     Power.b,Power.b.a,Power.b.a2,Power.b.p,Power.b.p2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==31){
     power.IL.c.F = data.frame(power.sim()$power[,c('N dyad','c.F','c.F.se')])
     data.IL.c.F = data.frame(Participants=c(paste(power.IL.c.F$N.dyad,sep=";")),
     Power=power.IL.c.F$c.F,se=power.IL.c.F$c.F.se)
     Power.c.F = ggplot(data = data.IL.c.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner A", x = "Number of dyads")

     power.IL.c.M = data.frame(power.sim()$power[,c('N dyad','c.M','c.M.se')])
     data.IL.c.M = data.frame(Participants=c(paste(power.IL.c.M$N.dyad,sep=";")),
     Power=power.IL.c.M$c.M,se=power.IL.c.M$c.M.se)
     Power.c.M = ggplot(data = data.IL.c.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept of partner B")

     power.IL.rho.YF = data.frame(power.sim()$power[,c('N dyad','rho.YF','rho.YF.se')])
     data.IL.rho.YF = data.frame(Participants=c(paste(power.IL.rho.YF$N.dyad,sep=";")),
     Power=power.IL.rho.YF$rho.YF,se=power.IL.rho.YF$rho.YF.se)
     Power.rho.YF = ggplot(data = data.IL.rho.YF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner A", x = "Number of dyads")

     power.IL.rho.YM = data.frame(power.sim()$power[,c('N dyad','rho.YM','rho.YM.se')])
     data.IL.rho.YM = data.frame(Participants=c(paste(power.IL.rho.YM$N.dyad,sep=";")),
     Power=power.IL.rho.YM$rho.YM,se=power.IL.rho.YM$rho.YM.se)
     Power.rho.YM = ggplot(data = data.IL.rho.YM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.YM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect of partner B", x = "Number of dyads")

     power.IL.a.FF = data.frame(power.sim()$power[,c('N dyad','a.FF','a.FF.se')])
     data.IL.a.FF = data.frame(Participants=c(paste(power.IL.a.FF$N.dyad,sep=";")),
     Power=power.IL.a.FF$a.FF,se=power.IL.a.FF$a.FF.se)
     Power.a.FF = ggplot(data = data.IL.a.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner A", x = "Number of dyads")

     power.IL.a.FF2 = data.frame(power.sim()$power[,c('N dyad','a.FF2','a.FF2.se')])
     data.IL.a.FF2 = data.frame(Participants=c(paste(power.IL.a.FF2$N.dyad,sep=";")),
     Power=power.IL.a.FF2$a.FF2,se=power.IL.a.FF2$a.FF2.se)
     Power.a.FF2 = ggplot(data = data.IL.a.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner A", x = "Number of dyads")

     power.IL.p.MF = data.frame(power.sim()$power[,c('N dyad','p.MF','p.MF.se')])
     data.IL.p.MF = data.frame(Participants=c(paste(power.IL.p.MF$N.dyad,sep=";")),
     Power=power.IL.p.MF$p.MF,se=power.IL.p.MF$p.MF.se)
     Power.p.MF = ggplot(data = data.IL.p.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner A", x = "Number of dyads")

     power.IL.p.MF2 = data.frame(power.sim()$power[,c('N dyad','p.MF2','p.MF2.se')])
     data.IL.p.MF2 = data.frame(Participants=c(paste(power.IL.p.MF2$N.dyad,sep=";")),
     Power=power.IL.p.MF2$p.MF2,se=power.IL.p.MF2$p.MF2.se)
     Power.p.MF2 = ggplot(data = data.IL.p.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner A", x = "Number of dyads")

     power.IL.d.F = data.frame(power.sim()$power[,c('N dyad','d.F','d.F.se')])
     data.IL.d.F = data.frame(Participants=c(paste(power.IL.d.F$N.dyad,sep=";")),
     Power=power.IL.d.F$d.F,se=power.IL.d.F$d.F.se)
     Power.d.F = ggplot(data = data.IL.d.F, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.F$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner A", x = "Number of dyads")

     power.IL.d.FF = data.frame(power.sim()$power[,c('N dyad','d.FF','d.FF.se')])
     data.IL.d.FF = data.frame(Participants=c(paste(power.IL.d.FF$N.dyad,sep=";")),
     Power=power.IL.d.FF$d.FF,se=power.IL.d.FF$d.FF.se)
     Power.d.FF = ggplot(data = data.IL.d.FF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the actor effect of partner A", x = "Number of dyads")

     power.IL.d.FF2 = data.frame(power.sim()$power[,c('N dyad','d.FF2','d.FF2.se')])
     data.IL.d.FF2 = data.frame(Participants=c(paste(power.IL.d.FF2$N.dyad,sep=";")),
     Power=power.IL.d.FF2$d.FF2,se=power.IL.d.FF2$d.FF2.se)
     Power.d.FF2 = ggplot(data = data.IL.d.FF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the actor effect of partner A", x = "Number of dyads")

     power.IL.d.MF = data.frame(power.sim()$power[,c('N dyad','d.MF','d.MF.se')])
     data.IL.d.MF = data.frame(Participants=c(paste(power.IL.d.MF$N.dyad,sep=";")),
     Power=power.IL.d.MF$d.MF,se=power.IL.d.MF$d.MF.se)
     Power.d.MF = ggplot(data = data.IL.d.MF, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MF$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the partner effect of partner A", x = "Number of dyads")

     power.IL.d.MF2 = data.frame(power.sim()$power[,c('N dyad','d.MF2','d.MF2.se')])
     data.IL.d.MF2 = data.frame(Participants=c(paste(power.IL.d.MF2$N.dyad,sep=";")),
     Power=power.IL.d.MF2$d.MF2,se=power.IL.d.MF2$d.MF2.se)
     Power.d.MF2 = ggplot(data = data.IL.d.MF2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MF2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the partner effect of partner A", x = "Number of dyads")

     power.IL.a.MM = data.frame(power.sim()$power[,c('N dyad','a.MM','a.MM.se')])
     data.IL.a.MM = data.frame(Participants=c(paste(power.IL.a.MM$N.dyad,sep=";")),
     Power=power.IL.a.MM$a.MM,se=power.IL.a.MM$a.MM.se)
     Power.a.MM = ggplot(data = data.IL.a.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect of partner B", x = "Number of dyads")

     power.IL.a.MM2 = data.frame(power.sim()$power[,c('N dyad','a.MM2','a.MM2.se')])
     data.IL.a.MM2 = data.frame(Participants=c(paste(power.IL.a.MM2$N.dyad,sep=";")),
     Power=power.IL.a.MM2$a.MM2,se=power.IL.a.MM2$a.MM2.se)
     Power.a.MM2 = ggplot(data = data.IL.a.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect of partner B", x = "Number of dyads")

     power.IL.p.FM = data.frame(power.sim()$power[,c('N dyad','p.FM','p.FM.se')])
     data.IL.p.FM = data.frame(Participants=c(paste(power.IL.p.FM$N.dyad,sep=";")),
     Power=power.IL.p.FM$p.FM,se=power.IL.p.FM$p.FM.se)
     Power.p.FM = ggplot(data = data.IL.p.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect of partner B", x = "Number of dyads")

     power.IL.p.FM2 = data.frame(power.sim()$power[,c('N dyad','p.FM2','p.FM2.se')])
     data.IL.p.FM2 = data.frame(Participants=c(paste(power.IL.p.FM2$N.dyad,sep=";")),
     Power=power.IL.p.FM2$p.FM2,se=power.IL.p.FM2$p.FM2.se)
     Power.p.FM2 = ggplot(data = data.IL.p.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect of partner B", x = "Number of dyads")

     power.IL.d.M = data.frame(power.sim()$power[,c('N dyad','d.M','d.M.se')])
     data.IL.d.M = data.frame(Participants=c(paste(power.IL.d.M$N.dyad,sep=";")),
     Power=power.IL.d.M$d.M,se=power.IL.d.M$d.M.se)
     Power.d.M = ggplot(data = data.IL.d.M, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.M$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of dichotomous predictor D of partner B", x = "Number of dyads")

     power.IL.d.MM = data.frame(power.sim()$power[,c('N dyad','d.MM','d.MM.se')])
     data.IL.d.MM = data.frame(Participants=c(paste(power.IL.d.MM$N.dyad,sep=";")),
     Power=power.IL.d.MM$d.MM,se=power.IL.d.MM$d.MM.se)
     Power.d.MM = ggplot(data = data.IL.d.MM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the actor effect of partner B", x = "Number of dyads")

     power.IL.d.MM2 = data.frame(power.sim()$power[,c('N dyad','d.MM2','d.MM2.se')])
     data.IL.d.MM2 = data.frame(Participants=c(paste(power.IL.d.MM2$N.dyad,sep=";")),
     Power=power.IL.d.MM2$d.MM2,se=power.IL.d.MM2$d.MM2.se)
     Power.d.MM2 = ggplot(data = data.IL.d.MM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.MM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the actor effect of partner B", x = "Number of dyads")

     power.IL.d.FM = data.frame(power.sim()$power[,c('N dyad','d.FM','d.FM.se')])
     data.IL.d.FM = data.frame(Participants=c(paste(power.IL.d.FM$N.dyad,sep=";")),
     Power=power.IL.d.FM$d.FM,se=power.IL.d.FM$d.FM.se)
     Power.d.FM = ggplot(data = data.IL.d.FM, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FM$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the partner effect of partner B", x = "Number of dyads")

     power.IL.d.FM2 = data.frame(power.sim()$power[,c('N dyad','d.FM2','d.FM2.se')])
     data.IL.d.FM2 = data.frame(Participants=c(paste(power.IL.d.FM2$N.dyad,sep=";")),
     Power=power.IL.d.FM2$d.FM2,se=power.IL.d.FM2$d.FM2.se)
     Power.d.FM2 = ggplot(data = data.IL.d.FM2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.FM2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the partner effect of partner B", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c.F,Power.c.M,Power.rho.YF,Power.rho.YM,Power.a.FF,Power.a.FF2,
     Power.p.MF,Power.p.MF2,Power.d.F,
     Power.d.FF,Power.d.FF2,Power.d.MF,Power.d.MF2,
     Power.a.MM,Power.a.MM2,Power.p.FM,Power.p.FM2,
     Power.d.M,Power.d.MM,Power.d.MM2,Power.d.FM,Power.d.FM2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
     }}

     if (length(as.numeric(unlist(strsplit(input$N.dyad,","))))!=1){
     if (input$Model==32){
     power.IL.c = data.frame(power.sim()$power[,c('N dyad','c','c.se')])
     data.IL.c = data.frame(Participants=c(paste(power.IL.c$N.dyad,sep=";")),
     Power=power.IL.c$c,se=power.IL.c$c.se)
     Power.c = ggplot(data = data.IL.c, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.c$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed intercept", x = "Number of dyads")

     power.IL.rho.Y = data.frame(power.sim()$power[,c('N dyad','rho.Y','rho.Y.se')])
     data.IL.rho.Y = data.frame(Participants=c(paste(power.IL.rho.Y$N.dyad,sep=";")),
     Power=power.IL.rho.Y$rho.Y,se=power.IL.rho.Y$rho.Y.se)
     Power.rho.Y = ggplot(data = data.IL.rho.Y, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.rho.Y$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed autoregressive effect", x = "Number of dyads")

     power.IL.a = data.frame(power.sim()$power[,c('N dyad','a','a.se')])
     data.IL.a = data.frame(Participants=c(paste(power.IL.a$N.dyad,sep=";")),
     Power=power.IL.a$a,se=power.IL.a$a.se)
     Power.a = ggplot(data = data.IL.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed actor effect", x = "Number of dyads")

     power.IL.a.2 = data.frame(power.sim()$power[,c('N dyad','a.2','a.2.se')])
     data.IL.a.2 = data.frame(Participants=c(paste(power.IL.a.2$N.dyad,sep=";")),
     Power=power.IL.a.2$a.2,se=power.IL.a.2$a.2.se)
     Power.a.2 = ggplot(data = data.IL.a.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.a.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed actor effect", x = "Number of dyads")

     power.IL.p = data.frame(power.sim()$power[,c('N dyad','p','p.se')])
     data.IL.p = data.frame(Participants=c(paste(power.IL.p$N.dyad,sep=";")),
     Power=power.IL.p$p,se=power.IL.p$p.se)
     Power.p = ggplot(data = data.IL.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed partner effect", x = "Number of dyads")

     power.IL.p.2 = data.frame(power.sim()$power[,c('N dyad','p.2','p.2.se')])
     data.IL.p.2 = data.frame(Participants=c(paste(power.IL.p.2$N.dyad,sep=";")),
     Power=power.IL.p.2$p.2,se=power.IL.p.2$p.2.se)
     Power.p.2 = ggplot(data = data.IL.p.2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.p.2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed partner effect", x = "Number of dyads")

     power.IL.d = data.frame(power.sim()$power[,c('N dyad','d','d.se')])
     data.IL.d = data.frame(Participants=c(paste(power.IL.d$N.dyad,sep=";")),
     Power=power.IL.d$d,se=power.IL.d$d.se)
     Power.d = ggplot(data = data.IL.d, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the fixed effect of continuous predictor D", x = "Number of dyads")

     power.IL.d.a = data.frame(power.sim()$power[,c('N dyad','d.a','d.a.se')])
     data.IL.d.a = data.frame(Participants=c(paste(power.IL.d.a$N.dyad,sep=";")),
     Power=power.IL.d.a$d.a,se=power.IL.d.a$d.a.se)
     Power.d.a = ggplot(data = data.IL.d.a, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.a$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the actor effect", x = "Number of dyads")

     power.IL.d.a2 = data.frame(power.sim()$power[,c('N dyad','d.a2','d.a2.se')])
     data.IL.d.a2 = data.frame(Participants=c(paste(power.IL.d.a2$N.dyad,sep=";")),
     Power=power.IL.d.a2$d.a2,se=power.IL.d.a2$d.a2.se)
     Power.d.a2 = ggplot(data = data.IL.d.a2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.a2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the actor effect", x = "Number of dyads")

     power.IL.d.p = data.frame(power.sim()$power[,c('N dyad','d.p','d.p.se')])
     data.IL.d.p = data.frame(Participants=c(paste(power.IL.d.p$N.dyad,sep=";")),
     Power=power.IL.d.p$d.p,se=power.IL.d.p$d.p.se)
     Power.d.p = ggplot(data = data.IL.d.p, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.p$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the linear fixed moderation effect of dichotomous predictor D on the partner effect", x = "Number of dyads")

     power.IL.d.p2 = data.frame(power.sim()$power[,c('N dyad','d.p2','d.p2.se')])
     data.IL.d.p2 = data.frame(Participants=c(paste(power.IL.d.p2$N.dyad,sep=";")),
     Power=power.IL.d.p2$d.p2,se=power.IL.d.p2$d.p2.se)
     Power.d.p2 = ggplot(data = data.IL.d.p2, aes(x=Participants, y=Power, group=1)) + 
     geom_line(size=1) + geom_point() + scale_x_discrete(limits=data.IL.d.p2$Participants) +
     geom_errorbar(aes(ymin=Power-se, ymax=Power+se), width=.2,
                 position=position_dodge(.9)) +
     labs(title = "Power curve for the quadratic fixed moderation effect of dichotomous predictor D on the partner effect", x = "Number of dyads")

     coef.plot = grid.arrange(Power.c,Power.rho.Y,Power.a,Power.a.2,Power.p,Power.p.2,
     Power.d,Power.d.a,Power.d.a2,Power.d.p,Power.d.p2,ncol=1,heights=c(4,4,4,4,4,4,4,4,4,4,4))
     }}
   })
   

  output$replicates = renderText({power.sim()$`replicates`})

  output$power = renderFormattable({formattable(power.sim()$`coef.sim`)})

  output$covariance = renderFormattable({formattable(power.sim()$`cov.sim`)})

})

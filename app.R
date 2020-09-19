#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

banner <- "Leslie matrix model for marine mammal populations"

# Define UI for application that draws a histogram
ui <- fluidPage(theme="styles.css",

    # Application title
    titlePanel(banner),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            h4("Leslie matrix parameters"),
            fluidRow(
                column(4,"First year mortality /y"),
                column(8,numericInput("mu1",NULL,0.5,min=0,max=1,step=0.01))
            ),
            fluidRow(
                column(4,"Older mortality /y"),
                column(8,numericInput("mu2",NULL,0.1,min=0,max=1,step=0.01))
            ),
            fluidRow(
                column(4,"Maximum age (y)"),
                column(8,numericInput("n",NULL,20,min=2,max=40,step=1))
            ),
            fluidRow(
                column(4,"Age at first breeding (y)"),
                column(8,numericInput("b1",NULL,2,min=0,max=40,step=1))
            ),
            fluidRow(
                column(4,"Fecundity /y"),
                column(8,numericInput("fec",NULL,0.5,min=0,max=2,step=0.1))
            ),
            actionButton("updateLeslie","Update"),
            h4("Population initialisation"),
            fluidRow(
                column(4,"First year population size"),
                column(8,numericInput("x1",NULL,1,min=0,step=1))
            ),
            fluidRow(
                column(4,"Later years population size"),
                column(8,numericInput("x2",NULL,1,min=0,step=1))
            ),
            fluidRow(
                column(4,"Current year"),
                column(8,numericInput("year",NULL,2020,min=1900,max=30000,step=10))
            ),
            actionButton("updatePopulation","Update"),
            h4("Model operation"),
            fluidRow(
                column(4,"Number of years"),
                column(8,numericInput("years",NULL,1,min=1,step=100))
            ),
            actionButton("yearStep","Run model!")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Overview",
                         h4(banner),
                         p("The dominant eigenvalue of the Leslie matrix is ",
                           textOutput("ev",inline=T),
                           "."),
                         p("The total population size is ",
                            textOutput("N",inline=T),
                           "."),
                         p("The current year in the simulation is ",
                           textOutput("year",inline=T),
                           "."),
                         p("The previous year's R ratio was ",
                           textOutput("r",inline=T),
                           "."),
                         img(src='ioa_logo.png',style="width: 256px; align: left; margin-right: 2em"),
                         "Darren Green (2020)",
                         img(src='parasite_2.png',style="width: 64px; align: right; margin-left: 2em")
                ),
                tabPanel("Leslie matrix",tableOutput('leslie')),
                tabPanel("Leslie matrix heatmap",plotOutput('plot_leslie')),
                tabPanel("Population",tableOutput('population')), 
                tabPanel("Population graph",plotOutput('plot_population')),
                tabPanel("Time series",tableOutput("time_series")),
                tabPanel("Time series graph",plotOutput("plot_time_series"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    D <- reactiveValues()
    n <- 20
    D$n <- n
    D$Leslie <- matrix(0,n,n,dimnames=list(1:n,1:n))
    D$Population  <- array(0,n)
    D$Year <- 2020
    D$TimeSeries <- rbind(c(2020,0))

    observeEvent(input$updateLeslie,{
        D$n <- input$n
        D$Leslie <- matrix(0,D$n,D$n,dimnames=list(1:D$n,1:D$n))
        D$Leslie[1,2] <- 1-input$mu1
        for (i in 2:(D$n-1)) D$Leslie[i,i+1] <- 1-input$mu2
        if (input$b1<input$n)
        for (i in input$b1:input$n) D$Leslie[i,1] <- input$fec
        length(D$Population) <- D$n
        D$Population[is.na(D$Population)] <- 0
    })
    
    observeEvent(input$updatePopulation,{
        D$n <- input$n
        D$Population <- array(0,D$n)
        D$Population[1] <- input$x1
        for (i in 2:D$n) D$Population[i] <- input$x2
        D$Year <- input$year
        D$TimeSeries <- rbind(c(D$Year,sum(D$Population)))
    })
    
    observeEvent(input$yearStep,{
        for (i in 1:input$years) {
           D$Population <- D$Population %*% D$Leslie
            D$Year <- D$Year+1
            D$TimeSeries <- rbind(D$TimeSeries,c(D$Year,sum(D$Population)))
        }
    })
    
    output$leslie <- renderTable(t(D$Leslie))
    
    output$population <- renderTable({
        T <- cbind(1:D$n,as.vector(D$Population))
        colnames(T) <- c("age class","population size")
        return(T)
    })
    
    output$time_series <- renderTable({
        T <- D$TimeSeries
        colnames(T) <- c("year","population size")
        return(T)
    })
    
    output$plot_population <- renderPlot({
        par(cex.axis=1.5, cex.main=1.5, cex.lab=1.5)
        barplot(D$Population,xlab="age class",ylab="population size",col="#9d1e65")
    })
    
    output$ev <- renderText(Re(eigen(D$Leslie)$values[1]))
    
    output$N <- renderText(sum(D$Population))
    
    output$year <- renderText(sum(D$Year))
    
    output$plot_time_series <- renderPlot({
        par(cex.axis=1.5, cex.main=1.5, cex.lab=1.5)
        plot(D$TimeSeries[,1],D$TimeSeries[,2],type="o",xlab="year",ylab="population size",col="#9d1e65",
             xlim=c(min(D$TimeSeries[,1]-1),max(D$TimeSeries[,1]+1)),
             ylim=c(0,max(D$TimeSeries[,2]))
             )
    })
    
    output$plot_leslie <- renderPlot({
        image(D$Leslie[,D$n:1], axes=F)
    })
    
    output$r <- renderText({
        d <- dim(D$TimeSeries)
        if (d[1]<2) return("undefined")
        x2 <- D$TimeSeries[d[1],2]
        x1 <- D$TimeSeries[d[1]-1,2]
        return(x2/x1)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

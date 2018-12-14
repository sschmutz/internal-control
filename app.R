library(shiny)
library(tidyverse)
library(ggforce)
library(tibble)
library(patchwork)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Internal control"),
    
    sidebarLayout(
        sidebarPanel(width = 2,
               sliderInput("prob_background", "background", 
                           min = 0, max = 1, value = 0.9, 
                           step = 0.05, round = -2),
               br(),
               sliderInput("prob_internal_control", "internal control", 
                           min = 0, max = 1, value = 0.1, 
                           step = 0.05, round = -2),
               br(),
               hr(),
               br(),
               sliderInput("sequences_n", "sequencing reads", 
                           min = 0, max = 500, value = 80, 
                           step = 10, round = 0)
        ),
        mainPanel(
            plotOutput("plot", width = "100%", height = "500px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    

    output$plot <- renderPlot({
        
        set.seed(seed = 42)
        
        sequences_n <- input$sequences_n
        prob_background <- input$prob_background
        prob_internal_control <- input$prob_internal_control
        
        # define parameters of circle
        radius_circle <- 80
        center_x <- 100
        center_y <- 100
        
        # sample coordinates and amount of background/internal_control sequences
        x_coordinate <-
            sample(1:200, sequences_n, replace = TRUE)
        
        y_coordinate <-
            sample(1:200, sequences_n, replace = TRUE)
        
        type <-
            sample(c("background", "internal_control"), sequences_n, replace = TRUE, prob = c(prob_background, prob_internal_control))
        
        # put everything in a tibble
        coordinates <-
            tibble(x_coordinate, y_coordinate, type) %>%
            mutate(within_circle = if_else((x_coordinate - center_x)^2 + (y_coordinate - center_y)^2 <= radius_circle^2, TRUE, FALSE))
        
        plot_sequencing_pool <-
            coordinates %>%
            ggplot(aes(x = x_coordinate, y = y_coordinate)) +
            geom_point(aes(colour = type, alpha = within_circle), size = 5) +
            scale_colour_manual(values = c("#999999", "#56B4E9")) +
            scale_alpha_manual(values = c(0.2, 0.8)) +
            geom_circle(aes(x0 = center_x, y0 = center_y, r = radius_circle), inherit.aes = FALSE, colour = "grey") +
            labs(colour = "", caption = "points within ring represent sequenced DNA") +
            guides(alpha = FALSE) +
            theme_void() +
            theme(legend.position = "left", text = element_text(size = 20))
        
        sequencing_count <-
            coordinates %>%
            filter(within_circle == TRUE) %>%
            ggplot(aes(x = type)) +
            geom_bar(aes(fill = type)) +
            scale_fill_manual(values = c("#999999", "#56B4E9")) +
            labs(x = "", y = "sequence count \n") +
            guides(fill = FALSE) +
            theme_minimal() +
            theme(text = element_text(size = 20))
        
        plot_sequencing_pool + sequencing_count + plot_layout(ncol = 2)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

library(shiny)
library(tidyverse)
library(ggforce)
library(tibble)
library(patchwork)
library(shinyWidgets)

# Define UI for application that draws a histogram
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h2("Metagenomic sequencing"),
      h3("internal control"),
      p("Select the amount of",
        strong(span("background", style = "color:#999999")),
        "or",
        strong(span("internal control", style = "color:#56B4E9")),
        "points (DNA fragments) that should be displayed."),
      p("See how the different ratios affect the internal control reads per 1000 reads (rpk)."),
      em("Note: the sequence numbers don't have a specific meaning, it's also rpk instead of rpm."),
      width = 2,
      p(),
      chooseSliderSkin("Flat"),
      setSliderColor(c("#999999", "#56B4E9"), c(1, 2)),
      sliderInput("background_n", "background",
                  min = 10, max = 200, value = 80,
                  step = 10, round = 0),
      br(),
      sliderInput("internal_control_n", "internal control",
                  min = 5, max = 50, value = 10,
                  step = 5, round = 0)
      ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot", width = "100%", height = "1000px")),
                  tabPanel("Documentation",
                           p(),
                           "Using an internal control (IC) for metagenomic sequencing
                           could lead to a better monitoring of wet lab performance.",
                           br(),
                           "Goal of this app is to show how the sequenced
                           reads of the IC are dependant not only
                           of the amount spiked in but also the amount of background
                           DNA.",
                           h4("How does it work"),
                           "The amount of background and internal control DNA fragments
                           can be chosen in the sidebar panel.",
                           br(),
                           "Those numbers represent the total points per category
                           in the plot. Of those only a (random) subset are
                           sequenced (points within ring).",
                           br(),
                           "The bar plot shows how many DNA fragments of each
                           category are sequenced.")
                  )
      )
    )
  )

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$plot <- renderPlot({
    
    set.seed(seed = 42)
    background_n <- input$background_n
    internal_control_n <- input$internal_control_n
    sequences_n <- background_n + internal_control_n
    
    # define parameters of circle
    radius_circle <- 80
    center_x <- 100
    center_y <- 100
    
    # sample coordinates and amount of background/internal_control sequences
    x_coordinate <-
      runif(n = sequences_n, min = 1, max = 200)
    
    y_coordinate <-
      runif(n = sequences_n, min = 1, max = 200)
    
    internal_control_row <-
      sample(1:sequences_n, internal_control_n, replace = FALSE)
    
    # put everything in a tibble
    coordinates <-
      tibble(x_coordinate, y_coordinate) %>%
      mutate(type = if_else(row_number() %in% internal_control_row, "internal_control", "background")) %>%
      mutate(within_circle = if_else((x_coordinate - center_x)^2 + (y_coordinate - center_y)^2 <= radius_circle^2, TRUE, FALSE))
    
    plot_sequencing_pool <-
      coordinates %>%
      ggplot(aes(x = x_coordinate, y = y_coordinate)) +
      geom_point(aes(colour = type, alpha = within_circle), size = 5) +
      scale_colour_manual(values = c("#999999", "#56B4E9")) +
      scale_alpha_manual(values = c(0.2, 0.8)) +
      geom_circle(aes(x0 = center_x, y0 = center_y, r = radius_circle), inherit.aes = FALSE, colour = "grey") +
      labs(colour = "", caption = "points within ring represent sequenced DNA") +
      guides(colour = FALSE, alpha = FALSE) +
      theme_void() +
      theme(legend.position = "left", text = element_text(size = 20))
    
    internal_control_circle <-
      coordinates %>%
      filter(within_circle == TRUE) %>%
      filter(type == "internal_control") %>%
      nrow()
    
    total_circle <-
      coordinates %>%
      filter(within_circle == TRUE) %>%
      nrow()
    
    sequencing_count <-
      coordinates %>%
      filter(within_circle == TRUE) %>%
      ggplot(aes(x = type)) +
      geom_bar(aes(fill = type)) +
      scale_fill_manual(values = c("#999999", "#56B4E9")) +
      scale_y_continuous(limits = c(0, 100)) +
      scale_x_discrete(labels = c("background", "internal control")) +
      labs(x = "", y = "sequence count \n", caption = paste("rpk =", round(internal_control_circle / total_circle * 1000))) +
      guides(fill = FALSE) +
      theme_minimal() +
      theme(text = element_text(size = 20))
    
    plot_sequencing_pool + sequencing_count + plot_layout(ncol = 1)
  })
  }

# Run the application 
shinyApp(ui = ui, server = server)

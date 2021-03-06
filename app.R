require(shiny); require(shinyBS)
require(tidyverse); require(magrittr)
require(rlist); require(data.tree); require(collapsibleTree)
require(brapir)
source("utils.R")
parents <<- read.csv("parents.csv", stringsAsFactors = F) %>%
  rowwise() %>%
  mutate(Code = toGenericPedigreeCode(Code))


ui <- navbarPage(id = "tabs", collapsible = TRUE, title = "Pedigree Toolbox",
    tabPanel("Pedigree Viewer", 
      fluidPage(
        textInput("germplasmName", "Germplasm Name:"),
        bsCollapsePanel(title = "Advanced Options",
          checkboxInput("collapsed", "Show full pedigree by default?")
        ),
        tableOutput("percentParents"),
        collapsibleTreeOutput("pedigreeTree")
      )
    ),
    tags$head(tags$link(rel="shortcut icon", href="https://thumbs.dreamstime.com/b/sweet-potato-white-background-sweet-potato-batata-white-background-isolated-103677860.jpg"))
  )
  
server <- function(input, output, session) {
  
  output$percentParents <- renderTable({
    validate(need((input$germplasmName != ""), ""))
    validate(need(!isTerminal(input$germplasmName), "Terminal parent or no data"))
    withProgress(message = "Getting terminal parents",{
      input$germplasmName %>% getPercentParents() %>%
      scales::label_percent()(.) %>% as.list %>% data.frame(check.names = F)
    })
  })
    
  
  output$pedigreeTree <- renderCollapsibleTree({
    validate(need((input$germplasmName != ""), ""))
    withProgress(message = "Building pedigree tree",{
      validate(need(!isTerminal(input$germplasmName), ""))
      collapsibleTree(fontSize = 20, linkLength = "200", collapsed = !input$collapsed,
        FromListSimple(pedigree(input$germplasmName)[[1]], nodeName = input$germplasmName)
      )
    })
  })
  
}

shinyApp(ui = ui, server = server)
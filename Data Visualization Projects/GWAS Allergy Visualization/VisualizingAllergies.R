library(shiny)
library(DT)
library(readr)
library(ggplot2)

# Read data
GeneralAllergies <- read_tsv("Asthma and Dermatitis.tsv")
DrugAllergies <- read_tsv("Drug Allergies.tsv")
FoodAllergies <- read_tsv("Food Allergies.tsv")


#Change riskFrequency and orValues to numerical

if (!is.numeric(FoodAllergies$riskFrequency)) {
  FoodAllergies$riskFrequency <- as.numeric(as.character(FoodAllergies$riskFrequency))
}

if (!is.numeric(GeneralAllergies$riskFrequency)) {
  #GeneralAllergies$riskFrequency <- ifelse(GeneralAllergies$riskFrequency == 'NR', -1, as.numeric(as.character(GeneralAllergies$riskFrequency)))
  GeneralAllergies$riskFrequency <- as.numeric(as.character(GeneralAllergies$riskFrequency))
}

if (!is.numeric(DrugAllergies$riskFrequency)) {
  DrugAllergies$riskFrequency <- as.numeric(as.character(DrugAllergies$riskFrequency))
}



if (!is.numeric(FoodAllergies$orValue)) {
 # FoodAllergies$orValue <- as.numeric(ifelse(FoodAllergies$orValue == "-" | is.na(as.numeric(FoodAllergies$orValue)), NA, FoodAllergies$orValue))
  FoodAllergies$orValue <- as.numeric(as.character(FoodAllergies$orValue))
}

if (!is.numeric(GeneralAllergies$orValue)) {
  #GeneralAllergies$orValue <- as.numeric(ifelse(GeneralAllergies$orValue == "-" | is.na(as.numeric(GeneralAllergies$orValue)), NA, GeneralAllergies$orValue))
  GeneralAllergies$orValue <- as.numeric(as.character(GeneralAllergies$orValue))
}

if (!is.numeric(DrugAllergies$orValue)) {
  #DrugAllergies$orValue <- as.numeric(ifelse(DrugAllergies$orValue == "-" | is.na(as.numeric(DrugAllergies$orValue)), NA, DrugAllergies$orValue))
  DrugAllergies$orValue <- as.numeric(as.character(DrugAllergies$orValue))
  
}

DrugAllergies[DrugAllergies=='-']<-NA

GeneralAllergies[GeneralAllergies=='-']<-NA

FoodAllergies[FoodAllergies=='-']<-NA

# Define UI
ui <- fluidPage(
  titlePanel("Allergies Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("General Allergies", 
                 selectInput("general_color", "Bar Color For General Allergies:",
                             choices = c("skyblue", "red", "green", "blue"))),
        tabPanel("Drug Allergies", 
                 selectInput("drug_color", "Bar Color For Drug Allergies:",
                             choices = c("skyblue", "red", "green", "blue"))),
        tabPanel("Food Allergies", 
                 selectInput("food_color", "Bar Color For Food Allergies:",
                             choices = c("skyblue", "red", "green", "blue")))
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Project Description", 
                 # Content for the text tab
                 textOutput("ProjectDescription")
        ),
        tabPanel("General Allergies", 
                 tags$div(style = "width: 1000px; overflow-x: scroll;",plotOutput("general_gene_plot")),
                 tags$div(style = "width: 1000px; overflow-x: scroll;",plotOutput("general_allele_plot")),
                 tags$div(style = "width: 1000px; overflow-x: scroll;",plotOutput("general_allergen_plot")),
                 tags$div(style = "width: 1000px; overflow-x: scroll;",plotOutput("general_gene_bubbleplot")),
                 dataTableOutput("general_table")),
        tabPanel("Drug Allergies", 
                 plotOutput("drug_gene_plot"),
                 plotOutput("drug_allele_plot"),
                 plotOutput("drug_allergen_plot"),
                 plotOutput("drug_gene_bubbleplot"),
                 dataTableOutput("drug_table")),
        tabPanel("Food Allergies", 
                 plotOutput("food_gene_plot"),
                 plotOutput("food_allele_plot"),
                 plotOutput("food_allergen_plot"),
                 plotOutput("food_gene_bubbleplot"),
                 dataTableOutput("food_table"))
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  output$ProjectDescription <- renderText({
    
  "This project performs EDA on a dataset that contains data regarding \n
  the relationship between Genes and allergies. The data was obtained from GWAS. \n
  The Allergies studied include General Allergies, Food Allergies and Drug Allergies. \n
  The EDA in each tab for all 3 allergies consists of the following: \n
  1.- Barcharts showing the top 5 genes, alleles and allergens with the highest frequency in each dataset. \n
  2.- A bubble chart showing how alleles are concentrated in relation to allergies, thereby giving an understanding of which allergies are more common."
  
  })


  output$general_gene_plot <- renderPlot({
    gene_freq <- table(GeneralAllergies$mappedGenes)
    top_5_genes <- names(sort(gene_freq, decreasing = TRUE))[1:5]
    barplot(gene_freq[top_5_genes], main = "Top 5 Gene Frequency - General Allergies",
            xlab = "Genes", ylab = "Frequency", col = input$general_color)
  })
  
  
  output$general_allele_plot <- renderPlot({
    allele_freq <- table(GeneralAllergies$riskAllele)
    top_5_alleles <- names(sort(allele_freq, decreasing = TRUE))[1:5]
    barplot(allele_freq[top_5_alleles], main = "Top 5 Allele Frequency - General Allergies",
            xlab = "Alleles", ylab = "Frequency", col = input$general_color)
  })
  
  output$general_allergen_plot <- renderPlot({
    allergen_freq <- table(GeneralAllergies$traitName)
    top_5_allergens <- names(sort(allergen_freq, decreasing = TRUE))[1:5]
    barplot(allergen_freq[top_5_allergens], main = "Top 5 Allergen Frequency - General Allergies",
            xlab = "Allergens", ylab = "Frequency", col = input$general_color)
  })

  output$general_gene_bubbleplot <- renderPlot({
    ggplot(GeneralAllergies, aes(x = riskAllele, y = traitName, size = orValue)) +
      geom_point() +
      scale_size_continuous(name = "Odds Ratio") +
      scale_color_manual(values = input$general_color) +
      labs(x = "Risk Allele", y = "Trait Name", 
           title = "Allele Concentration In General Allergy Types",
           col = input$general_color)
    
  })
  
  output$general_table <- renderDataTable({
    datatable(GeneralAllergies)
  })
  
  output$drug_gene_plot <- renderPlot({
    gene_freq <- table(DrugAllergies$mappedGenes)
    top_5_genes <- names(sort(gene_freq, decreasing = TRUE))[1:5]
    barplot(gene_freq[top_5_genes], main = "Top 5 Gene Frequency - Drug Allergies",
            xlab = "Genes", ylab = "Frequency", col = input$drug_color)
  })
  
  output$drug_allele_plot <- renderPlot({
    allele_freq <- table(DrugAllergies$riskAllele)
    top_5_alleles <- names(sort(allele_freq, decreasing = TRUE))[1:5]
    barplot(allele_freq[top_5_alleles], main = "Top 5 Allele Frequency - Drug Allergies",
            xlab = "Alleles", ylab = "Frequency", col = input$drug_color)
  })
  
  output$drug_allergen_plot <- renderPlot({
    allergen_freq <- table(DrugAllergies$traitName)
    top_5_allergens <- names(sort(allergen_freq, decreasing = TRUE))[1:5]
    barplot(allergen_freq[top_5_allergens], main = "Top 5 Allergen Frequency - Drug Allergies",
            xlab = "Allergens", ylab = "Frequency", col = input$drug_color)
  })
  
  output$drug_gene_bubbleplot <- renderPlot({
    ggplot(DrugAllergies, aes(x = riskAllele, y = traitName, size = orValue)) +
      geom_point() +
      scale_size_continuous(name = "Odds Ratio") +
      labs(x = "Risk Allele", y = "Trait Name", 
           title = "Allele Concentration In Drug Allergy Types",
           col = input$drug_color)
    
  })
  
  output$drug_table <- renderDataTable({
    datatable(DrugAllergies)
  })
  
  output$food_gene_plot <- renderPlot({
    gene_freq <- table(FoodAllergies$mappedGenes)
    top_5_genes <- names(sort(gene_freq, decreasing = TRUE))[1:5]
    barplot(gene_freq[top_5_genes], main = "Top 5 Gene Frequency - Food Allergies",
            xlab = "Genes", ylab = "Frequency", col = input$food_color)
  })
  
  output$food_allele_plot <- renderPlot({
    allele_freq <- table(FoodAllergies$riskAllele)
    top_5_alleles <- names(sort(allele_freq, decreasing = TRUE))[1:5]
    barplot(allele_freq[top_5_alleles], main = "Top 5 Allele Frequency - Food Allergies",
            xlab = "Alleles", ylab = "Frequency", col = input$food_color)
  })
  
  output$food_allergen_plot <- renderPlot({
    allergen_freq <- table(FoodAllergies$traitName)
    top_5_allergens <- names(sort(allergen_freq, decreasing = TRUE))[1:5]
    barplot(allergen_freq[top_5_allergens], main = "Top 5 Allergen Frequency - Food Allergies",
            xlab = "Allergens", ylab = "Frequency", col = input$food_color)
  })
  
  
  output$food_gene_bubbleplot <- renderPlot({
    ggplot(FoodAllergies, aes(x = riskAllele, y = traitName, size = orValue)) +
      geom_point() +
      scale_size_continuous(name = "Odds Ratio") +
      labs(x = "Risk Allele", y = "Trait Name", title = "Allele Concentration In Food Allergy Types",
           col = input$food_color)
    
  })
  
  output$food_table <- renderDataTable({
    datatable(FoodAllergies)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

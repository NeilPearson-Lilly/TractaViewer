
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(plotly)

shinyUI(fluidPage(

  # Application title
  titlePanel("TractaViewer mined data visualiser"),

  fixedRow(
    column(3, selectizeInput('gene_selection', NULL, width = '100%', choices = NULL, options = list(maxItems = 1))),
    column(1, submitButton("Go", icon("refresh"))), 
    column(5, h4(textOutput("ensembl_id")),
              h5(textOutput("full_gene_name")),
              h5(textOutput("basic_function"))),
    column(3, h5(textOutput("protein_class")),
              h5(textOutput("basic_location")))
  ),
  
  fluidRow(
    column(4, 
           h2("Small-molecule druggability"),
           # h4("Known interactions"),
           # tableOutput("interactions"),
           h4("Recommendations"),
           tableOutput("modality_recommendations"),
           h4("Bucketing overview"),
           plotOutput('buckets_spiderplot')
           ),
    column(5, 
           h2("Safety"),
           h4("Withdrawn drugs"),
           tableOutput("withdrawn_drug_target_table"),
           h4("Tissue-wise differential expression (Human Protein Atlas)"),
           # plotlyOutput("tissue_expression_table", height = 210),
           plotlyOutput("tissue_expression_table", height="auto"),
           h4("Brain tissue-wise differential expression (Barres mouse)"),
           # plotlyOutput("barres_tissue_expression_table", height = 270, width = 450),
           plotlyOutput("barres_tissue_expression_table", height="auto"),
           h4("Associations with genetic disorders"),
           tableOutput("genetic_disorders_table"),
           h4("Tissues affected by associated disorders"),
           # plotlyOutput("hpo_phen_table", height=170),
           plotlyOutput("hpo_phen_table", height="auto"),
           h4("Gene essentiality"),
           tableOutput("essentiality_table"),
           h4("Phenotypic trait associations in mouse"),
           # plotlyOutput("mouse_phen_table", height = 240)
           plotlyOutput("mouse_phen_table", height="auto")
           ),
    column(3, 
           h2("Feasibility"),
           tableOutput("pharos_table"),
           tableOutput("tool_compounds_table"),
           tableOutput("orthology_table"),
           tableOutput("assays"),
           tableOutput("literature"),
           tableOutput("protein_structures"),
           tableOutput("gene_families"),
           plotOutput("antitargets", height = 250)
           )
  )
))

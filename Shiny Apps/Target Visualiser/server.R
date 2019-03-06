
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(tidyverse)
library(fmsb)
library(plotly)
library(DBI)
library(RSQLite)
library(data.table)
library(stringr)
library(stringi)
library(DT)

con = dbConnect(SQLite(), dbname="proteome.db")
barreslab_enrichments = read.csv("BarresLab_Enrichments.csv")

# interactions = read.csv("IUPHAR_ligand-gene_interactions.csv")
# interactions = interactions %>% filter(target_species == "Human")

# This gives us exactly the thing we'll usually want to display
flip_df_round <- function(x) {
  tdata = rownames_to_column(as.data.frame(t(x)))
  colnames(tdata) = c("Property", "Value")
  tdata
}


shinyServer(function(input, output, session) {
  allgenes = sort(unname(unlist(dbFetch(dbSendQuery(con, 'SELECT DISTINCT `HGNC.Name` FROM `Basic information`')))))
  updateSelectizeInput(session, 'gene_selection', 
                       choices = allgenes,
                       server = TRUE)

  output$basic_table <- renderTable({
    data = xl$`Basic information` %>% filter(HGNC.Name == input$gene_selection)
    subdata = data[,c("HGNC.Name", "GeneID", "Uniprot.ID", "Entrez.ID", "Mouse.Ensembl.ID", "Mouse.Uniprot.ID", "Is.protein", "RNA.class")]
    flip_df_round(subdata)
  })
  
  output$feasibility_table <- renderTable({
    data = xl$Feasibility %>% filter(HGNC.Name == input$gene_selection)
    subdata = subdata = data[5:length(colnames(data))]
    flip_df_round(subdata)
  })
  
  output$subloc_table <- renderTable({
    data = xl$`SM Druggability` %>% filter(HGNC.Name == input$gene_selection)
    subdata = data[,c("Main.location")]
    flip_df_round(subdata)
  })
  
  # output$interactions <- renderTable({
  #   req(input$gene_selection)
  #   data = interactions %>% filter(target_gene_symbol == input$gene_selection)
  #   validate(
  #     need(nrow(data) != 0, "No interactions on record.")
  #   )
  #   unique(data[,c("ligand", "type", "action", "selectivity")])
  # })
  
  output$druggability_table <- renderTable({
    data = xl$`SM Druggability` %>% filter(HGNC.Name == input$gene_selection)
    subdata = data[5:length(colnames(data))]
    flip_df_round(subdata)
  })
  
  output$omim_table <- renderTable({
    data = xl$`Disease associations` %>% filter(HGNC.Name == input$gene_selection)
    subdata = data[6:length(colnames(data))]
  })
  
  output$drug_table <- renderTable({
    data = xl$`Existing drugs` %>% filter(HGNC.Name == input$gene_selection)
    subdata = data[4:length(colnames(data))]
  })
  
  output$buckets_spiderplot <- renderPlot({
    # https://www.r-graph-gallery.com/142-basic-radar-chart/
    req(input$gene_selection)
    buckets = dbFetch(dbSendQuery(con, 'SELECT * FROM Buckets 
                                  WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    buckets = buckets[grep("bucket", colnames(buckets))]
    # I think we need 3 rows: max, min, and actual value.
    max_column_values = c(6, 14, 8, 7, 4, 5)
    max_column_values = max_column_values - 1
    subdata = rbind(
      # rep(9,length(colnames(buckets))), 
      # Column names, in order: Safety, SM.Druggability, Feasibility, AB.ability, Modality
      max_column_values,
      rep(0,length(colnames(buckets))) , 
      max_column_values - (buckets[1,] - 1))
    radarchart(subdata, axistype=1)
  })
  
  output$ligandability <- renderTable({
    req(input$gene_selection)
    bucket_name = c(
      "Bucket 1: Druggable, by precedent",
      "Bucket 2: Targetable, by homology",
      "Bucket 3: Targetable, structurally enabled",
      "Bucket 4: Targetable, by homology",
      "Bucket 5: Probably targetable",
      "Bucket 6: Probably targetable, by homology", 
      "Bucket 7: Potentially targetable by family structure", 
      "Bucket 8: Endogenous ligand", 
      "Bucket 9: Potentially targetable by family ligand",
      "Bucket 10: Potentially targetable by low activity family ligand", 
      "Bucket 11: Targetable, by class",
      "Bucket 12: Low SM druggability",
      "Bucket 13: Unknown druggability",
      "Bucket 14: Non-protein target"
    )
    bucket_info = c(
      "Target has known small molecule ligand passing activity threshold.", 
      "Targetable by homology – homolog has known small molecule ligand passing activity threshold.",
      "Target has structural indications of druggability.",
      
      "Targetable by homology – homolog has structural indications of druggability.",
      "Protein with a SM ligand identified from ChEMBL data, but the ligand does not meeting TCRD activity criteria.",
      "Targetable by homology – homolog has known small molecule ligand below activity threshold.",
      
      "Targetable by gene family – family member has structural indications of druggability.",
      "Target has endogenous ligand but no known small molecule ligands.",
      "Member of a gene family which has a member with an SM ligand meeting activity criteria.",
      
      "Targetable by homology – homolog has known small molecule ligand, but with low activity.",
      "Potentially druggable protein class; no other information available.",
      "High-resolution 3D structure available but no evidence of a druggable pocket.",
      
      "No information available.",
      "Non-protein target."
    )
    buckets = unlist(dbFetch(dbSendQuery(con, 'SELECT * FROM Buckets 
                                          WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    df = as.data.frame(c(bucket_info[as.integer(buckets['SM.Druggability.bucket'])]))
    colnames(df) = c(bucket_name[as.integer(buckets['SM.Druggability.bucket'])])
    df
  })
  
  output$antibodyability <- renderTable({
    req(input$gene_selection)
    bucket_name = c(
      "Bucket 1: Most accessible",
      "Bucket 2: Highly accessible",
      "Bucket 3: Highly accessible",
      "Bucket 4: Accessible",
      "Bucket 5: Less accessible",
      "Bucket 6: Less accessible",
      "Bucket 7: Unknown AB-bility"
    )
    bucket_info = c(
      "Secreted protein. Highly accessible to antibody-based therapies.", 
      "Component of the extracellular matrix (ECM). Highly accessible to antibody-based therapies, but potentially less so than secreted proteins.",
      "Cell membrane-bound proteins. Highly accessible to antibody-based therapies, but potentially less so than secreted proteins or ECM components.",
      
      "Limited evidence that target is a secreted protein, ECM component or cell membrane-bound protein. ",
      
      "Protein located in the cytosol. Not practically accessible to antibody-based therapies, but may be more easily accessible to other modalities.",
      "Protein located in intracellular compartment.",
      "Dark target. Paucity of biological knowledge means progress will be difficult. Encourage others to investigate. "
    )
    buckets = unlist(dbFetch(dbSendQuery(con, 'SELECT * FROM Buckets 
                                          WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    df = as.data.frame(c(bucket_info[as.integer(buckets['AB.ability.bucket'])]))
    colnames(df) = c(bucket_name[as.integer(buckets['AB.ability.bucket'])])
    df
  })
  
  output$new_modality <- renderTable({
    req(input$gene_selection)
    bucket_name = c(
      "Bucket 1: New modalities required",
      "Bucket 2: New modalities viable",
      "Bucket 3: New modalities viable but other options preferred",
      "Bucket 4: New modalities not requested"
    )
    bucket_info = c(
      "Target is a non-protein-coding RNA; or has been specified suitable for degradation/inhibition by user, and is unsuitable for both small-molecule and antibody modalities.", 
      "Specified suitable for degradation/inhibition by user, but one of either small-molecule or antibody modalities may be preferred.",
      "Specified suitable for degradation/inhibition by user, but both small-molecule and antibody modalities may be preferred.",
      "Not specified suitable for degradation/inhibition by user."
    )
    buckets = unlist(dbFetch(dbSendQuery(con, 'SELECT * FROM Buckets 
                                         WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    df = as.data.frame(c(bucket_info[as.integer(buckets['New.modality.bucket'])]))
    colnames(df) = c(bucket_name[as.integer(buckets['New.modality.bucket'])])
    df
  })
  
  output$modality_recommendations <- renderTable({
    req(input$gene_selection)
    ligandability_buckets = c(
      "Bucket 1: Druggable, by precedent",
      "Bucket 2: Targetable, by homology",
      "Bucket 3: Targetable, structurally enabled",
      "Bucket 4: Targetable, by homology",
      "Bucket 5: Probably targetable",
      "Bucket 6: Probably targetable, by homology", 
      "Bucket 7: Potentially targetable by family structure", 
      "Bucket 8: Endogenous ligand", 
      "Bucket 9: Potentially targetable by family ligand",
      "Bucket 10: Potentially targetable by low activity family ligand", 
      "Bucket 11: Targetable, by class",
      "Bucket 12: Low SM druggability",
      "Bucket 13: Unknown druggability",
      "Bucket 14: Non-protein target"
    )
    antibodyability_buckets = c(
      "Bucket 1: Most accessible",
      "Bucket 2: Highly accessible",
      "Bucket 3: Highly accessible",
      "Bucket 4: Accessible",
      "Bucket 5: Less accessible",
      "Bucket 6: Less accessible",
      "Bucket 7: Unknown AB-bility"
    )
    new_modality_buckets = c(
      "Bucket 1: New modalities required",
      "Bucket 2: New modalities viable",
      "Bucket 3: New modalities viable but other options preferred",
      "Bucket 4: New modalities not requested or unknown"
    )
    buckets = unlist(dbFetch(dbSendQuery(con, 'SELECT * FROM Buckets 
                                         WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    
    druggability = as.integer(buckets['SM.Druggability.bucket'])
    abability = as.integer(buckets['AB.ability.bucket'])
    new_mod = as.integer(buckets['New.modality.bucket'])
    
    # unitick = "☑"
    # unicross = "☒"
    
    df = t(data.frame(
      c("Ligandability", 
        ifelse(druggability < length(ligandability_buckets), 
               ifelse(druggability <= 11, 
                      ifelse(druggability <= 6, "Yes", "Yes (low priority)"), "No"), "Unknown"), 
        ligandability_buckets[druggability]),
      c("AB-ability", 
        ifelse(abability < length(antibodyability_buckets), ifelse(abability <= 4, "Yes", "No"), "Unknown"), 
        antibodyability_buckets[abability]),
      c("New modalities (ASOs siRNA etc.)", 
        ifelse(new_mod <= 3, "Yes", "Unknown"), 
        new_modality_buckets[new_mod])
    ))
    colnames(df) = c("Therapeutic class", "Suitable", "Status")
    df
  })
  
  output$full_gene_name <- renderText({
    req(input$gene_selection)
    fullname = unlist(dbFetch(dbSendQuery(con, 'SELECT `Approved.Name`, Synonyms FROM `Basic information` 
                                          WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    
    paste0(fullname[1], " (Synonyms: ", fullname[2], ")")
  })
  
  output$basic_function <- renderText({
    req(input$gene_selection)
    funcsum = unlist(dbFetch(dbSendQuery(con, 'SELECT `Functional.summary` FROM `Basic information` 
                                         WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    if (identical(funcsum, character(0))) { funcsum = "(No functional summary available)" } 
    if (is.na(funcsum)) { funcsum = "(No functional summary available)" } 
    funcsum
  })
  
  output$protein_class <- renderText({
    req(input$gene_selection)
    funcsum = unlist(dbFetch(dbSendQuery(con, 'SELECT `Top.level.protein.classes` FROM `Basic information` 
                                         WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    if (identical(funcsum, character(0))) { funcsum = "No subcellular location knowb" } 
    if (is.na(funcsum)) { funcsum = "No protein class known" } 
    paste("Protein Class:", gsub(";", "; ", funcsum))
  })
  
  output$basic_location <- renderText({
    req(input$gene_selection)
    funcsum = unlist(dbFetch(dbSendQuery(con, 'SELECT `Main.location` FROM `SM druggability` 
                                         WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    if (identical(funcsum, character(0))) { funcsum = "No subcellular location known" } 
    if (is.na(funcsum)) { funcsum = "No subcellular location known" } 
    paste("Human Protein Atlas main location:", gsub(";", "; ", funcsum))
  })
  
  output$ensembl_id <- renderText({
    req(input$gene_selection)
    ensembl = unlist(dbFetch(dbSendQuery(con, 'SELECT GeneID FROM `Basic information` 
                                         WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1))
    if (identical(ensembl, character(0))) { ensembl = "(No Ensembl ID available)" } 
    if (is.na(ensembl)) { ensembl = "(No Ensembl ID available)" } 
    ensembl
  })
  
  output$withdrawn_drug_target_table <- renderTable({
    req(input$gene_selection)
    withdrawn = dbFetch(dbSendQuery(con, 'SELECT * FROM `Existing drugs` 
                                    WHERE `HGNC.Name` == :x AND Withdrawn == 1', params = list(x = input$gene_selection)), n = 1)
    wd = data.frame(c(ifelse(nrow(withdrawn) > 0, "Yes", "No")))
    colnames(wd) = c("Intended target for withdrawn drug:")
    wd
  })
  
  output$tissue_expression_table <- renderPlotly({
    req(input$gene_selection)
    tex = dbFetch(dbSendQuery(con, 'SELECT * FROM `HPA Enrichment` 
                              WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    tex = tex[,5:ncol(tex)]
    
    # validate(
    #   need(!all(is.na(tex)), "No tissues present significant differential expression of this gene.")
    # )
    
    diffexp_classes = c("Tissue enriched", "Group enriched", "Tissue enhanced")
    df = as.data.frame(do.call(rbind, lapply(diffexp_classes, function (i) {
      lapply(tex, function (x) {as.integer(grepl(i, x))})
    })))
    rownames(df) = diffexp_classes
    df = as.data.frame(t(df)) %>% rownames_to_column("tissue")
    # Right - now need to melt that.
    df = melt(as.data.table(df), id.vars = "tissue")
    colnames(df) = c("Tissue", "Expression", "value")
    df$Tissue = gsub("HPA.", "", df$Tissue)
    df$value = as.integer(df$value)
    
    colpal = c(ifelse(min(df$value) == 0, "white", "blue"),
               ifelse(max(df$value) == 0, "white", "blue"))
    
    plot_ly(x=df$Tissue,
            y=df$Expression,
            z=df$value,
            type = "heatmap",
            height = 210,
            showscale=FALSE,
            colors = colorRamp(colpal)) %>%
      layout(margin = list(l = 150, r = 70, b = 30, t = 120, pad = 8),
             xaxis = list(side ="top", tickangle=-45))
  })
  
  output$barres_tissue_expression_table <- renderPlotly({
    # Here, it looks like I'll have to query for the mouse homolog gene name first, then filter the CSV file based on that. 
    req(input$gene_selection)
    tex = dbFetch(dbSendQuery(con, 'SELECT `MGI.Symbol`, `HCOP.Mouse` FROM `Basic information` 
                              WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    # Mouse name can be in two places. MGI should be preferred, but if not, go with HCOP Mouse.
    # Due to evolutionary processes, we may get more than one name come out here. We'll deal with that shortly.
    mousenames = NULL
    if (!is.na(tex$MGI.Symbol[1])) {
      mousenames = tex$MGI.Symbol[1]
    } else {
      mousenames = unlist(strsplit(tex$HCOP.Mouse[1], ","))
    }
    relevant_enrichments = barreslab_enrichments[barreslab_enrichments$Gene %in% mousenames,]
    
    # validate(
    #   need(nrow(relevant_enrichments) != 0, "No tissues present significant differential expression of this gene.")
    # )
    # That's what we need! Now slot it into the existing stuff somehow.
    # Pro tip: this is actually in the melted format we need already - but we need to supply 0/1 (true/false) values, and
    # we also need to fill in any combinations that aren't currently there. 
    # Easiest way would be to create this, tissues x groups, and apply in values accordingly. 
    df = data.frame(matrix(ncol = 3, nrow = 0))
    colnames(df) = c("Tissue", "Expression", "value")
    for (x in unique(barreslab_enrichments$Tissue)) {
      for (y in unique(barreslab_enrichments$Group)) {
        v = 0
        if (nrow(filter(relevant_enrichments, Tissue == x & Group == y)) > 0) {
          v = 1
        }
        df[nrow(df) + 1,] = list(x, y, v)
      }
    }
    
    colpal = c(ifelse(min(df$value) == 0, "white", "blue"),
               ifelse(max(df$value) == 0, "white", "blue"))
    
    plot_ly(x=df$Tissue,
            y=df$Expression,
            z=df$value,
            type = "heatmap",
            height = 270,
            width = 450,
            showscale=FALSE,
            colors = colorRamp(colpal)) %>%
      layout(margin = list(l = 150, r = 145, b = 30, t = 180, pad = 8),
             xaxis = list(side ="top", tickangle=-45))
  })
  
  output$pharos_table <- renderTable({
    req(input$gene_selection)
    pharos = dbFetch(dbSendQuery(con, 'SELECT Pharos FROM `Pharos` 
                                 WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    colnames(pharos) = c("Pharos categorisation")
    pharos
  }, digits=0, align='l')
  
  output$tool_compounds_table <- renderTable({
    req(input$gene_selection)
    pharos = dbFetch(dbSendQuery(con, 'SELECT `ChEMBL.ligand` FROM `Pharos` 
                                 WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    colnames(pharos) = c("Tool compounds meeting activity thresholds")
    pharos
  }, digits=0, align='l')
  
  output$assays <- renderTable({
    req(input$gene_selection)
    assays = dbFetch(dbSendQuery(con, 'SELECT Assays FROM `Feasibility` 
                                 WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
  }, digits=0, align='l')
  
  output$literature <- renderTable({
    req(input$gene_selection)
    literature = dbFetch(dbSendQuery(con, 'SELECT Literature FROM `Feasibility` 
                                     WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    colnames(literature) = c("Publications in literature")
    literature
  }, digits=0, align='l')
  
  output$protein_structures <- renderTable({
    req(input$gene_selection)
    structures = dbFetch(dbSendQuery(con, 'SELECT COUNT(*) FROM `Protein structures` 
                                     WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)))
    colnames(structures) = c("PDB protein structure models")
    structures
  }, digits=0, align='l')
  
  output$gene_families <- renderTable({
    req(input$gene_selection)
    families = dbFetch(dbSendQuery(con, 'SELECT `Gene.family`, `Gene.family.ID` FROM `Basic information` 
                                   WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)))
    families = families %>% separate_rows(colnames(families), sep="[|]")
    
    toLink <- function(gf, gf_id) {
      paste0('<a href="https://www.genenames.org/cgi-bin/genefamilies/set/', gf_id, '/downoad/node">', gf, '</a>')
    }
    
    df = as.data.frame(apply(families, 1, function (x) toLink(x['Gene.family'], x['Gene.family.ID'])))
    colnames(df) = c("Gene families")
    df
  }, align='l', sanitize.text.function = function(x) x)
  
  output$essentiality_table <- renderTable({
    req(input$gene_selection)
    cancer = dbFetch(dbSendQuery(con, 'SELECT `Mutational.cancer.driver.genes`, `COSMIC.somatic.mutations.in.cancer.genes` 
                                 FROM `Risk factors` 
                                 WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    core_fitness = dbFetch(dbSendQuery(con, 'SELECT * FROM `Core fitness` 
                                       WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)), n = 1)
    core_fitness = core_fitness[,5:ncol(core_fitness)]
    df = as.data.frame(t(cbind(cancer, core_fitness))) %>% rownames_to_column("stat")
    colnames(df) = c("Group membership", " ")
    df$`Group membership` = gsub('\\.', " ", df$`Group membership`)
    df
  }, digits=0, align='l')
  
  output$genetic_disorders_table <- renderTable({
    req(input$gene_selection)
    disorders = dbFetch(dbSendQuery(con, 'SELECT `Disease.name`, `Disease.ID`, `Genetic.association`
                                    FROM `Rare disease associations` 
                                    WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)))
    validate(need(!all(is.na(disorders)), "No genetic disorders associated with this gene."))
    validate(need(nrow(disorders) != 0, "No genetic disorders associated with this gene."))
    disorders
  }, digits=0)
  
  output$hpo_phen_table <- renderPlotly({
    req(input$gene_selection)
    tex = dbFetch(dbSendQuery(con, 'SELECT * FROM `Risk factors` 
                              WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)))
    tex = tex[,17:ncol(tex)]
    
    validate(
      # need(!all(is.na(tex)), "No data on affected tissues available for this gene.")
      need(nrow(tex) != 0, "No data on affected tissue classes available for this gene.")
    )
    
    df = t(tex)
    df[is.na(df)] <- 0
    df = as.data.frame(df) %>% rownames_to_column("Tissue")
    colnames(df) = c("Tissue", "Present")
    df$Placeholder = " "
    
    # validate(
    #   # need(!all(is.na(tex)), "No data on affected tissues available for this gene.")
    #   need(!all(df$Present == 0), "No tissue classes affected by diseases associated with this gene.")
    # )
    
    # Dynamic colours that allow us to accurately represent data, no matter what's coming up!
    colpal = c(ifelse(min(df$Present) == 0, "white", "blue"),
               ifelse(max(df$Present) == 0, "white", "blue"))
    
    plot_ly(x=df$Tissue,
            y=df$Placeholder,
            z=df$Present,
            type = "heatmap",
            height = 170,
            width = 400,
            showscale = FALSE,
            colors = colorRamp(colpal)) %>%
      layout(margin = list(l = 50, r = 70, b = 30, t = 120, pad = 8),
             xaxis = list(side ="top", tickangle=-45))
  })
  
  output$antitargets <- renderPlot({
    req(input$gene_selection)
    antitgts = dbFetch(dbSendQuery(con, 'SELECT `Antitarget.gene.name`, `Protein.seq..identity....` FROM `Antitargets` 
                                   WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)))
    validate(need(!all(is.na(antitgts$Antitarget)), "No antitargets found for this gene."))
    validate(need(nrow(antitgts) != 0, "No antitargets found for this gene."))
    colnames(antitgts) = c("Antitarget", "identity")
    antitgts = antitgts[order(antitgts$identity, decreasing = TRUE),]
    ggplot(antitgts, aes(x=factor(Antitarget, level = antitgts$Antitarget), y=identity)) +
      geom_bar(stat="identity", fill="blue") +
      # ylim(40, 100) +
      xlab("Antitarget") + 
      ylab("% protein sequence identity") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      ggtitle("Antitargets passing protein sequence identity threshold (>40%)") +
      coord_cartesian(ylim=c(40,100))
  })
  
  output$mouse_phen_table <- renderPlotly({
    req(input$gene_selection)
    # if(!is.null(input$gene_selection)) {
      high_level_phenotypes = c('adipose tissue', 'behaviour/neurological', 'cardiovascular system', 'cellular',
                                'craniofacial', 'digestive/alimentary system', 'embryo', 'endocrine/exocrine glands',
                                'growth/size/body', 'hearing/vestibular/ear', 'hematopoietic system', 'homeostasis/metabolism',
                                'integument', 'immune system', 'limbs/digits/tail', 'liver/biliary system', 'mortality/aging',
                                'muscle', 'nervous system', 'pigmentation', 'renal/urinary system', 'reproductive system',
                                'respiratory system', 'skeleton', 'taste/olfaction', 'neoplasm', 'vision/eye')
      # We want ALL of those columns in the plot, even if 0.
      # data = xl$`MGI mouse` %>% filter(HGNC.Name == input$gene_selection)
      df = dbFetch(dbSendQuery(con, 'SELECT * FROM `MGI Mouse` WHERE `HGNC.Name` == :x', params = list(x = input$gene_selection)))
      validate(
        need(nrow(df) != 0, "No mouse phenotypes on record for this gene."),
        need(df$MGI.Gene.Marker.ID[1] != "No phenotypes found", "No mouse phenotypes on record for this gene.")
      )
      validate(
        need(!all(is.na(df)), "No mouse phenotypes on record for this gene.")
      )
      # The top-level phenotypes may occasionally appear as comma-delimited lists, due to a given phenotype appearing in
      # more than one top-level group. This code straightens that out into having one row per top-level group for easier counting.
      df = df %>% mutate(Top.level.phenotype = strsplit(Top.level.phenotype, ", ")) %>% unnest(Top.level.phenotype)
      # NOW we can start counting stuff up.
      counts = as.data.frame(table(df$`Top.level.phenotype`))
      cmon = as.data.frame(setdiff(high_level_phenotypes, counts$Var1))
      cmon$Freq = 0
      colnames(cmon)[1] = "Var1"
      counts = rbind(counts, cmon)
      counts$Var1 = as.character(counts$Var1)
      counts = counts[order(counts$Var1),]
      counts$Present = as.integer(as.logical(counts$Freq))
      counts$Placeholder = " "
      # That's got it. Now do the plot. This might be a job for Plotly, given what it can do with hover.
      
      # Dynamic colours that allow us to accurately represent data, no matter what's coming up!
      colpal = c(ifelse(min(counts$Present) == 0, "white", "blue"),
                 ifelse(max(counts$Present) == 0, "white", "blue"))
      
      plot_ly(x=counts$Var1,
              y=counts$Placeholder,
              z=counts$Present,
              type = "heatmap",
              hoverinfo = "text",
              text = ~paste('N:', counts$Freq, sep = ""),
              height = 240,
              showscale=FALSE,
              colors = colorRamp(colpal)) %>%
        layout(margin = list(l = 50, r = 50, b = 20, t = 200, pad = 4),
               xaxis = list(side ="top", tickangle=-45))
      # Can I get some hover text? It seems to not want to play along.
    # }
  })
  
})

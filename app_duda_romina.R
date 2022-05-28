
# Llamar librerias


library(pubmed.mineR)
library(easyPubMed)
library(lsa)
library(bibliometrix)
library(pubmedR)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(rgl)

# Información MD previa

Disease = c("Diabetes", "Obesity", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing’s syndrome", "Hyperuricemia", "Hemochromatosis", "Metabolic syndrome", "Fatty liver disease", "Hypercholesterolemia", "Metabolic disease")

Disease_one_word = c("Diabetes", "Obesity", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing", "Hyperuricemia", "Hemochromatosis", "Metabolic", "Fatty liver", "Hypercholesterolemia")

Disease_one_word_without_Obesity = c("Diabetes", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing", "Hyperuricemia", "Hemochromatosis", "Metabolic", "Fatty liver", "Hypercholesterolemia")

Mesh = c("Diabetes Mellitus", "Obesity", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing's syndrome", "Hyperuricemia", "Hemochromatosis", "Metabolic syndrome", "Fatty liver", "Hypercholesterolemia", "Metabolic disease")

MD_list = data.frame(Disease, Mesh)
MD_list = data.frame(lapply(MD_list, as.character), stringsAsFactors = FALSE)

choices = setNames(MD_list$Mesh, MD_list$Disease)

Microbiota = c("Dysbiosis", "Short chain fatty acids", "Bacteria-related molecules", "Serratia", "Enterobacter", "Morganella", "Skunalikevirus", "Phifllikevirus", "Roseburia", "Blautia", "Clostridium", "Akkermansia", "Ruminococcus", "Lactobacillus", "Decreased microbial richness", "Microbial diversity", "Firmicutes:Bacteroidetes ratio", "Microbial metabolites", "Metagenomic analysis", "Butyrate producing bacteria", "Firmicutes", "Bacteroidetes", "Butyrate-producing bacteria", "Bacteroides", "Prevotella", "Xylanibacter", "Proteobacteria", "Fecal microbiota", "Propionate", "Acetate", "Christensenellaceae", "Tenericutes", "Prevotella:Bacteroides ratio", "Gut microbes", "Microbiota", "Gastrointestinal microbiome", "Verrucomicrobia", "Lachnospiraceae", "Bifidobacterium", "Fusobacterium", "Faecalibacterium", "Roseburia", "Eubacterium", "Bilophila", "Desulfovibrio", "Blautia", "Turicibacter", "Bilophila", "Adlercreutzia", "Actinobacteria", "Streptococcus", "Lactic acid bacteria", "Lachnospiraceae", "Rikenellaceae", "Parasutterella", "Sutterella", "Lachnospiraceae", "Veillonellaceae", "Alistipes", "Actinobacteria", "Enterobacteriaceae", "Staphylococcus aureus", "Escherichia coli", "Methabacteriodes", "Betaproteobacteria", "Firmicutes", "Clostridia", "Proteobacteria", "Tenericutes", "Actinobacteria", "Mollicutes", "Negativicutes", "Bacteroidia", "Erysipelotrichales", "Selenomonadales", "Bacteroidales", "Coriobacteriales", "Clostridiales", "Prevotellaceae", "Lachnospiraceae", "Peptostreptococcaceae", "Ruminococcaceae", "Coriobacteriaceae", "Collinsella")

Microbiota_one_word = c("Dysbiosis", "acids", "Serratia", "Enterobacter", "Morganella", "Skunalikevirus", "Phifllikevirus", "Roseburia", "Blautia", "Clostridium", "Akkermansia", "Ruminococcus", "Lactobacillus", "microbial ", "diversity", "metabolites", "Metagenomic", "Butyrate", "Firmicutes", "Bacteroidetes", "Bacteroides", "Prevotella", "Xylanibacter", "Proteobacteria", "microbiota", "Propionate", "Acetate", "Christensenellaceae", "Tenericutes", "Prevotella", "microbes", "microbiome", "Verrucomicrobia", "Lachnospiraceae", "Bifidobacterium", "Fusobacterium", "Faecalibacterium", "Roseburia", "Eubacterium", "Bilophila", "Desulfovibrio", "Blautia", "Turicibacter", "Bilophila", "Adlercreutzia", "Actinobacteria", "Streptococcus", "Lactic acid bacteria", "Lachnospiraceae", "Rikenellaceae", "Parasutterella", "Sutterella", "Lachnospiraceae", "Veillonellaceae", "Alistipes", "Actinobacteria", "Enterobacteriaceae", "Staphylococcus", "Escherichia", "Methabacteriodes", "Betaproteobacteria", "Firmicutes", "Clostridia", "Proteobacteria", "Tenericutes", "Actinobacteria", "Mollicutes", "Negativicutes", "Bacteroidia", "Erysipelotrichales", "Selenomonadales", "Bacteroidales", "Coriobacteriales", "Clostridiales", "Prevotellaceae", "Lachnospiraceae", "Peptostreptococcaceae", "Ruminococcaceae", "Coriobacteriaceae", "Collinsella")

# UI

header = dashboardHeader(title = "GUTMetabolic | Explore genes and microbiota in metabolic diseases", titleWidth = 650)
sidebar = dashboardSidebar(width = 300, sidebarMenu(id = "tabs", menuItem("Publication Search", tabName = "a", icon = icon("search")), conditionalPanel("input.tabs == 'a'", selectInput(inputId = "busqueda", label = "Metabolic disease:", choices = choices, selected = NULL, multiple = FALSE, selectize = TRUE, width = NULL, size = NULL), dateInput(inputId = "fechadesde", "Since:", value = "2019-01-01", format = "yyyy-mm-dd"), dateInput("fechahasta", "Until:", value = "2021-01-01", format = "yyyy-mm-dd"), tags$div(id = "espacio1", br()), actionButton(inputId = "downPub", label = "Search", color = "success", style = "gradient", block = FALSE, size = "md"), tags$div(id = "espacio3", br()), hidden(actionButton(inputId = "nuevaBusq", label = "New Search", color = "success", style = "gradient", block = FALSE, size = "md"))), menuItem("Explore Publications", tabName = "b", icon = icon("chart-bar")), conditionalPanel("input.tabs == 'b'"), menuItem("Gene-Disease Relation", tabName = "c", icon = icon("dna")), menuItem("Microbiota-Disease Relation", tabName = "d", icon = icon("bacterium"))))

body = dashboardBody(
  useShinyjs(),
  tags$style(type = "text/css", 
             ".shiny-output-error {visibility:hidden;}", 
             ".shiny-output-error:before {visibility:hidden;"
             ),
  fluidRow(
    tabItems(
      tabItem(
        tabName = "a", h2("Find publications of metabolic diseases in PubMed", align = "center"),
        br(),
        box(id = "resultado",
            title = "Search result", solidHeader = TRUE, status = "success", collapsible = TRUE, width = 11,
            tags$div(style = 'cursor:pointer', dataTableOutput("pmidresult")),
            br(),
            wellPanel(htmlOutput("rows_selected"))
        )
      )
      ,
      tabItem(
        tabName = "b", h2("Genes and word clouds", align = "center"),
        br(),
        fluidRow(         
          column(6,
                 box(id= "Word Cloud", title = "Word Cloud", width = 11,  solidHeader = TRUE, status = "success", plotOutput("nubepalabras"), hidden(dataTableOutput("tablapalabras")))),
          column(6,
                 box(id = "Gene Cloud", title = "Gene Cloud", width = 11, solidHeader = TRUE, status = "success", plotOutput("nubegenes"), hidden(dataTableOutput("tablagenes"))))),
        fluidRow(      
        column(6,br(),
               actionButton(inputId = "VerTabla", label = "See Table", color= "success", style = "gradient", icon = icon("table"),
                              block = FALSE, size = "md"), br(), hidden(actionButton(inputId = "VolverGrafico", label = "Back", color = "success", style = "gradient", icon = icon("redo"), block = FALSE, size = "md"))),
        column(6,br(),  
               actionButton(inputId = "VerTabla", label = "See Table", color= "success", style = "gradient", icon = icon("table"),
                            block = FALSE, size = "md"), br(), hidden(actionButton(inputId = "VolverGrafico", label = "Back", color = "success", style = "gradient", icon = icon("redo"), block = FALSE, size = "md")))),
        fluidRow( 
        column(6, br(),
                 box(id="parametroswordcloud", title = "Parameters Word Cloud", width = 11, solidHeader = TRUE, status = "primary", sliderInput("minpalabras", label = "Word frequency Min.:", min = 1, max = 15, value = 1),
                 sliderInput("maxpalabras", label = "Word frequency Max.:", min = 1, max = 500, value = 100))),
        column(6,br(),
                 box(id="parametroswordcloud", title = "Parameters Gene Cloud", width = 11, solidHeader = TRUE, status = "primary", sliderInput("mingenes", label = "Gene frequency Min.:", min = 1, max = 15, value = 1),
                                 sliderInput("maxgenes", label = "Gene frequency Max.:", min = 1, max = 500, value = 100))))
                                         )
        ,
        tabItem(
          tabName = "c", h2("Relation Genes - Obesity", align = "center"),
          fluidRow(br(),
                   column(4,
                          box(id = "seleccMD", title = "Parameters", solidHeader = TRUE, status = "primary", width = 11, 
                              dateInput(inputId = "fechadesde2", "Since:", value = "2019-01-01", format = "yyyy-mm-dd"),
                              dateInput(inputId = "fechahasta2", "Until:", value = "2021-01-01", format = "yyyy-mm-dd"),
                              actionButton(inputId = "downPub2", label = "Search", color="success", style = "gradient", icon= icon("download")),
                              tipify(actionButton(inputId = "Instrucciones", label = "", color = "success", style = "gradient", icon = icon("question-circle")),
                                     "xxxxxxxxxxxx", trigger = "click")
                              ,
                              tags$div(id="espacio5", br()),
                              actionButton(inputId = "nuevaBusqObesity", label = "Find Publications about Gene", color = "success", style = "gradient", icon=icon("download")),
                              tags$div(id="espacio6", br()),
                              hidden(actionButton(inputId = "VolverObesity", label = "Back to List", color = "success", style = "gradient", icon=icon("redo")))
                              )),
                   column(7,
                          box(id="relacionesgenes3D", title = "Obesity-Gene Relation Graph", solidHeader = TRUE, status = "success", width = 30, rglwidgetOutput("relacionesgenes3D"))),br(),
                   fluidRow(br(),
                     column(12,
                     box(id = "Related Genes Table",title = "Related Gene Table", solidHeader = TRUE, status = "success", width = 11, tags$div(style='cursor:pointer', dataTableOutput("pmidresult2"), br(), hidden(textOutput("nTextObesity_new")), hidden(textOutput("genes_warn")),
                                                              br(), hidden(dataTableOutput("pmidresultObesitynew")), br(), hidden(htmlOutput("rows_selected_obesitynew")))))),
                     
          ))
        ,
        tabItem(
          tabName = "d", h2("Relation Microbiota - Obesity", align = "center"),
          fluidRow(br(),
                   column(4,
                          box(id = "seleccMD", title = "Parameters", solidHeader = TRUE, status = "primary", width = 11, 
                              dateInput(inputId = "fechadesde2", "Since:", value = "2019-01-01", format = "yyyy-mm-dd"),
                              dateInput(inputId = "fechahasta2", "Until:", value = "2021-01-01", format = "yyyy-mm-dd"),
                              actionButton(inputId = "downPub2", label = "Search", color="success", style = "gradient", icon= icon("download")),
                              tipify(actionButton(inputId = "Instrucciones", label = "", color = "success", style = "gradient", icon = icon("question-circle")),
                                     "xxxxxxxxxxxx", trigger = "click")
                              ,
                              tags$div(id="espacio5", br()),
                              actionButton(inputId = "nuevaBusqObesity", label = "Find Publications about Gene", color = "success", style = "gradient", icon=icon("download")),
                              tags$div(id="espacio6", br()),
                              hidden(actionButton(inputId = "VolverObesity", label = "Back to List", color = "success", style = "gradient", icon=icon("redo")))
                          )),
                   column(7,
                          box(id="relacionesgenes3D", title = "Obesity-Gene Relation Graph", solidHeader = TRUE, status = "success", width = 30, rglwidgetOutput("relacionesgenes3D"))),br(),
                          fluidRow(br(),
                                   column(12,
                                          box(id = "Related Microbiota Table",title = "Related Microbiota Table", solidHeader = TRUE, status = "success", width = 11, tags$div(style='cursor:pointer', dataTableOutput("pmidresult2"), br(), hidden(textOutput("nTextObesity_new")), hidden(textOutput("genes_warn")),
                                                                                                                                                                    br(), hidden(dataTableOutput("pmidresultObesitynew")), br(), hidden(htmlOutput("rows_selected_obesitynew")))))),
                          
                   ))
          )))
        

    
ui = dashboardPage(title = "GUTMetabolic", header, skin = "green", sidebar, body)


# Server

server = function(input, output){
  # Buscar publicaciones
    nTexto <- eventReactive(input$createPub,{paste(c(input$busqueda,"[MH] AND", as.character(input$fechadesde, "%Y/%m/%d"),":", as.character(input$fechahasta, "%Y/%m/%d", "[DP]", collapse = "")))})
                                                 output$nText = renderText(nTexto())
                                                 pubmedResult = eventReactive(input$downPub, 
                                                                              {withProgress(message = "Downloading information from PubMed", detail = "Please, wait :)", value = 0, {incProgress(2/10, detail = "Please, wait :)")
                                                                                corpus = batch_pubmed_download(pubmed_query_string = paste(c(input$busqueda,"[MH] AND", as.character(input$fechadesde, "%Y/%m/%d"),
                                                                                                                                             ":", as.character(input$fechahasta, "%Y/%m/%d", "[DP]", collapse = ""),
                                                                                                                                             format("abstract"),
                                                                                                                                             batch_size = 1500)))
                                                                                                                                           incProgress(5/10, detail = "Consolidating files...")
                                                                                                                                           file.create("pubmed_result.txt")
                                                                                                                                           for(i in 1:length(corpus)){
                                                                                                                                             file.append("pubmed_result.txt", corpus[i])}
                                                                                                                                           corpus_output = readabs("pubmed_result.txt")
                                                                                                                                           incProgress(10/10, detail = "Ready!")
                                                                                                                                           return(corpus_output)
                                                                              })
                                                                              }
                                                                                )
                                                                                output$pmidresult = renderDataTable({corpus = pubmedResult()
                                                                                pmidAbs = data.frame(abstracts=corpus@Abstracts)
                                                                                titulos = c()
                                                                                for(i in 1:length(corpus@PMID)){
                                                                                abstractnro = i
                                                                                text = unlist(strsplit(corpus@Abstract[abstractnro], "\\."))
                                                                                op = which(unlist(lapply(unlist(strsplit(corpus@Abstract[abstractnro], "\\.")), wordcount))>=6)[1]
                                                                                           titulos = c(titulos, text[op])
                                                                                           }
                                                                                Intro = unlist(lapply(pmidAbs$abstracts, substr, start=1, stop=300))
                                                                                pmidRes = data.frame(PMID = corpus@PMID, Titulo = titulos)
                                                                                datatable(pmidRes, selection=list(mode = 'single', selected = 1), options = list(lenghtMenu = c(5,30,50), pageLenght=5)) %>% formatStyle(names(pmidRes), cursor = "hand", target = "cell")
                                                                                })
observeEvent(input$downPub,{
  hide("nText")
  hide("nText2")
  hide("boxdebusqueda")
  shinyjs::toggle(id = "nuevaBusq")
})
observeEvent(input$nuevaBusq,{
  shinyjs::toggle(id = "nText")
  shinyjs::toggle(id = "nText2")
  shinyjs::toggle(id = "boxdebusqueda")
  shinyjs::toggle(id = "nuevaBusq")
  hide("nuevaBusq")
  reset("busqueda")
  reset("fechahasta")
  reset("fechadesde")
})
observeEvent(input$VerTabla,{
  hide("nuvegenes")
  hide("nuvepalabras")
  hide("maxpalabras")
  hide("minpalabras")
  shinyjs::toggle(id = "tablagenes")
  shinyjs::toggle(id = "tablapalabras")
  shinyjs::toggle(id = "VolverGrafico")
}) 
observeEvent(input$VolverGrafico,{
  shinyjs::toggle(id = "nuvegenes")
  shinyjs::toggle(id = "nuvepalabras")
  shinyjs::toggle(id = "nText2")
  shinyjs::toggle(id = "minpalabras")
  shinyjs::toggle(id = "maxpalabras")
  hide("tablagenes")
  hide("VolverGrafico")
})
observeEvent(input$nuevaBusqObesity,{
  hide("instrucciones")
  hide("pmidresult2")
  shinyjs::toggle(id = "VolverObesity")
  shinyjs::toggle(id = "pmidresultObesitynew")
  shinyjs::toggle(id = "rows_selected_obesitynew")
  shinyjs::toggle(id = "nTextObesity_new")
  shinyjs::toggle(id = "genes_warn")
})                                                                                
observeEvent(input$VolverObesity,{
  hide("genes_warn")
  hide("nTextObesity_new")
  hide("pmidresultObesitynew")
  hide("rows_selected_obesitynew")
  shinyjs::toggle(id = "instrucciones")
  shinyjs::toggle(id = "pmidresult2")
  hide("VolverObesity")
})    
output$rows_selected = renderPrint({
  corpus = pubmedResult()
  text = unlist(strsplit(corpus@Abstract[input$pmidresult_rows_selected],"\\."))
  op = which(unlist(lapply(unlist(strsplit(corpus@Abstract[input$pmidresult_rows_selected],"\\.")), wordcount))>=6)[1]
  titulos = text[op]
  for (i in 1:length(text)){
    if (i==op)cat(paste('<p><h4>','<font color=\"#B404F6\"><b>',text[i], '</b></font>', '</h4>', '\n', '<p>'), fill = TRUE)
    else cat(paste('<SPAN class=sentence><i>', text[i],".", '</i></SPAN>'), fill = TRUE)
    }
  cat(paste("<p><a href='https://pubmed.ncbi.nlm.nih.gov/", corpus@PMID[pmidresult_rows_selected]," 'target= blank>", "Link to publication in Pubmed","</a><p>"))
  })

# Explorar publicaciones:

output$nubepalabras = renderPlot(width = "auto", height = "auto", {
  withProgress(message = "Processing corpus...", detail = "Please, wait :)",
               value = 0, {incProgress(2/10, detail = "Please, wait :)")
                 corpus = pubmedResult()
                 words = word_associations(corpus)
                 row.names(words) = NULL
                 incProgress(8/10, detail = "Generating word cloud...")
                 wordcloud(words$words, words$Freq, max.words = input$maxpalabras, min.freq = input$minpalabras, random.color = TRUE, colors = brewer.pal(8, "Dark2"), scale = c(4, 0.5))
                 })})

getFrecGene = reactive({
  withProgress(message = "Processing corpus...", detail = "Please, wait :)",
               value = 0, {incProgress(2/10, detail = "Please, wait :)")
  corpus = pubmedResult()
  incProgress(8/10, detail = "Generating gene cloud...")
  frecGenes = data.frame(gene_atomization(corpus), stringsAsFactors = FALSE)
  return(frecGenes)
  })})
  
wordcloudRep = repeatable(wordcloud)
output$nubegenes = renderPlot(width = "auto", height = "auto", {
  freqGenes = getFrecGene()
  wordcloudRep(words = freqGenes$Gene_symbol, freq = as.numeric(freqGenes$Freq), max.words = input$maxpalabras, min.freq = input$minpalabras, random.color = TRUE, colors = brewer.pal(8, "Dark2"), scale = c(8,0.6))
  })
 
output$tablagenes = renderDataTable({
  freqGenes = getFrecGene()
  colnames(freqGenes) = c("Symbol", "Name", "Frequency")
  datatable(freqGenes, rownames = FALSE)
})

output$tablapalabras = renderDataTable({
  corpus = pubmedResult()
  words = word_atomizations(corpus)
  colnames(words) = c("Word", "Frequency")
  datatable(words, rownames = FALSE)
})

# Relación gen - Enfermedad:

pubmedResult2 = eventReactive(input$downPub2,
                              {withProgress(message = "Downloading information from PubMed", detail = "Please, wait :)",
                                            value = 0, {incProgress(0.1, detail = "Please, wait :)")
                                              corpusMD = batch_pubmed_download(pubmed_query_string = paste(c("Metabolic Diseases [MH], AND", as.character(input$fechadesde2, "%Y/%m/%d"), ":", as.character(input$fechahasta2, "%Y/%m/%d"), "[DP]", collapse = ""),
                                                                                                           format = "abstract", batch_size = 1500, dest_file_prefix = "corpusMD_"))
                               incProgress(0.1, detail = "Consolidating files...")
                               file.create("pubmed_result2.txt")
                                                                  for (i in 1:lenght(corpusMD)){
                                                                                 file.append("pubmed_result2.txt", corpusMD[i])}
                               corpus_query_comb = readabs("pubmed_result2.txt")
                                                                               incProgress(0.2, detail = "Extracting genes...")
                                                                               genes_MD_comb = data.frame(gene_atomization(corpus_query_comb, stringAsFactor = FALSE))
                                                                               genes_MD_comb$Freq = as.numeric(genes_MD_comb$Freq)
                                                                               genes_MD_list_comb = genes_MD_comb[genes_MD_comb$Freq>10]
                                                                               terms_obesity = c("bmi", "metabolic", "diabetes", "glucose", "liver", "adipose")
                                                                               incProgress(0.1, detail = "LSA")
                                                                               tdm_main = tdm_for_lsa(corpus_query_comb, c(genes_MD_list_comb, "obesity", terms_obesity, Disease_one_word_without_Obesity))
                                                                               lsa_corpus = lsa(tdm_main, dims = dimcalc_share())
                                                                               matriz = as.textmatrix(lsa_corpus)
                                                                               aa = associate(matriz, "obesity", measure = "cosine", threshold = 1e-399)
                                                                               coseno_out = data.frame(genes = names(aa), coseno= aa, stringsAsFactors = FALSE)
                                                                               incProgress(0.1, detail = "Obesity publications...")
                                                                               output_obesity = batch_pubmed_download(pubmed_query_string = paste(c("Obesity [MH] AND", as.character(input$fechadesde2, "%Y/%m/%d"),
                                                                                                                                                    ":", as.character(input$fechahasta2, "%Y/%m/%d"), "[DP]", collapse = ""),
                                                                                                                                                  format = "abstract", batch_size = 1500, dest_file_prefix = "pubmedobesity"))
                                                                                })})
file.create("pubmed_result_obesity.txt")
for(i in 1:length(output_obesity)){
  file.append("pubmed_result_obesity.txt", output_obesity[i])}
corpus_obesity = readabs("pubmed_result_obesity.txt")
incProgress(0.1, detail = "Obesity-associated genes")
genes_MD = gene_atomization(corpus_obesity)
genes_MD_list = genes_MD[,1]
genes_MD_list_obesity = data.frame(genes = genes_MD_list, obesity = "yes", stringsAsFactors = FALSE)
result_genes = merge(coseno_out, genes_MD_list_obesity, by.x = "genes", by.y = "genes", all.x = TRUE)
result_genes = result_genes[!result_genes$genes %in% c(terms_obesity, Disease_one_word_without_Obesity),]
funcion_relMD = function(x){combine_words(names(associate(matrix, x, measure = "cosine", threshold = 0.4))[which(names(associate(matriz,x, measure = "cosine", threshold = 0.4)) %in% Disease_one_word_without_Obesity)], sep = ",", and = "and")}
result_genes$other_MD = lapply(result_genes$genes, function_relMD)
result_genes[result_genes$other_MD == "character(0)", 4] = NA
result_genes$coseno = round(result_genes$coseno, 5)
incProgress(0.3, detail = "Complete")
row.names(result_genes) = seq(1:length(result_genes$genes))
mylist = list("matriz" = matriz, "result_genes" = result_genes, "terms" = terms_obesity)
return(mylist)

output$pmidresult2 = renderDataTable({
  result_genes = pubmedResult2()$result_genes
  datatable(result_genes, colnames = c("Gene", "Cosine similarity with obesity", "Mentioned in obesity publications", "Other MD related"),
            selection = list(mode='single', selected = 1), options = list(lenghtMenu = c(5,30,50), pageLength = 5, order = list(list(1, 'desc'))),
            rownames = FALSE) %>% formatStyle(names(result_genes), cursor ="hand", target = "cell")%>% formatStyle(columns = 4, fontSize = '80%')%>% formatStyle(c("Genes", "Cosine", "Obesity", "Other MD"), textAlign = 'center')
})

output$relacionesgenes3D = renderRglwidget({
  matriz = pubmedResult2()$matriz
  terms_obesity = pubmedResult2()$terms
  termMatrix2 = matriz%*%t(matriz)
  exclude = which(row.names(termMatrix2) %in% terms_obesity)
  include = which(row.names(termMatrix2)%in% c((pubmedResult2()$result_genes)$genes, "obesity"))
  plot_neighbors("obesity", n=length(include), tvectors = matriz[include,include], legend = T, connect.lines = 0, alpha = c(0.5,0.5), col = c("lemonchiffon", "orange", "darkred"))
  rglwidget()
})

nTextoObesity_new = eventReactive(input$nuevaBusqObesity, {
  result_genes = pubmedResult2()$result_genes$genes
  selected_gene = result_genes[input$pmidresult2_rows_selected]
  paste(c("Your search is:", as.character(selected_gene[1]), "AND Obesity"), collapse = "")})

output$nTextObesity_new = renderText(nTextoObesity_new())

pubmedResultObesity_and_gene = eventReactive(input$nuevaBusqObesity, {
  withProgress(message = "Downloading information from PubMed...", detail = "Please wait :)", value = 0, {incProgress(2/10, detail = "Please wait :)")
    result_genes = pubmedResult2()$result_genes$genes
    selected_gene = result_genes[input$pmidresult2_rows_selected]
    corpus_obesity_new = batch_pubmed_download(pubmed_query_string = paste(c(as.character(selected_gene[1]), "AND Obesity"), collapse = ""), 
                                               format = "abstract", batch_size = 1500, dest_file_prefix = "corpus_Obesity_new")
    incProgress(5/10, detail = "Consolidating files...")
    file.create("pubmed_Obesity_new.txt")
    for (i in 1:lenght(corpus_obesity_new)){
      file.append("pubmed_Obesity_new.txt", corpus_obesity_new[i])}
    
    corpus_output_Obesity_new = readabs("pubmed_Obesity_new.txt")
    incProgress(10/10, detail = "Complete")
    mylist = list("corpus_output_Obesity_new"=corpus_output_Obesity_new, "corpus_Obesity_new" = corpus_Obesity_new)
    return(mylist)
    })
})

output$pemidresultObesitynew = renderDataTable({
  corpusObesitynew = pubmedResultObesity_and_gene()$corpus_output_Obesity_new
  pmidAbs_Obesitynew = data.frame(abstracts = corpusObesitynew@Abstract)
  titulos = c()
  for (i in 1:lenght(corpusObesitynew@PMID)) {
    abstractnro = i
    text = unlist (strsplit(corpusObesitynew@Abstract[abstractnro], "\\."))
    op = which(unlist(lapply(unlist(strsplit(corpusObesitynew@Abstract[abstractnro], "\\.")), wordcount))>=6)[1]
    titulos = c(titulos, text[op])
  }
  Intro = unlist(lapply(pmidAbs_Obesitynew$abstracts, substr, start = 1, stop = 300))
  pmidRes = data.frame(PMID = corpusObesitynew@PMID, Titulos = titulos)
  datatable ( pmidRes, selection = list(mode='single', selected = 1), options = list(lenghtMenu = c(5, 30, 50), pageLenght = 5),
              rownames = FALSE) %>% formatStyle(names(pmidRes), cursor = 'hand', target = "cell")
  
})

genes_warn = eventReactive(input$nuevaBusqObesity, {
  result_genes = pubmedResult2()$result_genes$genes
  selected_gene = result_genes[input$pmidresult2_rows_selected]
  corpus_Obesity_new = batch_pubmed_download(pubmed_query_string = paste(c(as.character(selected_gene[1]), "AND Obesity"), collapse = ""),
                                             format = "abstract", batch_size = 1500, dest_file_prefix = "corpus_Obesity_new")
  paste(ifelse(length(corpus_Obesity_new)==0, "Publications not found"))
  
})

output$genes_warn = renderText(genes_warn())
output$rows_selected_obesitynew = renderPrint({
  corpusObesitynew = pubmedResultObesity_and_gene()$corpus_output_Obesity_new
  text = unlist(strsplit(corpusObesitynew@Abstract[input$pmidresultObesitynew_rows_selected], "\\."))
  op = which(unlist(lapply(unlist(strsplit(corpusObesitynew@Abstract[input$pmidresultObesitynew_rows_selected], "\\.")), wordcount))>=6)[1]
  titulos = text[op]
  for(i in 1:length(text)){
    if (i == op)cat(paste('<p><h4>', '<font color=\"#4B04F6\"><b>', text[i],'</b></font>','</h4>','\n','<p>'), fill = TRUE)
    else cat(paste('<SPAN class=sentence><i>', text[i], ".",'</i></SPAN>'), fill = TRUE)
  }
  cat(paste("<p><a
            href='https://pubmed.ncbi.nlm.nih.gov/", corpusObesitynew@PMID[input$pmidresultObesitynew_rows_selected], "'taget=_blank>", "Link to PubMed publications", "</a><p>"))
})

# Relación Microbiota - Enfermedad: Pendiente.

}

# App                                                                              

shinyApp(ui, server)


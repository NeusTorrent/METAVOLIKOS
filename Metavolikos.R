
# Packages

if (!require("pubmed.mineR")) install.packages("pubmed.mineR")
if (!require("easyPubMed")) install.packages("easyPubMed")
if (!require("lsa")) install.packages("lsa")
if (!require("bibliometrix")) install.packages("bibliometrix")
if (!require("shiny")) install.packages("shiny")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinyBS")) install.packages("shinyBS")
if (!require("rgl")) install.packages("rgl")
if (!require("ngram")) install.packages("ngram")
if (!require("dplyr")) install.packages("dplyr")
if (!require("DT")) install.packages("DT")
if (!require("wordcloud")) install.packages("wordcloud")
if (!require("V8")) install.packages("V8")
if (!require("RISmed")) install.packages("RISmed")
if (!require("LSAfun")) install.packages("LSAfun")
if (!require("igraph")) install.packages("igraph")
if (!require("knitr")) install.packages("knitr")
if (!require("stringr")) install.packages("stringr")

# Preliminary arrangements

Disease = c("Diabetes", "Obesity", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing’s syndrome", "Hyperuricemia", "Hemochromatosis", "Metabolic syndrome", "Fatty liver disease", "Hypercholesterolemia", "Metabolic disease")

Disease_one_word = c("Diabetes", "Obesity", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing", "Hyperuricemia", "Hemochromatosis", "Metabolic", "Fatty liver", "Hypercholesterolemia")

Disease_one_word_without_Obesity = c("Diabetes", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing", "Hyperuricemia", "Hemochromatosis", "Metabolic", "Fatty liver", "Hypercholesterolemia")

Mesh = c("Diabetes Mellitus", "Obesity", "Dyslipidemia", "Hypolipidemia", "Hyperlipidemia", "Hyperthyroidism", "Hypothyroidism", "Hypoparathyroidism", "Hyperparathyroidism", "Cushing's syndrome", "Hyperuricemia", "Hemochromatosis", "Metabolic syndrome", "Fatty liver", "Hypercholesterolemia", "Metabolic disease")

MD_list = data.frame(Disease, Mesh)
MD_list = data.frame(lapply(MD_list, as.character), stringsAsFactors = FALSE)

choices = setNames(MD_list$Mesh, MD_list$Disease)

Microbiota_one_word = c("dysbiosis", "acids", "serratia", "enterobacter", "morganella", "skunalikevirus", "phifllikevirus", "roseburia", "blautia", "clostridium", "akkermansia", "ruminococcus", "lactobacillus", "microbial ", "diversity", "metabolites", "metagenomic", "butyrate", "firmicutes", "bacteroidetes", "bacteroides", "prevotella", "xylanibacter", "proteobacteria", "microbiota", "propionate", "acetate", "christensenellaceae", "tenericutes", "prevotella", "microbes", "microbiome", "verrucomicrobia", "lachnospiraceae", "bifidobacterium", "fusobacterium", "faecalibacterium", "roseburia", "eubacterium", "bilophila", "desulfovibrio", "blautia", "turicibacter", "bilophila", "adlercreutzia", "actinobacteria", "streptococcus", "lactic acid bacteria", "lachnospiraceae", "rikenellaceae", "parasutterella", "sutterella", "lachnospiraceae", "veillonellaceae", "alistipes", "actinobacteria", "enterobacteriaceae", "staphylococcus", "escherichia", "methabacteriodes", "betaproteobacteria", "firmicutes", "clostridia", "proteobacteria", "tenericutes", "actinobacteria", "mollicutes", "negativicutes", "bacteroidia", "erysipelotrichales", "selenomonadales", "bacteroidales", "coriobacteriales", "clostridiales", "prevotellaceae", "lachnospiraceae", "peptostreptococcaceae", "ruminococcaceae", "coriobacteriaceae", "collinsella")

# Home

header = dashboardHeader(title = "METAVOLIKOS | Explore genes and microbiota in metabolic diseases", titleWidth = 650)
sidebar = dashboardSidebar(
  width = 300, 
  sidebarMenu(id = "tabs", 
              menuItem("Publication Search", 
                       tabName = "a", 
                       icon = icon("search")), 
              conditionalPanel("input.tabs == 'a'", 
                               selectInput(inputId = "busqueda", 
                                           label = "Metabolic disease:", 
                                           choices = choices, 
                                           selected = NULL, 
                                           multiple = FALSE, 
                                           selectize = TRUE, 
                                           width = NULL, 
                                           size = NULL), 
                               dateInput(inputId = "fechadesde", 
                                         "Since:", 
                                         value = "2019-01-01", 
                                         format = "yyyy-mm-dd"), 
                               dateInput("fechahasta", 
                                         "Until:", 
                                         value = "2021-01-01", 
                                         format = "yyyy-mm-dd"), 
                               tags$div(id = "espacio1", br()), 
                               actionButton(inputId = "downPub1", 
                                            label = "Apply Search", 
                                            color = "success", 
                                            style = "gradient", 
                                            block = FALSE, 
                                            size = "md"),
                               tags$div(id = "espacio2", br()),
                               actionButton(inputId = "downPub2", 
                                            label = "Download Information", 
                                            color = "success", 
                                            style = "gradient", 
                                            block = FALSE, 
                                            size = "md"), 
                               tags$div(id = "espacio3", br()), 
                               hidden(
                                 actionButton(inputId = "nuevaBusq", 
                                              label = "New Search", 
                                              color = "success", 
                                              style = "gradient", 
                                              block = FALSE, 
                                              size = "md"))),
              menuItem("Explore Publications", 
                       tabName = "b", 
                       icon = icon("chart-bar")), 
              conditionalPanel("input.tabs == 'b'"),
              menuItem("Gene-Disease Relation", 
                       tabName = "c", 
                       icon = icon("dna")),
              conditionalPanel("input.tabs == 'C'"),
              menuItem("Microbiota-Disease Relation", 
                       tabName = "d", 
                       icon = icon("bacterium")),
              conditionalPanel("input.tabs == 'd'")))

# Body tabs

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
        box(id = "boxdebusqueda",
            title = "Applied search",
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            width = 11,
            textOutput("nText"),
            htmlOutput("nText2")),
        box(id = "resultado",
            title = "Search result", 
            solidHeader = TRUE, 
            status = "success", 
            collapsible = TRUE, 
            width = 11,
            tags$div(style = 'cursor:pointer', 
                     dataTableOutput("pmidresult")),
            br(),
            wellPanel(htmlOutput("rows_selected"))
        )
      ),
      tabItem(
        tabName = "b", h2("Genes and word clouds", align = "center"),
        br(),
        fluidRow(         
          column(6,
                 box(id= "Word Cloud", 
                     title = "Word Cloud", 
                     width = 11,  
                     solidHeader = TRUE, 
                     status = "success", 
                     plotOutput("nubepalabras"))),
          column(6,
                 box(id = "Gene Cloud", 
                     title = "Gene Cloud", 
                     width = 11, 
                     solidHeader = TRUE, 
                     status = "success", 
                     plotOutput("nubegenes")))),
        fluidRow(      
          column(6,br(),
                 box(id= "Wordtable", 
                     title = "Word Table", 
                     width = 11,  
                     solidHeader = TRUE, 
                     status = "success", 
                     dataTableOutput("tablapalabras"),
                     collapsible = TRUE,
                     collapsed = TRUE)),
          column(6,br(),  
                 box(id= "Genetable", 
                     title = "Gene Table", 
                     width = 11,  
                     solidHeader = TRUE, 
                     status = "success", 
                     dataTableOutput("tablagenes"),
                     collapsible = TRUE,
                     collapsed = TRUE))),
        fluidRow( 
          column(6, br(),
                 box(id="parametroswordcloud", 
                     title = "Parameters Words", 
                     width = 11, 
                     solidHeader = TRUE, 
                     status = "primary", 
                     sliderInput("minpalabras", 
                                 label = "Word frequency Min.:", 
                                 min = 1, 
                                 max = 15, 
                                 value = 1),
                     sliderInput("maxpalabras", 
                                 label = "Word frequency Max.:", 
                                 min = 1, 
                                 max = 500, 
                                 value = 100))),
          column(6, br(),
                 box(id="parametrosgenecloud", 
                     title = "Parameters Genes", 
                     width = 11, 
                     solidHeader = TRUE, 
                     status = "primary", 
                     sliderInput("mingenes", 
                                 label = "Gene frequency Min.:", 
                                 min = 1, 
                                 max = 15, 
                                 value = 1),
                     sliderInput("maxgenes", 
                                 label = "Gene frequency Max.:", 
                                 min = 1, 
                                 max = 500, 
                                 value = 100))))
      ),
      tabItem(
        tabName = "c", h2("Relation Genes - Obesity", align = "center"),
        fluidRow(br(),
                 column(4,
                        box(id = "selectMD", 
                            title = "Parameters", 
                            solidHeader = TRUE, 
                            status = "primary", 
                            width = 11, 
                            dateInput(inputId = "fechadesde2", 
                                      "Since:", 
                                      value = "2019-01-01", 
                                      format = "yyyy-mm-dd"),
                            dateInput(inputId = "fechahasta2", 
                                      "Until:", 
                                      value = "2021-01-01", 
                                      format = "yyyy-mm-dd"),
                            actionButton(inputId = "downPub3", 
                                         label = "Search", 
                                         color="success", 
                                         style = "gradient",
                                         block = FALSE),
                            tipify(actionButton(inputId = "Instrucciones", 
                                                label = "", 
                                                color = "success", 
                                                style = "gradient", 
                                                icon = icon("question-circle")),
                                   "A bibliographic search of genes related to metabolic diseases is carried out on PubMed. The search uses the MeSH tag: Metabolic Diseases. The result specifies those genes that have been found under the MeSH tag: Obesity.", trigger = "click")
                                                    )),
                 column(7,
                        box(id="relacionesgenes3D", 
                            title = "Obesity-Gene Relation Graph", 
                            solidHeader = TRUE, 
                            status = "success", 
                            width = 30
                            , 
                            rglwidgetOutput("relacionesgenes3D")
                            ))),
        fluidRow(br(),
                 column(12, 
                        box(id = "Related Genes Table",
                            title = "Related Gene Table", 
                            solidHeader = TRUE, 
                            status = "success", 
                            width = 11, 
                            tags$div(style='cursor:pointer', 
                                     dataTableOutput("pmidresult2"), br(), 
                                     hidden(textOutput("nTextObesity_new")), 
                                     hidden(textOutput("genes_warn"))))))
        ),
    tabItem(
      tabName = "d", h2("Relation Microbiota - Obesity", align = "center"),
      fluidRow(br(),
               column(4,
                      box(id = "selectMD2", 
                          title = "Parameters", 
                          solidHeader = TRUE, 
                          status = "primary", 
                          width = 11, 
                          dateInput(inputId = "fechadesde3", 
                                    "Since:", 
                                    value = "2019-01-01", 
                                    format = "yyyy-mm-dd"),
                          dateInput(inputId = "fechahasta3", 
                                    "Until:", 
                                    value = "2021-01-01", 
                                    format = "yyyy-mm-dd"),
                          actionButton(inputId = "downPub4", 
                                       label = "Search", 
                                       color="success", 
                                       style = "gradient",
                                       block = FALSE),
                          tipify(actionButton(inputId = "Instrucciones2", 
                                              label = "", 
                                              color = "success", 
                                              style = "gradient", 
                                              icon = icon("question-circle")),
                                 "A bibliographic search of microbiota related to metabolic diseases is carried out on PubMed. The search uses the MeSH tag: Metabolic Diseases. The result specifies those bacteria population that have been found under the MeSH tag: Obesity. A secondary corpus including the terms microbiota or microbiome is performed.", trigger = "click")
                          )),
               column(7,
                      box(id="relacionesmicrobiota3D", 
                          title = "Obesity-Microbiota Relation Graph", 
                          solidHeader = TRUE, 
                          status = "success", 
                          width = 30
                          , 
                          rglwidgetOutput("relacionesmicrobiota3D")
                      ))),
      fluidRow(br(),
               column(12, 
                      box(id = "Related Microbiota Table",
                          title = "Related Microbiota Table", 
                          solidHeader = TRUE, 
                          status = "success", 
                          width = 11, 
                          tags$div(style='cursor:pointer', 
                                   dataTableOutput("pmidresult3"), br(), 
                                   hidden(textOutput("nTextObesity_new2")), 
                                   hidden(textOutput("genes_warn2"))))))
    )
    )))
    
# User interface
        
ui = dashboardPage(title = "METAVOLIKOS", header, skin = "green", sidebar, body)


# Server

server = function(input, output){
  
  # Publications search
  
  nTexto <- eventReactive(input$downPub1,{paste(c(input$busqueda," [MH] AND ", as.character(input$fechadesde, "%Y/%m/%d"),":", as.character(input$fechahasta, "%Y/%m/%d"), " [DP]"), collapse = "")})
  output$nText = renderText(nTexto())
  pubmedResult = eventReactive(input$downPub2, 
                               {withProgress(message = "Downloading information from PubMed...", detail = "Please, wait :)", value = 0, {incProgress(2/10, detail = "Please, wait :)")
                                 corpus = batch_pubmed_download(pubmed_query_string = paste(c(input$busqueda,"[MH] AND ", as.character(input$fechadesde, "%Y/%m/%d"),
                                                                                              ":", as.character(input$fechahasta, "%Y/%m/%d"), " [DP]"), collapse = ""),
                                                                                              format = "abstract",
                                                                                              batch_size = 500,
                                                                dest_file_prefix = "corpus_",
                                                                dest_dir = paste(getwd()))
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
  output$pmidresult <- renderDataTable({
    corpus <- pubmedResult()
  pmidAbs = data.frame(abtracts=corpus@Abstract)
  titulos = c()
  for(i in 1:length(corpus@PMID)){
    abstractnro = i
    text = unlist(strsplit(corpus@Abstract[abstractnro], "\\."))
    op = which(unlist(lapply(unlist(strsplit(corpus@Abstract[abstractnro], "\\.")), wordcount)) >=6)[1]
    titulos = c(titulos, text[op])
  }
  Intro = unlist(lapply(pmidAbs$abtracts, substr, start=1, stop=300))
  pmidRes = data.frame(PMID = corpus@PMID, Title = titulos)
  datatable(pmidRes, 
            selection=list(mode = 'single', selected = 1), 
            options = list(lenghtMenu = c(5,30,50), pageLenght=5)) %>% formatStyle(names(pmidRes), 
                                                                                   cursor = "hand", 
                                                                                   target = "cell")
  
  })                                                               
  observeEvent(input$downPub2,{
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
  
  observeEvent(input$pmidresult_rows_selected, {
   output$rows_selected <- renderPrint({
    corpus = pubmedResult()
    text = unlist(strsplit(corpus@Abstract[input$pmidresult_rows_selected],"\\."))
    op = which(unlist(lapply(unlist(strsplit(corpus@Abstract[input$pmidresult_rows_selected],"\\.")), wordcount))>=6)[1]
    titulos = text[op]
    for (i in 1:length(text)){
      if(i == op)cat(paste('<p><h4>','<font color=\"#0c820e\"><b>',text[i], '</b></font>', '</h4>', '\n', '<p>'), fill = TRUE)
      else cat(paste('<SPAN class=sentence><i>', text[i],".", '</i></SPAN>'), fill = TRUE)
    }
    cat(paste("<p><a href='https://pubmed.ncbi.nlm.nih.gov/",corpus@PMID[input$pmidresult_rows_selected],"/","'target=_blank>","Link to PubMed publication","</a><p>", sep = ""))
    })})
  
  # Explore publications
  
  output$nubepalabras = renderPlot(width = "auto", height = "auto", {
    withProgress(message = "Processing corpus...", detail = "Please, wait :)",
                 value = 0, {incProgress(2/10, detail = "Please, wait :)")
                   corpus = pubmedResult()
                   words = word_atomizations(corpus)
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
    wordcloudRep(words = freqGenes$Gene_symbol, freq = as.numeric(freqGenes$Freq), max.words = input$maxgenes, min.freq = input$mingenes, random.color = TRUE, colors = brewer.pal(8, "Dark2"), scale = c(8,0.6))
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
 
  # Relation Gene-Disease
  
  pubmedResult2 <- 
    eventReactive(input$downPub3, 
                                   {withProgress(message = "Downloading information from PubMed...", detail = "Please, wait :)",
                                              value = 0, {incProgress(0.1, detail = "Please, wait :)")
                                                corpusMD = batch_pubmed_download(pubmed_query_string = 
                                                                                   paste(c(
                                                                                     "Metabolic disease [MH] AND ", 
                                                                                                               as.character(input$fechadesde2, "%Y/%m/%d"),
                                                                                                               ":", 
                                                                                                               as.character(input$fechahasta2, "%Y/%m/%d"),
                                                                                                               " [DP]"), collapse = ""),
                                                                                 format = "abstract", 
                                                                                 batch_size = 500, 
                                                                                 dest_file_prefix = "corpusMD_", 
                                                                                 dest_dir = paste(getwd()))
                                                
 incProgress(0.1, detail = "Consolidating files...")
 file.create("pubmed_result2.txt")
 for (i in 1:length(corpusMD)){
 file.append("pubmed_result2.txt", corpusMD[i])}
 corpus_query_comb = readabs("pubmed_result2.txt")
  
 incProgress(0.2, detail = "Extracting genes...")
 genes_MD_comb = data.frame(gene_atomization(corpus_query_comb))
 genes_MD_comb$Freq = as.numeric(genes_MD_comb$Freq)
 genes_MD_list_comb = genes_MD_comb[genes_MD_comb$Freq>10,1]
 
 terms_obesity = c("bmi", "metabolic", "diabetes", "glucose", "liver", "adipose")
 incProgress(0.1, detail = "LSA")
 tdm_main = tdm_for_lsa(corpus_query_comb, c(genes_MD_list_comb, "obesity", terms_obesity, Disease_one_word_without_Obesity))
 lsa_corpus = lsa(tdm_main, dims = dimcalc_share())
 matriz = as.textmatrix(lsa_corpus)
 aa = associate(matriz, "obesity", measure = "cosine", threshold = 1e-399)
 coseno_out = data.frame(genes = names(aa), coseno= aa, stringsAsFactors = FALSE)
 incProgress(0.1, detail = "Obesity publications...")
 output_obesity = batch_pubmed_download(pubmed_query_string = 
                                          paste(c(
                                            "Obesity [MH] AND ", 
                                        as.character(input$fechadesde2, "%Y/%m/%d"),
                                                                      ":", 
                                        as.character(input$fechahasta2, "%Y/%m/%d"), "[DP]"), collapse = ""),
                                                                    format = "abstract", 
                                        batch_size = 500, 
                                        dest_file_prefix = "pubmedobesity",
                                        dest_dir = paste(getwd()))
 
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
 funcion_relMD = function(x){combine_words(names(associate(matriz, x, measure = "cosine", threshold = 0.4))[which(names(associate(matriz,x, measure = "cosine", threshold = 0.4)) %in% Disease_one_word_without_Obesity)], sep = ",", and = "and")}
 result_genes$other_MD = lapply(result_genes$genes, funcion_relMD)
 result_genes[result_genes$other_MD == "character(0)", 4] = NA
 result_genes$coseno = round(result_genes$coseno, 5)
 incProgress(0.3, detail = "Complete")
 row.names(result_genes) = seq(1:length(result_genes$genes))
 mylist = list("matriz" = matriz, "result_genes" = result_genes, "terms" = terms_obesity)
 return(mylist)          
 
                                              })})
 
 output$pmidresult2 = 
   renderDataTable({
   result_genes = pubmedResult2()$result_genes
   datatable(result_genes, colnames = c("Gene", "Cosine similarity with obesity", "Mentioned in obesity publications", "Other MD related"),
             selection = list(mode='single', selected = 1), options = list(lenghtMenu = c(5,30,50), pageLength = 5, order = list(list(1, 'desc'))),
             rownames = FALSE) %>% formatStyle(names(result_genes), cursor ="hand", target = "cell")%>% formatStyle(columns = 4, fontSize = '80%')%>% formatStyle(c("genes", "coseno", "obesity", "other_MD"), textAlign = 'center')
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
 

 # Relation Microbiota - Disease
 
 pubmedResult3 <- 
   eventReactive(input$downPub4, 
    {withProgress(message = "Downloading information from PubMed...", detail = "Please, wait :)",
    value = 0, {incProgress(0.1, detail = "Please, wait :)")
   corpusMD2 = batch_pubmed_download(pubmed_query_string = 
                                       paste(c(
                                         "Metabolic disease [MH] AND ", 
                                     as.character(input$fechadesde3, "%Y/%m/%d"),
                                     ":", 
                                     as.character(input$fechahasta3, "%Y/%m/%d"),
                                     " [DP]"), collapse = ""),
                                     format = "abstract", 
                                     batch_size = 500, 
                                     dest_file_prefix = "corpusMD2_", 
                                     dest_dir = paste(getwd()))
 
 incProgress(0.1, detail = "Consolidating files...")
 file.create("pubmed_result3.txt")
 for (i in 1:length(corpusMD2)){
   file.append("pubmed_result3.txt", corpusMD2[i])}
 corpus_query_comb2 = readabs("pubmed_result3.txt")
 corpus_query_comb2 = searchabsL(corpus_query_comb2, include = "microbiota | microbiome")
 
 incProgress(0.2, detail = "Extracting microbiota...")
 
 micro_MD_comb = data.frame(word_atomizations(corpus_query_comb2))

 colour_match <- str_c(Microbiota_one_word, collapse = "|")
 has_colour <- str_subset(micro_MD_comb$words, colour_match)
 matches <- str_extract(has_colour, colour_match)
 micro_MD_list_comb <- unique(matches)

 terms_obesity = c("bmi", "metabolic", "diabetes", "glucose", "liver", "adipose")
 incProgress(0.1, detail = "LSA")
 tdm_main2 = tdm_for_lsa(corpus_query_comb2, c(micro_MD_list_comb, "obesity", terms_obesity, Disease_one_word_without_Obesity))
 lsa_corpus2 = lsa(tdm_main2, dims = dimcalc_share())
 matriz2 = as.textmatrix(lsa_corpus2)
  aa2 = associate(matriz2, "obesity", measure = "cosine", threshold = 1e-399)
 coseno_out2 = data.frame(microbiota = names(aa2), coseno= aa2, stringsAsFactors = FALSE)
 incProgress(0.1, detail = "Obesity publications...")
 output_obesity2 = batch_pubmed_download(pubmed_query_string = 
                                           paste(c(
                                             "Obesity [MH] AND ",
                                         as.character(input$fechadesde3, "%Y/%m/%d"),
                                         ":", 
                                         as.character(input$fechahasta3, "%Y/%m/%d"), "[DP]"), collapse = ""),
                                         format = "abstract", 
                                         batch_size = 500, 
                                         dest_file_prefix = "pubmedobesity2",
                                         dest_dir = paste(getwd()))
 
 file.create("pubmed_result_obesity2.txt")
 for(i in 1:length(output_obesity2)){
   file.append("pubmed_result_obesity2.txt", output_obesity2[i])}
 corpus_obesity2 = readabs("pubmed_result_obesity2.txt")
 corpus_obesity2 = searchabsL(corpus_obesity2, include = "microbiota | microbiome")
 incProgress(0.1, detail = "Obesity-associated microbiota")
 
 micro_MD_comb2 = data.frame(word_atomizations(corpus_obesity2))
 colour_match2 <- str_c(Microbiota_one_word, collapse = "|")
 has_colour2 <- str_subset(micro_MD_comb2$words, colour_match2)
 matches2 <- str_extract(has_colour2, colour_match2)
 micro_MD_list_comb2 <- unique(matches2)

 micro_MD_list_obesity = data.frame(microbiota = micro_MD_list_comb2, obesity = "yes", stringsAsFactors = FALSE)
 
 result_micro = merge(coseno_out2, micro_MD_list_obesity, by.x = "microbiota", by.y = "microbiota", all.x = TRUE)
 result_micro = result_micro[!result_micro$microbiota %in% c(terms_obesity, Disease_one_word_without_Obesity),]
 funcion_relMD2 = function(x){combine_words(names(associate(matriz2, x, measure = "cosine", threshold = 0.4))[which(names(associate(matriz2,x, measure = "cosine", threshold = 0.4)) %in% Disease_one_word_without_Obesity)], sep = ",", and = "and")}
 
 result_micro$other_MD = lapply(result_micro$microbiota, funcion_relMD2)
 
 result_micro[result_micro$other_MD == "character(0)", 4] = NA
 result_micro$coseno = round(result_micro$coseno, 5)
 incProgress(0.3, detail = "Complete")
 row.names(result_micro) = seq(1:length(result_micro$microbiota))
 mylist2 = list("matriz" = matriz2, "result_micro" = result_micro, "terms" = terms_obesity)
 return(mylist2)                            
 
    })})
 
 output$pmidresult3 = 
   renderDataTable({
     result_micro = pubmedResult3()$result_micro
 tabla2 = datatable(result_micro, colnames = c("Microbiota", "Cosine similarity with obesity", "Mentioned in obesity publications", "Other MD related"),
                    selection = list(mode='single', selected = 1), options = list(lenghtMenu = c(5,30,50), pageLength = 5, order = list(list(1, 'desc'))),
                    rownames = FALSE) %>% formatStyle(names(result_micro), cursor ="hand", target = "cell")%>% formatStyle(columns = 4, fontSize = '80%')%>% formatStyle(c("microbiota", "coseno", "obesity", "other_MD"), textAlign = 'center')
   })  

 output$relacionesmicrobiota3D = renderRglwidget({
 matriz2 = pubmedResult3()$matriz
 terms_obesity = pubmedResult3()$terms
 termMatrix3 = matriz2%*%t(matriz2)
 exclude2 = which(row.names(termMatrix3) %in% terms_obesity)
 include2 = which(row.names(termMatrix3)%in% c((pubmedResult3()$result_micro)$micro, "obesity"))
 plot_neighbors("obesity", n=length(include2), tvectors = matriz2[include2,include2], legend = T, connect.lines = 0, alpha = c(0.5,0.5), col = c("lemonchiffon", "orange", "darkred"))
 rglwidget()
 })

}
                                        
#Shiny app  
  
shinyApp(ui, server)



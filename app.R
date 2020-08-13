# genes do not always show on regional plot when they overlap the border. 
print(Sys.getenv("USER"))
# # LOCAL:
if (Sys.getenv("USER")=="niek"){
  setwd("/Users/niek/repos/ecgenetics/")
} else if (Sys.getenv("USER")=="benzhi") {
  setwd("/mnt/linux_data/Repos/bitbucket_itc/itc-main-repo/Niek/ECG/exerciseecg_data/ECG_averaging/shinyapp/ecgenetics7")
} else{
  setwd("/srv/shiny-server/") # <- for strato
}
source("global.R")
source("helpers.R")
source("helpers.heatmap.R")
source("helpers.regionalplot.R")
source("helpers.mendelianrandomization.R")

source("get_nearest_gene.r")
load("data/tsnedata.Rdata")

library(DT)
library(shinyWidgets)
####

ui <- navbarPage(title="ECGenetics Browser",
                 # footer =  div(
                 #     id = "img-id",
                 #     img(src = "https://www.opciweb.nl/wp-content/uploads/2017/08/umcg-rug-onderzoek-ci-leren-engels-horizontaal.gif",width="10%")
                 #   ),
   tabPanel("ECG Plotter",
            titlePanel(h3("ECG Plotter")),
            hr(),includeCSS("style.css"),
      sidebarLayout(
        sidebarPanel(
          width=3,
          
#          hr(),
          p("Enter a list of SNPs (rs572474770, rs776293589, 2:179698596), a region (2:179381323-179405807) or a Gene (e.g. SCN5A) and click the Go button to extract variants."), #, Region (e.g. 10:30000-40000) or Gene name (e.g. TTN) 
          textAreaInput("rsid", "RSids, region or gene",value = paste(dftsne$rsid,collapse = ", "),#TTN independent
                        width="100%",rows = 4),
          radioButtons("subset", "Search Dataset",
                       c("Tophits (140k+, fast)" = "tophits",
                         "GW (19M+ slow)" = "all")),
          radioButtons("phenotype", "PhenoType",
                       c("Unadjusted" = "unadjusted",
                         "RR-adjusted" = "stretch")),
          sliderInput("slider_window","Extend region (Kb) in case you provide a single region or gene.",
                      min=0,max=250, value=50, step = 10),
          checkboxInput(inputId="include_betas", label="Include the coefficients and 95% CI!", value = FALSE, width = NULL),
          useShinyjs(),
          shinyjs::hidden(
            checkboxInput(inputId="plot_adjusted_means", label="> plot the predicted effect on the means. This is intented to be an illustration, as the effect is exagragated.", value = TRUE, width = NULL)
          ),
          actionButton("goButton", "Go!"), #downloadButton("go_downloadall", "SNP Data"),
          hr(),
          img(src = "https://www.opciweb.nl/wp-content/uploads/2017/08/umcg-rug-onderzoek-ci-leren-engels-horizontaal.gif",width="100%"),
          p("under development, n.verweij@umcg.nl. 2019" ,align = "right", style = "font-size:10px;  font-style: italic ")
          #p("n.verweij@umcg.nl", style = "font-size:11px")
        ),
        mainPanel(
          plotlyOutput("oPlot",width="100%",height="auto"),
          tabsetPanel(id="analysistabs",
            tabPanel("Table",
              DT::dataTableOutput("oTable"),
              downloadButton("downloadData_table", "Download"),
              verbatimTextOutput("oText")
            ),
            tabPanel("Regional Plot",
                     actionButton("regionalprefreshButton", "Refresh"),
                     plotlyOutput("regionalPlot",width="100%",height="auto"),
                     plotlyOutput("regionalPlot_genes",width="100%",height="125"),
                     p("Colors indicate LD (r2); filled circles indicate the significance of the lead ECG trait; hollow circles indicate the minimum value across ECG all traits.")
                     ),
            tabPanel("Heatmap", 
                     div(style="display: inline-block;vertical-align:top; width: 200px;", ## divs for inline. 
                     shinyWidgets::sliderTextInput("slider_hmpvalue_treshold","P-Value treshold:",
                                                   choices=c(0.00000000001,0.0000000001,0.000000001,0.00000001,0.00000005,0.0000001,0.0000005 ,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005, 0.01,0.05,1),
                                                   selected=0.00000005, grid = T)),
                     div(style="display: inline-block;vertical-align:top; padding-top: 27px; width: 125px;",  verbatimTextOutput("oTextHeatmap")),
                     div(style="display: inline-block;vertical-align:top; padding-top: 27px; width: 150px;",  actionButton("hmrefreshButton", "Refresh")),   # https://github.com/dreamRs/shinyWidgets/issues/91
                     p("Please limit to <150 for now, else it will probably fail to plot due to memory issues."),
                      
                     plotlyOutput("hmPlot",width="100%",height="auto")
            )
          )
        )
      )
    ),
   # tabPanel("MATTERHORN",
   #          titlePanel(h3("MATTERHORN: THE REGION GENErator")),
   #          hr(),
   #          p("v0.00 - Under construction.  2019" ,align = "left"),
   #          verbatimTextOutput("Summary"),
   #          img(src='https://db-service.toubiz.de/var/plain_site/storage/images/orte/zermatt/matterhorn/cr-pascal-gertschen-2018-zt_sommer_web_dscf0814/3206064-1-ger-DE/cr-Pascal-Gertschen-2018-ZT_Sommer_Web_DSCF0814_front_large.jpg', align = "left")
   # ),
   tabPanel("Interactive t-SNE plot",
            titlePanel(h3("Interactive t-SNE plot (manuscript)")),
            plotlyOutput("oPlotTsne",width="100%",height="auto"),
            splitLayout(width = "1005px" , 
                        cellArgs = list(style =  "width: 500px"), #cellArgs = list(style =  "width: 500px, float:left; display:inline"),
                        plotlyOutput("oPlot_ecg_unadjusted"), 
                        plotlyOutput("oPlot_ecg_stretch")),
            
            # plotlyOutput("oPlot_ecg_unadjusted",width="100%",height="auto"),
            # plotlyOutput("oPlot_ecg_stretch",width="100%",height="auto"),
            
            hr()
   ),
  tabPanel("ECG-wide Mendelian randomization",
           titlePanel(h3("ECG-wide Mendelian randomization")),
           
           fileInput('fmrexposures', 'Choose file. '),
           "Provide a summary statistics file with up to 1000 variants. Include the following columns with headers: CHR, BP,  EFAL, NEFAL, EAF, BETA, SE. Ensure positions are on Build 37. This is a try-out. It may be very slow since it needs to query the full 75 GB gz-compressed data and the VPS has no superpowers.",
           tags$a(href="https://raw.githubusercontent.com/niekverw/ecgenetics/master/data/QTinterval.tsv", "Download an example (QT interval) here. "),
           # splitLayout(width = "1005px" , 
           #             cellArgs = list(style =  "width: 500px"), #cellArgs = list(style =  "width: 500px, float:left; display:inline"),
           #             plotlyOutput("omrPlot_ecg_unadjusted"), 
           #             plotlyOutput("omrPlot_ecg_stretch")),
           splitLayout(width = "1005px" , 
                       cellArgs = list(style =  "width: 500px"), #cellArgs = list(style =  "width: 500px, float:left; display:inline"),
                       plotlyOutput("Omrplot_unadjusted"), 
                       plotlyOutput("Omrplot_rradjusted")),
           textOutput("txt_num_mr_variants"),
           plotOutput("Oplot_freqcheck",width=400,height=400),
           tableOutput('tbl_dfmrexposures'),
           hr()
  ),
  tabPanel("Terms and Data Information",
           titlePanel(h3("Terms and Data Information")),
           hr(),
           verbatimTextOutput("About"),
           includeMarkdown("about.md")
  )

#     tabPanel("Downloads",
#          titlePanel(h3("Data")),
#          hr(),
#          includeMarkdown("download.md")
# )

)

server <- function(input, output, session) {
  hideTab(inputId = "analysistabs", target = "Regional Plot")
  rsids=c("rs201106462")
  makeReactiveBinding("df.snpinfo")
  futureData <- reactiveValues(data=NULL)

  ###### first input
  observeEvent(input$goButton,  ignoreInit=TRUE, {
    #print(input$rsid)
    #showTab(inputId = "analysistabs", target = "Table") # oPlot panelid=analysistabs , tabpanel=Table
    updateTabsetPanel(session, "analysistabs",
                      selected = "Table"
    )
    shinyjs::hide("plot_adjusted_means")
    
    # load preloaded data quicly (all lead snps) if nothing on input is changed
    if(input$rsid == paste(dftsne$rsid,collapse = ", ")){
      if (input$phenotype =="unadjusted" & !input$include_betas ){
        futureData$data <<- datatsne.unadjusted
        return(NULL)
      }
      if (input$phenotype =="stretch" & !input$include_betas){
        futureData$data <<- datatsne.stretch
        return(NULL)
      }
      
    }
    # record query settings
    query <- process_user_input(input$rsid,
                                mapping.proteincoding,
                                window=input$slider_window*1000,
                                subset=input$subset,
                                phenotype=input$phenotype)
    cat(file=stderr(),"QUERY: ",session$clientData$url_hostname,paste(";",paste(unlist(query),collapse=";"),"\n"))

    if(query$error){shiny::showNotification("invalid input", type = "error")} # catch error
    validate(need(!query$error,"error with query.")) # catch error
    
    tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)
    if(tabix_query$error){shiny::showNotification("invalid query", type = "error")} # catch error
    validate(need(!tabix_query$error,"error with query.")) # catch error
    
    if (query$entry == "positional"){
      showTab(inputId = "analysistabs", target = "Regional Plot")
    } else {
      hideTab(inputId = "analysistabs", target = "Regional Plot")
    }

    #print(df.static_phenotype)
    f.data_beta=""
    f.data_se=""
    if( query$subset == "tophits"){
      f.data_p=paste0(datadir,"/tophits_data/",query$phenotype,".logP.outfile.tsv.gz.tophits.gz")
      f.data.index=paste0(datadir,"/tophits_data/","unadjusted.logP.outfile.index.tsv.gz.tophits.gz")
      if(input$include_betas){
        f.data_beta=paste0(datadir,"/tophits_data/",query$phenotype,".BETA.outfile.tsv.gz.tophits.gz")
        f.data_se=paste0(datadir,"/tophits_data/",query$phenotype,".SE.outfile.tsv.gz.tophits.gz")
        shinyjs::show("plot_adjusted_means")
      }
    } else {
      f.data_p = paste0(datadir,"/full_data_combined/",query$phenotype,".logP.outfile.tsv.gz")
      f.data.index = paste0(datadir,"/full_data_combined/","unadjusted.logP.outfile.index.tsv.gz")
      if(input$include_betas){
        f.data_beta=paste0(datadir,"/full_data_combined/",query$phenotype,".BETA.outfile.tsv.gz")
        f.data_se=paste0(datadir,"/full_data_combined/",query$phenotype,".SE.outfile.tsv.gz")
        shinyjs::show("plot_adjusted_means")
      } 
    }
    showModal(modalDialog("Please wait.", footer=NULL))
    
    myFuture <- future(  {
      data <- extract_multiple_variants(tabix_query,
                                        f.data_p,
                                        f.data_beta,f.data_se,
                                        f.data.index)
      return(data)
    },  globals = list(user_input = input$rsid,
                       mapping.proteincoding= mapping.proteincoding,
                       tabix_query=tabix_query,
                       f.data_p =f.data_p,f.data_beta=f.data_beta,f.data_se=f.data_se,
                       f.data.index = f.data.index))
    
    promises::then(
      myFuture,
      onFulfilled = function(value) {
        futureData$data <<- value
      },
      onRejected = function(value) {
        shiny::showNotification("No data found", type = "error")
        removeModal()
      }
      )
    
    #### REPLACE ABOVE WITH THIS ON MAC???
    # futureData$data <<- extract_multiple_variants(tabix_query,
    #                                               dir_data,
    #                                               f.data_p,
    #                                               f.data_beta,f.data_se,
    #                                               f.data.index)
      return(NULL)
    }) ## end observe event

  
    #https://yihui.shinyapps.io/DT-proxy/
    ##### LOAD TABLE 
    output$oTable <- DT::renderDataTable({
        req(futureData$data)
        removeModal()
        df=futureData$data$df_snp_info
        
        df=df[,c("SNP", "CHR", "BP", "EFAL", "NEFAL", "EAF",
                                 "INFO", "minP","band","Gene")] 
        colnames(df) = c("SNP", "Chr", "Pos", "Effect Allele", "Non-Effect Allele", "Freq",
                           "Info", "minP","band","Gene") 
        df #<- datatable(df) %>% formatRound(c('Info','Freq',"maxLogP"),digits = 3)  <- selection=.. options are not working when this is added.  
        
        
      },
        options = list(order = list(list(7, 'asc'))),
        selection=list(mode='single',  selected=which.min(futureData$data$df_snp_info$minP)), # TODO: replace with most significant index from the intial query. 
        rownames= FALSE,filter = 'bottom') #single or multiple. 
   
    
    #### TEXT UNDER TABLE, DEBUGGING
    output$oText <- renderText({
			  req(futureData$data)
        
        shiny::validate(need(input$oTable_rows_selected, 'Please select SNP from the table - or wait a bit')) # add multiple by , 
        Clicked <- input$oTable_rows_selected
			  paste0("SNP: ",futureData$data$df_snp_info[Clicked,"SNP"],
			         "\n","NEARBY GENES:",futureData$data$df_snp_info[Clicked,"Gene"],
 			                                "\n","#snps found: ",futureData$data$num_snps_found,
 			                                "\n","#snps missing: ",futureData$data$n.missingsnps,
                                      "\n","---- DEBUGGING",
                                      "\n","time query: user(",round(unname(futureData$data$query_time_tabix[1]),3),"), system(",round(unname(futureData$data$query_time_tabix[2]),3),") elapsed(",round(unname(futureData$data$query_time_tabix[3]),3),")",
			                                "\n","w: ",input$slider_window*1000,  #, dput(futureData$data$df_snp_p[Clicked(),])
			                                "\n",futureData$data$query$chr,":",futureData$data$query$startpos,"-",futureData$data$query$endpos
                                     )
    })
  
    #
    
    #### CAPTURE CLICK EVENTS
    Clicked <- reactiveValues(i=1,BP="NA") # initialize
    # click on table
    observeEvent(input$oTable_rows_selected, {
      Clicked$i <- input$oTable_rows_selected
      Clicked$BP <- "NA"
    })
    # click on regionplot; this is generating a warning on initialization because The 'plotly_click' event tied a source ID of 'source_rp' is not registered.'
    
    observeEvent( event_data(event = "plotly_click",source = "source_rp"), {
      req(futureData$data)
      s <- event_data(event = "plotly_click",source = "source_rp")
      Clicked$i <- "NA"
      Clicked$BP <- s$x
    })
    # click on heatmap  (( called source="A" because "A" is plotly's default and can't change it. ))
    observeEvent(event_data(event = "plotly_click",source = "A"), {
      req(futureData$data)
      s <- event_data(event = "plotly_click",source = "A")
      Clicked$i <- rev(hmPlot()$order)[s$y] 
      Clicked$BP <- "NA"
      
    })
    
    #### ECG plot 
    output$oPlot <- renderPlotly({
    	 req(futureData$data)
       if(Clicked$BP!="NA"){ # if clicked on region plot, we don't know the index
         print("Clicked$BP")
         print(Clicked$BP)
         Clicked$i <- which(futureData$data$df_snp_info$BP == Clicked$BP)
       # } else { # else...?
       #   print("Clicked$i")
       #   print(Clicked$i)
       }
       if (futureData$data$query$phenotype =="unadjusted") {
         df_ecg_stats = df_ecg_unadjusted
       }
       if (futureData$data$query$phenotype =="stretch") { 
         df_ecg_stats = df_ecg_stretch
       }
       ecg_plot <- suppressWarnings(  make_ecg_plot(vct_snp_p=futureData$data$df_snp_p[Clicked$i,],
                                 vct_snp_beta=if(futureData$data$df_snp_beta !=""){futureData$data$df_snp_beta[Clicked$i,]}else{""}, # not used//todo
                                 vct_snp_se=if(futureData$data$df_snp_se !=""){futureData$data$df_snp_se[Clicked$i,]}else{""}, # not used//todo
                                 vct_snp_info=futureData$data$df_snp_info[Clicked$i,],
                                 df_ecg_stats = df_ecg_stats,
                                 invert=FALSE,plot_adjusted_means = input$plot_adjusted_means) # average ecg signal.
                  )
       ggplotly(ecg_plot,tooltip = "text",height = 400, width = 500,dynamicTicks=TRUE,source="source_ecgplot" )
       
    })
    
    #### REGIONAL PLOT
    regionalPlot <- reactiveVal(NULL)
    
    observeEvent(c(input$regionalprefreshButton), { #input$goButton,
      req(futureData$data)
      showModal(modalDialog("Please wait. LD is currently very slow and might take > 30sec, will make this faster soon.", footer=NULL))
      # showModal(modalDialog(p("Region is Plotting, This might take a while"),
      #                       title = "Region Plotting started"
      # ))
      data = futureData$data
      f <- future({
        #https://www.r-bloggers.com/shiny-1-2-0-plot-caching/
        #regional_plotly(data)
        regional_plot(data,for_plotly = T,LD=TRUE)
        
      },  globals = list(data=data))
      
      f %...>% regionalPlot()
    })
    output$regionalPlot <- renderPlotly({
      removeModal()
      req(regionalPlot())
      ggplotly(regionalPlot()$plt_regional,source="source_rp")
    })
    output$regionalPlot_genes <- renderPlotly({
      removeModal()
      req(regionalPlot())
      ggplotly(regionalPlot()$plt_genes,source="source_rp")
    })
    
    #### HEATMAP
    #rs55851300, rs1051375, rs34081637,rs12541595 , rs10774625, 
    ## https://stackoverflow.com/questions/53079904/shiny-future-error-in-ctxoninvalidate-reactive-context-was-created-in-one
    hm_snps_selected <- reactiveValues(number=0)
    
    slider_hmpvalue_treshold_text <- reactive({
      req(futureData$data)
      hm_snps_selected$number <- sum(as.numeric(futureData$data$df_snp_info$minP<input$slider_hmpvalue_treshold))
      paste("#SNPs: ",hm_snps_selected$number)
    })
    output$oTextHeatmap <- renderText({
      slider_hmpvalue_treshold_text()
    })
    
    hmPlot <- reactiveVal(NULL)
    
    observeEvent(c(input$hmrefreshButton), { #input$goButton,
      req(futureData$data)
      if (hm_snps_selected$number<=1){shiny::showNotification('Need more than 1 variant, provide more variants or change P value treshold', type = "error")}
      validate(need(hm_snps_selected$number>1, 'Need more than 1 variant, provide more variants or change P value treshold'))
      
      if (hm_snps_selected$number>150){shiny::showNotification("too many snps, the max is 150 right now.", type = "error")}
      validate(need(hm_snps_selected$number<=150,"too many snps, max is 150, contact us if more is needed")) # catch error
      
      showModal(modalDialog("Please wait.", footer=NULL))
      data = futureData$data
      pvaltreshold = input$slider_hmpvalue_treshold # <- needs to be here
      f <- future({
        #https://www.r-bloggers.com/shiny-1-2-0-plot-caching/
        make_heatmap_plot(df_snp_p = data$df_snp_p,
                          vct_snp_info = data$df_snp_info,
                          df_ecg_unadjusted, # does nothing yet.
                          pvaltreshold = pvaltreshold)
        
      },  globals = list(df_snp_p = data$df_snp_p,
                         vct_snp_info = data$df_snp_info,
                         df_ecg_unadjusted=df_ecg_unadjusted,
                         pvaltreshold=pvaltreshold))
      
      f %...>% hmPlot()
    })
    output$hmPlot <- renderPlotly({
      #print(length(hmap$x$layout$yaxis$ticktext))
      removeModal()
      hmPlot()
    })



    #################################################
    # Downloadable csv of selected dataset ----
    # #sudo chown -R shiny:shiny-apps data
    #################################################
    output$downloadData_table <- downloadHandler(
      filename = function() {
        paste(gsub("stretch","RRadjusted",futureData$data$query$phenotype), "_ecgenetics_data.zip", sep = "")
      },
      #filename="download.zip",
      content = function(file) {
        phenotypename=gsub("stretch","RRadjusted",futureData$data$query$phenotype)

        req(futureData$data)

        df_snp_info <- futureData$data$df_snp_info
        df_snp_p <- futureData$data$df_snp_p

        names(df_snp_p) <- paste0("pval_",1:500)
        dfecg <- cbind(df_snp_info,df_snp_p)

        if (input$include_betas){
          df_snp_beta <- futureData$data$df_snp_beta
          df_snp_se <- futureData$data$df_snp_se
          names(df_snp_beta) <- paste0("beta_",1:500)
          names(df_snp_se) <- paste0("se_",1:500)

          dfecg <- cbind(dfecg,df_snp_beta,df_snp_se)

        }
        dfecg$phenotype <- futureData$data$query$phenotype
        #write.csv(dfecg, file, row.names = FALSE)
        

        df_averagevolts <- cbind(df_ecg_unadjusted,df_ecg_stretch)
        names(df_averagevolts) <- c(paste0(names(df_ecg_unadjusted),"_unadjusted"),paste0(names(df_ecg_stretch),"_RRadjusted"))

        d.tmp=tempfile(pattern=paste0(phenotypename,"_"),tmpdir="data_tmp/") #paste0(getwd(),"/data"))
        dir.create(d.tmp,showWarnings=TRUE,recursive=TRUE)
        fecg <- paste(d.tmp,"/",phenotypename, ".snpstats.csv", sep = "")
        faveragevoltage <- paste(d.tmp,"/average_voltages.csv", sep = "")
        files <- c(fecg, faveragevoltage)

        future({
         write.csv(dfecg, fecg)
         write.csv(df_averagevolts, faveragevoltage)
         zip(zipfile = file, files = files)
       },contentType = "application/zip") 

      }
    )

    
  ######### TSNE: 
    
    output$oPlotTsne <- renderPlotly({
      
      ptsne <- ggplot(dftsne, aes(x = tsne1, y = tsne2, colour = clusters,text=snps)) + 
        geom_point(alpha = 0.3,show.legend = FALSE ) + theme_bw() 
      plotly::hide_legend(plotly::ggplotly(ptsne,tooltip="text",height = 600, width = 900))

    })
    
    #### CAPTURE CLICK EVENTS
    Clicked_tsne <- reactiveValues(i="NA",rsid="rs116015634") # initialize
    observeEvent( event_data(event = "plotly_click",source = "A"), {
      s <- event_data(event = "plotly_click",source = "A")
      Clicked_tsne$i <- "NA"
      Clicked_tsne$rsid <- dftsne[round(dftsne$tsne1,5) %in%  round(s$x,5),'rsid']  
      
      #plotlyOutput("oPlot_ecg_unadjusted",width="100%",height="auto"),
      #plotlyOutput("oPlot_ecg_stretch",width="100%",height="auto"),

    })

    output$oPlot_ecg_unadjusted <- renderPlotly({
      if(Clicked_tsne$rsid!="NA"){ # if clicked on region plot, we don't know the index
        print("Clicked_tsne$rsid")
        print(Clicked_tsne$rsid)
        Clicked_tsne$i <- which(datatsne.unadjusted$df_snp_info$SNP == Clicked_tsne$rsid)
        # } else { # else...?
        #   print("Clicked$i")
        #   print(Clicked$i)
      
        df_ecg_stats = df_ecg_unadjusted


      ecg_plot <- suppressWarnings(  make_ecg_plot(vct_snp_p=datatsne.unadjusted$df_snp_p[Clicked_tsne$i,],
                                                   vct_snp_beta=if(datatsne.unadjusted$df_snp_beta !=""){datatsne.unadjusted$df_snp_beta[Clicked_tsne$i,]}else{""}, # not used//todo
                                                   vct_snp_se=if(datatsne.unadjusted$df_snp_se !=""){datatsne.unadjusted$df_snp_se[Clicked_tsne$i,]}else{""}, # not used//todo
                                                   vct_snp_info=datatsne.unadjusted$df_snp_info[Clicked_tsne$i,],
                                                   df_ecg_stats = df_ecg_stats,
                                                   invert=FALSE) + ggtitle(paste0(Clicked_tsne$rsid,"(",datatsne.unadjusted$df_snp_info$Gene[Clicked_tsne$i],") - unadjusted")) # average ecg signal.
      )
      ggplotly(ecg_plot,tooltip = "text",height = 400, width = 500,dynamicTicks=TRUE,source="source_ecgplot" )
      }
    })
    
    output$oPlot_ecg_stretch <- renderPlotly({
      if(Clicked_tsne$rsid!="NA"){ # if clicked on region plot, we don't know the index
        print("Clicked_tsne$rsid")
        print(Clicked_tsne$rsid)
        Clicked_tsne$i <- which(datatsne.stretch$df_snp_info$SNP == Clicked_tsne$rsid)
        # } else { # else...?
        #   print("Clicked$i")
        #   print(Clicked$i)
        
        df_ecg_stats = df_ecg_stretch
        
        
        ecg_plot <- suppressWarnings(  make_ecg_plot(vct_snp_p=datatsne.stretch$df_snp_p[Clicked_tsne$i,],
                                                     vct_snp_beta=if(datatsne.stretch$df_snp_beta !=""){datatsne.stretch$df_snp_beta[Clicked_tsne$i,]}else{""}, # not used//todo
                                                     vct_snp_se=if(datatsne.stretch$df_snp_se !=""){datatsne.stretch$df_snp_se[Clicked_tsne$i,]}else{""}, # not used//todo
                                                     vct_snp_info=datatsne.stretch$df_snp_info[Clicked_tsne$i,],
                                                     df_ecg_stats = df_ecg_stats,
                                                     invert=FALSE) + ggtitle(paste0(Clicked_tsne$rsid,"(",datatsne.stretch$df_snp_info$Gene[Clicked_tsne$i],") - RR-adjusted")) # average ecg signal.
        )
        ggplotly(ecg_plot,tooltip = "text",height = 400, width = 500,dynamicTicks=TRUE,source="source_ecgplot" )
      }
    })
    
    
    #####################################################################
    #### MR: 
    ####################################################################
    rv_mr <- reactiveValues(dfmrexposures = NULL,mrplot_unadjusted=NULL,mrplot_rradjusted=NULL)
    
    observe({ # file upload
      req(input$fmrexposures)
      if (is.null(input$fmrexposures)) {return(NULL)}
      print(input$fmrexposures$datapath)
      dfmrexposures <- data.table::fread(input$fmrexposures$datapath, header=T, data.table = F)#, verbose = T)
      ### insert mendelian_randomization.R ()
      dfmrexposures$uniqid <- make_uniqID(dfmrexposures$CHR,dfmrexposures$BP,dfmrexposures$EFAL,dfmrexposures$NEFAL)

      output$tbl_dfmrexposures <- renderTable({
        dfmrexposures
      })
      rv_mr$dfmrexposures <<- dfmrexposures
    })
    
    observeEvent(rv_mr$dfmrexposures,{
      req(rv_mr$dfmrexposures)
      
      input= unique(rv_mr$dfmrexposures$uniqid)
      query <- process_user_input(input,mapping.proteincoding)
      tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)

      showModal(modalDialog("Please wait.", footer=NULL))
      myFuture <- future(  {
        data_unadjusted <- extract_multiple_variants(tabix_query,
                                                     f.data_p="unadjusted.logP.outfile.tsv.gz",
                                                     f.data_beta="unadjusted.BETA.outfile.tsv.gz",
                                                     f.data_se="unadjusted.SE.outfile.tsv.gz",
                                                     f.data.index="unadjusted.logP.outfile.index.tsv.gz"
        )
        
        data_unadjusted <- harmonizedfs(rv_mr$dfmrexposures,data_unadjusted)
        mr_unadjusted <- ecg_wide_ivw(data_unadjusted)
        data_unadjusted$df_snp_info[1,]$SNP <- "IVW-fixed effect unadjusted"
        data_unadjusted$df_snp_info[1,]$Gene <- "IVW-fixed effect unadjusted"
        mrplot_unadjusted <- make_ecg_plot(vct_snp_p=t(data.frame(mr_unadjusted$p)),
                                           vct_snp_beta=t(data.frame(mr_unadjusted$b)), #data$df_snp_beta,
                                           vct_snp_se=t(data.frame(mr_unadjusted$s)), #data$df_snp_se,
                                           vct_snp_info=data_unadjusted$df_snp_info[1,],
                                           df_ecg_stats=df_ecg_unadjusted,
                                           invert=FALSE) #+ ggtitle(paste0("unadjusted"))
        
        
        
        
        data_rradjusted <- extract_multiple_variants(tabix_query,dir_data,
                                                     f.data_p="stretch.logP.outfile.tsv.gz",
                                                     f.data_beta="stretch.BETA.outfile.tsv.gz",
                                                     f.data_se="stretch.SE.outfile.tsv.gz",
                                                     f.data.index="stretch.logP.outfile.index.tsv.gz"
        )
        
        data_rradjusted <- harmonizedfs(rv_mr$dfmrexposures,data_rradjusted)
        mr_rradjusted <- ecg_wide_ivw(data_rradjusted)
        data_rradjusted$df_snp_info[1,]$SNP <- "IVW-fixed effect RR-adjusted"
        data_rradjusted$df_snp_info[1,]$Gene <- "IVW-fixed effect RR-adjusted"
        
        mrplot_rradjusted <- make_ecg_plot(vct_snp_p=t(data.frame(mr_rradjusted$p)),
                                           vct_snp_beta=t(data.frame(mr_rradjusted$b)), #data$df_snp_beta,
                                           vct_snp_se=t(data.frame(mr_rradjusted$s)), #data$df_snp_se,
                                           vct_snp_info=data_rradjusted$df_snp_info[1,],
                                           df_ecg_stats=df_ecg_unadjusted,
                                           invert=FALSE) #+ ggtitle(paste0("RR-adjusted"))
        
        dfmrexposures_NA <- rv_mr$dfmrexposures[!rv_mr$dfmrexposures$uniqid %in% data_unadjusted$df_snp_info$uniqid,]
        dfmrexposures_nonNA <- rv_mr$dfmrexposures[rv_mr$dfmrexposures$uniqid %in% data_unadjusted$df_snp_info$uniqid,]
        
        return(list(dfmrexposures_NA=dfmrexposures_NA,
                    dfmrexposures_nonNA=dfmrexposures_nonNA,
                    data_unadjusted=data_unadjusted,
                    data_rradjusted=data_rradjusted,
                    mrplot_unadjusted=mrplot_unadjusted,
                    mrplot_rradjusted=mrplot_rradjusted))

      })

      promises::then(
        myFuture,
        onFulfilled = function(value) {
          #print(str(value))
          #print(value)
          #rv_mr <<- value
            
        rv_mr$mrplot_rradjusted <<- value$mrplot_rradjusted
        rv_mr$mrplot_unadjusted <<- value$mrplot_unadjusted
        rv_mr$dfmrexposures_NA <<- value$dfmrexposures_NA
        rv_mr$dfmrexposures_nonNA <<- value$dfmrexposures_nonNA
        rv_mr$data_unadjusted <<- value$data_unadjusted
        rv_mr$data_rradjusted <<- value$data_rradjusted
        
        futureData$data <<- value$data_unadjusted
        removeModal()
          #rv_mr$mrplot_rradjusted <<- value$mrplot_rradjusted
        },
        onRejected = function(value) {
          shiny::showNotification("No data found", type = "error")
          removeModal()
        }
      )
    })
    output$Omrplot_rradjusted <- renderPlotly({
      req(rv_mr$mrplot_rradjusted)
        ggplotly(rv_mr$mrplot_rradjusted,tooltip = "text",height = 400, width = 500,dynamicTicks=TRUE,source="source_ecgplot" )
        
    })

    output$Omrplot_unadjusted <- renderPlotly({
      req(rv_mr$mrplot_unadjusted)
        ggplotly(rv_mr$mrplot_unadjusted,tooltip = "text",height = 400, width = 500,dynamicTicks=TRUE,source="source_ecgplot" )
    })
    # 
    output$Oplot_freqcheck <- renderPlot({
      req(rv_mr$data_unadjusted)
      plot(rv_mr$data_unadjusted$dfexposure$EAF ,
           rv_mr$data_unadjusted$df_snp_info$EAF,
           xlab="EAF (Exposure data)",
           ylab="EAF (ECG data)")
    })



    output$txt_num_mr_variants <- renderText({ 
     
      req(rv_mr$dfmrexposures)
        paste0("number of variants in the MR: ",nrow(rv_mr$dfmrexposures_nonNA),"\n",
             "number of variants not found: ",nrow(rv_mr$dfmrexposures_NA))
      
    })

    
    

    
}

shinyApp(ui, server)

# "\n","Heart rate: <tbd>",
# "\n","Heart rate recovery 10ms: <tbd>",
# "\n","Heart rate recovery 50ms: <tbd>",
# "\n","Heart rate increase: <tbd>",
# "\n","QRS Duration: <tbd>",
# "\n","PR Interval: <tbd>",
# "\n","QT Interval: <tbd>",
# "\n","12-Lead Sum: <tbd>",
# "\n","Sokolow-Lyon: <tbd>",
# "\n","Cornel: <tbd>",
# "\n","ST-Wave Voltage: <tbd>",
# "\n","T-Wave Voltage: <tbd>",
# "\n","Atrium fibrilation: <tbd>"


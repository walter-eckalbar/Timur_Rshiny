
library(shinythemes)
library(shiny)
library(Biobase)

library(ggrepel)
library(reshape)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(cluster)
library(dendextend)

library(plyr)
library(dplyr)
library(ggplot2)
library(genefilter)
library("devtools")
library(Rtsne)


library(gplots)
library("calibrate")

library(R.utils)

options(shiny.maxRequestSize=1000*1024^2)

shinyServer(function(input, output, session) {


  condsTable = read.table("conditions.txt",sep='\t',header=T,row.names = 1)
  data = read.table("data.txt",sep = '\t',header = T,row.names = 1)
  annot = read.table("~/data/human_ens_GRCh38_annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
  
  observe({
    boxNameType <- input$nameType
    #print(boxNameType)
    if (boxNameType == "Ensembl_ID"){
      boxNames <- annotData$Ensembl_ID
    }
    if (boxNameType == "Gene") {
      boxNames <- annotData$Gene
    }
    
    condsNames = colnames(conds)
    
    boxNamesPlusCols = c(condsNames,boxNames)
    
    updateSelectInput(session, "geneLookUpY",
                      label = "Enter gene for the Y-axis",
                      choices = boxNames,
                      selected = boxNames[c(1)])
    
    updateSelectInput(session, "geneLookUpX",
                      label = "Enter condition or gene for the X-axis",
                      choices = boxNamesPlusCols,
                      selected = boxNamesPlusCols[c(1)])
    
    
    observe({
      xValue = input$geneLookUpX
      yValue = input$geneLookUpY
      nameType = input$nameType
      
      if (nameType == "Ensembl_ID"){                            # checks if both are genes 
        if (xValue %in% boxNames) {                              
          xValueGene = annotData[(annotData$Ensembl_ID == xValue),]$Gene
          yValueGene = annotData[(annotData$Ensembl_ID == yValue),]$Gene
          xValueEnsembl = xValue
          yValueEnsembl = yValue
        } else {
          yValueGene = annotData[(annotData$Ensembl_ID == xValue),]$Gene
          yValueEnsembl = yValue
        }
      }
      if (nameType == "Gene"){                            # checks if both are genes 
        if (xValue %in% boxNames) {                              
          xValueEnsembl = annotData[(annotData$Gene == xValue),]$Ensembl_ID
          yValueEnsembl = annotData[(annotData$Gene == yValue),]$Ensembl_ID
          xValueGene = xValue
          yValueGene = yValue
          ens_ids = c(xValueEnsembl,yValueEnsembl)
          
        } else {
          yValueEnsembl = annotData[(annotData$Gene == xValue),]$Ensembl_ID
          yValueGene = yValue
          ens_ids = yValueEnsembl
        }
      }

      selectData = data[(rownames(data) %in% ens_ids),]
      selectData = merge(annot,selectData,by.x='Ensembl_ID',by.y='row.names')
      maxAnnotCol = ncol(annot)
      maxCol = ncol(selectData)
      rownames(selectData) = selectData$Gene
      
      selectData = selectData[,(maxAnnotCol+1):maxCol]
      
      condsData = merge(conds, t(selectData),by= 'row.names')
      
    })
    
    box <- reactive({
      xValue = input$geneLookUpX
      yValue = input$geneLookUpY
      nameType = input$nameType
      
      if (nameType == "Ensembl_ID"){                            # checks if both are genes 
        if (xValue %in% boxNames) {                              
          xValueGene = annotData[(annotData$Ensembl_ID == xValue),]$Gene
          yValueGene = annotData[(annotData$Ensembl_ID == yValue),]$Gene
          xValueEnsembl = xValue
          yValueEnsembl = yValue
        } else {
          yValueGene = annotData[(annotData$Ensembl_ID == yValue),]$Gene
          yValueEnsembl = yValue
        }
      }
      if (nameType == "Gene"){                            # checks if both are genes 
        if (xValue %in% boxNames) {                              
          xValueEnsembl = annotData[(annotData$Gene == xValue),]$Ensembl_ID
          yValueEnsembl = annotData[(annotData$Gene == yValue),]$Ensembl_ID
          xValueGene = xValue
          yValueGene = yValue
          ens_ids = c(xValueEnsembl,yValueEnsembl)
          
        } else {
          yValueEnsembl = annotData[(annotData$Gene == yValue),]$Ensembl_ID
          yValueGene = yValue
          ens_ids = yValueEnsembl
        }
      }
      
      selectData = data[(rownames(data) %in% ens_ids),]
      selectData = merge(annot,selectData,by.x='Ensembl_ID',by.y='row.names')
      maxAnnotCol = ncol(annot)
      maxCol = ncol(selectData)
      rownames(selectData) = selectData$Gene
      
      selectData = selectData[,(maxAnnotCol+1):maxCol]
      print(selectData)
      condsData = merge(conds, t(selectData),by= 'row.names')
      rownames(condsData) = condsData[,c(1)]
      condsData = condsData[,-c(1)]
      
      if (xValue %in% boxNames){
        p <- ggplot(condsData, aes(x=eval(parse(text = xValueGene)),y=eval(parse(text = yValueGene)))) + 
          geom_point() + 
          labs(x=paste("Log Stabilized Expression of ", xValueGene, sep=""), 
                y=paste("Log Stabilized Expression of ", yValueGene, sep=""))
      } else{
        if (is.numeric(condsData[,xValue])){
          p <- ggplot(condsData, aes(x=eval(parse(text = xValue)),y=eval(parse(text = yValueGene)))) + 
            geom_point() + 
            labs(x=xValue, y=paste("Log Stabilized Expression of ", yValueGene, sep=""))
        } else {
          p <- ggplot(condsData, aes(x=factor(eval(parse(text = xValue))),y=eval(parse(text = yValueGene)))) + 
            geom_boxplot(outlier.colour = NA) + 
            geom_point(position = position_jitter(width=.2), size=1) +
            labs(x=xValue, y=paste("Log Stabilized Expression of ", yValueGene, sep=""))
        }
        p
      }
    })
    
    output$box <- renderPlot({
      print(box())
    })
    
    
    output$downloadBoxPlot <- downloadHandler(
      filename = "ShinyBox.png",
      content = function(file) {
        
        png(file, width = 700, height = 400)
        print(box())
        dev.off()
        
      }
    ) 
    
    ## render box plot using volcano function ##
    
    output$plotBox.ui <- renderUI({
      plotOutput("box", width = 700, height = 400)
    })
    
    
  })
})
  

#ensembl_ids = rownames(annot[(annot$Gene %in% genes),])
#selectData = data[(rownames(data) %in% ensembl_ids),]
#selectData = merge(annot,selectData,by.x='Ensembl_ID',by.y='row.names')
#
#
#maxAnnotCol = ncol(annot)
#maxCol = ncol(selectData)
#rownames(selectData) = selectData$Gene
#
#selectData = selectData[,(maxAnnotCol+1):maxCol]
#
#condsData = merge(conds, t(selectData),by= 'row.names')



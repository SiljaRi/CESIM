server <- function(input, output, session) {
    # Initialization #####
    # Hide all outputs
    hide_outputs <- function(){
        output$fishPlot <- renderPlot({plot.new()})
        hide("clonetable")
        hide("mutationtable")
        hide("updatetext1")
        hide("updatetext2")
        hide("updatebtn1")
        hide("updatebtn2")
        hide("uiaddclone")
        hide("uiaddsample")
        hide("uideletesamplevalue")
        hide("uideleteclonevalue")
        hide("uideleteclone")
        hide("uideletesample")
        hide("addcnvbtn")
        hide("deletecnvvalue")
        hide("deletecnvbtn")
        hide("cnvtable")
        hide("puritytable")
        output$button <- renderUI({h5("Choose Basic Settings first and simulate")})
    }
    #Hide outputs at beginning
    hide_outputs()
    
    #Empty mutation table that is diplayed before simulation
    output$mutationtable <- renderDataTable((data.frame(ID = 0, Chr = 0, Start = 0, End = 0, Ref = 0, Var = 0, VAF = 0)))
    show("mutationtable")
    
    #Load Pictures
    output$wwulogo <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/wwulogo.png",
            fileType = 'image/png',
            height = 150,
            alt = "Logo of WWU"))
    }, deleteFile = FALSE)
    
    output$imilogo <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/IMI-logo.png",
            contentType = 'image/png',
            height = 120,
            alt = "Logo of IMI"))
    }, deleteFile = FALSE)
    
    output$independentpng <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/independent.png",
            contentType = 'image/png',
            height = 120,
            width = 200,
            alt = "Logo of IMI"))
    }, deleteFile = FALSE)
    
    output$linearpng <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/linear.png",
            contentType = 'image/png',
            height = 120,
            width = 200,
            alt = "Logo of IMI"))
    }, deleteFile = FALSE)
    
    output$parallelpng <- renderImage({
        
        # Return a list containing the filename
        return(list(
            src = "images/parallel.png",
            contentType = 'image/png',
            height = 120,
            width = 200,
            alt = "Logo of IMI"))
    }, deleteFile = FALSE)
    
    
    # Update mutationcount if clones is modified so that Mutations >= Clones
    observe({
        req(input$clonecount)
        clon <- as.numeric(input$clonecount)
        val <- as.numeric(input$mutationcount)
        if(val < clon){
            val = ceiling(clon*2)
        }
        updateSelectInput(session, "mutationcount", label = "Number of Mutations (SNV+CNV)", choices = clon:200, selected = val)
    })
    
    # Update CNV slider so that max CNV == max Mutations
    observe({
        req(input$mutationcount)
        mut <- as.numeric(input$mutationcount)
        val <- as.numeric(input$cnvcount)
        if(val > mut){
            val = 0
        }
        updateSliderInput(session, "cnvcount", label = "Number of CNVs", min = 0, max = mut, value = val)
    })
    
    # Unselect Checkboxes if CNV = 0 
    observe({
        req(input$cnvcount)
        cnv <- input$cnvcount
        if(cnv == 0){
            updateCheckboxGroupInput(session, "cnvcheckbox", "Explicit", choices = c("Duplication","Deletion","LOH"), selected = NULL)
        }
    })
    
    inputvalues <- reactiveValues(data = NULL)
    
    #Control variables which influence outputs (bad for fishplot)
    rsamples <- eventReactive(input$simulatebtn, {input$samplevalue})
    rclones <- eventReactive(input$simulatebtn, {input$clonecount})
    
    
    # Simulation ####
    observeEvent(input$simulatebtn, {
        if(input$evolutionvalue == "parallel-dependent" && input$clonecount == 2){
            hide_outputs()
            showNotification("Please select at least 3 clones for parallel-dependant evolution.", duration = 5, closeButton = TRUE, type = "error")
        } else {
            req(input$mutationcount)
            req(input$clonecount)
            id <- showNotification("Simulation is running. This can take up to 2 minutes.", duration = NULL, closeButton = FALSE, type = "message")
            # Try to create a phylogeny with all chosen basic settings
            phylogeny <<- create_phylogeny(input$mutationcount, input$clonecount, input$samplevalue, input$evolutionvalue, input$detectionvalue, 0, input$distancevalue)
            if(is.null(phylogeny)){# no phylogeny found
                hide_outputs()
                removeNotification(id)
                showNotification("No resolution found. Try fewer clones or more samples.", duration = 5, closeButton = TRUE, type = "error")
                clones <<- NULL
                mutations <<-NULL
                cnvs <<- NULL
                purities <<-NULL
                showcnvs <<- NULL
            }else{# phylogeny successfully created
                showcnvs <<- NULL
                clones <<- phylogeny[-1,-1]# delete simulation artifacts (Number (column), Ancestor (row))
                #Create CNVs
                if (input$cnvcount > 0 ){
                    cnvs <<- create_cnvs(input$cnvcheckbox, input$cnvcount, clones)
                    showcnvs <<- create_showcnvs(cnvs)
                } else{# no CNVs chosen
                    cnvs <<- NULL
                    showcnvs <<- NULL
                }
                #Initial Puritytable
                purities <<- data.frame(matrix(data = as.numeric(input$purityvalue), nrow = 1, ncol = as.numeric(rsamples())))
                colnames(purities) <<- paste0("t", 1:rsamples())
                # Create SNVs
                mutations <<- create_snvs(clones, input$coveragevalue, purities, 
                                          input$distancevalue, cnvs)[[1]]
                removeNotification(id)
                # Update all tables and plots
                update_outputs()
            }
        }
    })
    
    
    # Modularized Updatefunction (Called by Update- & Delete Buttons) ####
    update  <- function(pclones, pcnvs = NULL){
        hide("updatetext1") # Delete old messages.
        hide("updatetext2")
        u_clones = pclones # needed parameter
        if(any(is.na(u_clones))){
            output$updatetext1 <- renderText("Please fill all cells.")
            output$updatetext2 <- renderText("Please fill all cells (Advanced Settings).")
            show("updatetext1")
            show("updatetext2")
            return()
        }
        if(!is.null(input$puritytable)){
            u_purities <- round(hot_to_r(input$puritytable,0))
            if(any(is.na(u_purities))){
                output$updatetext1 <- renderText("Please fill all cells (Purity).")
                output$updatetext2 <- renderText("Please fill all cells (Purity).")
                show("updatetext1")
                show("updatetext2")
                return()
            }
        }else{
            u_purities <- NULL
        }
        warning <-NULL
        if(!is.null(input$cnvtable)){
            if(is.null(pcnvs)){
                u_cnvs <- hot_to_r(input$cnvtable)
            } else {
                u_cnvs = pcnvs #If a CNV was deleted, the cnv table is given as parameter for checks before deletion.
                if(nrow(u_cnvs) < 1){#Last CNV was deleted
                    u_cnvs <- NULL
                    u_showcnvs <- NULL
                }
            }
            if(nrow(u_cnvs) > 0){
                if(any(is.na(u_cnvs[,-8]))){
                    output$updatetext1 <- renderText("Please fill all cells (CNVs).")
                    output$updatetext2 <- renderText("Please fill all cells (CNVs).")
                    show("updatetext1")
                    show("updatetext2")
                    return()
                }
                # Check CNVs for validity
                cnvcheck <- check_cnvs(u_cnvs, u_clones[,"Mutations"], (sum(u_clones[,"Mutations"])-nrow(u_cnvs)), u_clones[,-(1:2)])
                if(!cnvcheck[[1]]){# Something went wrong. Errormessage in cnvcheck[[2]]
                    output$updatetext1 <- renderText(cnvcheck[[2]])
                    output$updatetext2 <- renderText(cnvcheck[[2]])
                    show("updatetext1")
                    show("updatetext2")
                    return() # no update of intern variables happened
                } else {
                    # CNVs are checked successfully
                    warning <- cnvcheck[[2]]
                    u_cnvs <- cnvcheck[[3]]
                    # Create CNV input for mutation table
                    u_showcnvs <- create_showcnvs(u_cnvs)
                }
            } else{
                u_cnvs <- NULL
                u_showcnvs <- NULL
            }
        }else{
            u_cnvs <- NULL
            u_showcnvs <- NULL
        }
        # CNV and purity are checked. Now Phylogeny table.
        clonecheck <- check_phylogeny(u_clones, input$detectionvalue, input$mutationcount, input$distancevalue, u_cnvs)
        if(clonecheck[[1]]){# CLONES VALID
            clones <<- as.data.frame(clonecheck[[3]])
            if(!is.null(u_purities)){
                purities <<- u_purities# update intern variables
            }
            if(!is.null(u_cnvs)){
                cnvs <<- u_cnvs# update intern variables
            }
            if(!is.null(u_showcnvs)){
                showcnvs <<- u_showcnvs# update intern variables
            }
            id <- showNotification("Updating (creating new mutation data may take some seconds).", duration = NULL, closeButton = FALSE, type = "message")
            # Create new SNV data with updated clones, puritites, cnvs
            mutationresult <- create_snvs(clones, input$coveragevalue, purities, 
                                          input$distancevalue, cnvs)
            removeNotification(id)
            mutations <<- mutationresult[[1]] # update intern variables
            warning <- c(warning, mutationresult[[2]])
            update_outputs()
        }
        output$updatetext1 <- renderText(paste(warning,clonecheck[[2]])) # Show Warning or success message
        output$updatetext2 <- renderText(paste(warning,clonecheck[[2]]))
        show("updatetext1")
        show("updatetext2")
        print("UPDATE COMPLETE!") # intern print
    }
    
    # Update-Button ####
    observeEvent({
        input$updatebtn1 
        input$updatebtn2
    }, {
        if(!is.null(input$updatebtn1) || !is.null(input$updatebtn2)){
            if(!is.null(input$clonetable)){
                u_clones = round(hot_to_r(input$clonetable),2)
                update(u_clones) # modularized update function
            } else { 
                print("No Clonetable found") # intern print
            }
        }
    }, ignoreNULL = FALSE)
    
    
    
    
    # Add Clone ####
    observeEvent(input$addclonebtn,{
        phylodata <- round(hot_to_r(input$clonetable),2)
        new_row <-  data.frame(matrix(c(0,0,1, rep(0, times = ncol(phylodata)-3)), nrow = 1))
        colnames(new_row) <- colnames(phylodata)
        rownames(new_row) <- nrow(phylodata)+1
        DF <- rbind(phylodata, new_row)
        update_outputs(DF)
        output$updatetext1 <- renderText("Empty Clone added. Please fill in values.")
        show("updatetext")
    })
    
    
    # Add Sample ####
    observeEvent(input$addsamplebtn,{
        phylodata <- round(hot_to_r(input$clonetable),2)
        new_column <- data.frame(matrix(rep(0, times = nrow(phylodata)), ncol = 1))
        colnames(new_column) <- paste0("t",ncol(phylodata)-1)
        rownames(new_column) <- rownames(phylodata)
        clones <<- cbind(phylodata, new_column) 
        purities <<- cbind(purities, as.numeric(input$purityvalue))
        colnames(purities) <<- paste0("t",1:ncol(purities))
        update_outputs()
        output$updatetext1 <- renderText("Empty Sample added. Please fill in values.")
        show("updatetext1")
    })
    
    
    # Add CNV ####
    observeEvent(input$addcnvbtn,{
        cnvdata <- hot_to_r(input$cnvtable)
        clonedata <- hot_to_r(input$clonetable)
        samplecount <- ncol(clonedata)-2
        samples <- matrix(rep(0, times = samplecount), nrow = 1)
        names <- paste0("t", 1:samplecount)
        colnames(samples) <- names
        new_row <-  cbind(data.frame(ID = paste("CNV", nrow(cnvdata)+1, sep="_"), Chr = 1, Start = 1, End = 1, Type = NA,  Overlap = FALSE, SNVs = 0, Function = NA, in_Clone = 0), samples)
        rownames(new_row) <- nrow(cnvdata)+1
        DF <- rbind(cnvdata, new_row)
        update_outputs(pcnvs = DF)
        output$updatetext2 <- renderText("Empty CNV added. Please fill in values.")
        show("updatetext2")
    })
    
    
    # Delete Sample ####
    observeEvent(input$deletesamplebtn,{
        samplenumber <- input$deletesamplevalue
        maxsamples <- (ncol(hot_to_r(input$clonetable))-2)
        if(maxsamples < 2){
            output$updatetext1 <- renderText(paste("Last sample can not be deleted."))
            return()
        }
        if(samplenumber < 1 || samplenumber > maxsamples){
            range <- paste0(1,"-",maxsamples)
            output$updatetext1 <- renderText(paste("Invalid Sample Number. Range is", range))
        }else{
            DF <- round(hot_to_r(input$clonetable),2)
            DF <- DF[,-(samplenumber+2)]
            # Check if phylogeny would still be valid after deleting clone
            clonecheck <-  check_phylogeny(DF, input$detectionvalue, input$mutationcount, input$distancevalue, u_cnvs)
            if(clonecheck[[1]]){# Clone can be deleted
                clones <<- as.data.frame(clonecheck[[3]]) # update intern variables
                colnames(clones) <<- c("Parent", "Mutations", paste0("t",1:(ncol(clones)-2)))
                purities <<- data.frame(purities[,-(samplenumber)]) # update intern variables
                if(length(purities)>1){
                    colnames(purities) <<- paste0("t",1:ncol(purities))
                }else{
                    pur <- cbind(purities, 0)
                    colnames(pur) <- c("t1","t2")
                    purities <<- pur[,1]
                }
                cnvs <<- cnvs[,-(samplenumber+9)] # update intern variables
                colnames(cnvs) <<- c("ID", "Chr", "Start", "End", "Type", "Overlap", "SNVs",
                                     "Function", "in_Clone", paste0("t",1:(ncol(clones)-2)))
                showcnvs <<- create_showcnvs(cnvs) # update intern variables
                mutations <<- create_snvs(clones, input$coveragevalue, purities, 
                                          input$distancevalue, cnvs)[[1]] # update intern variables
                update_outputs() # update all tables and plots
                output$updatetext1 <- renderText("Sample deleted successfully.")
            } else { # Clone can not be deleted. Reason in clonecheck[[2]]
                output$updatetext1 <- renderText(clonecheck[[2]])
            }
            show("updatetext1")
        }
    })
    
    
    # Delete Clone ####
    observeEvent(input$deleteclonebtn,{
        clonenumber <- input$deleteclonevalue
        phylodata <- round(hot_to_r(input$clonetable),2)
        maxclones <- nrow(phylodata)
        if(maxclones < 2){
            output$updatetext1 <- renderText(paste("Last Clone can not be deleted."))
            return()
        }
        if(clonenumber < 1 || clonenumber > maxclones){
            range <- paste0(1,"-",maxclones)
            output$updatetext1 <- renderText(paste("Invalid clone number. Range is", range))
        }else{
            DF <- delete_clone(phylodata, clonenumber) # Helper Function in CheckTable.R
            update(DF) # updated with DF corrected by helper function
        }
    })
    
    # Delete CNV ####
    observeEvent(input$deletecnvbtn,{
        cnvnumber <- input$deletecnvvalue
        maxcnv <- nrow(hot_to_r(input$cnvtable))
        if(maxcnv == 0){
            output$updatetext2 <- renderText(paste("No CNV to be deleted."))
        } else if(cnvnumber < 1 || cnvnumber > maxcnv){
            range <- paste0(1,"-",maxcnv)
            output$updatetext2 <- renderText(paste("Invalid CNV Number. Range is", range))
        }else{
            DF <- hot_to_r(input$cnvtable)
            u_cnvs <- DF[-cnvnumber,]
            if(nrow(u_cnvs)>0){#New IDs for ongoing numeration
                rownames(u_cnvs) <- 1:nrow(u_cnvs)
                u_cnvs$ID <- paste("CNV", 1:nrow(u_cnvs), sep="_")
            } 
            DF <- hot_to_r(input$clonetable)
            update(DF, u_cnvs)
        }
    })
    
    # Round values ####
    observeEvent(input$clonetable, { # 2 digits clonefreqs
        clones <<- round(hot_to_r(input$clonetable),2)
        draw_colortable(clones)
    })
    
    observeEvent(input$puritytable, { # 0 digits purity 
        purities <<- round(hot_to_r(input$puritytable),0)
        output$puritytable <- renderRHandsontable({
            rhandsontable(data.frame(purities),rowHeaderWidth = 75, rowHeaders = NULL,
                          colHeaders = colnames(purities)) %>%
                hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
                hot_cols(format = "0", halign = "htCenter") %>%
                hot_validate_numeric(cols = 1:ncol(purities), min = 1, max = 100) 
        })
    })
    
    # Load session ####
    observeEvent(input$upload, { 
        infile <- input$upload
        if(is.null(infile)){
            return(NULL)
        } else {
            hide_outputs()
            numfiles <- nrow(infile)
        }
        if(numfiles != 6){
            output$updatetext1 <- renderText("Please choose exact 6 RDS files (inputs, clones, mutations, cnvs, showcnvs, purities).")
            output$updatetext2 <- renderText("Please choose exact 6 RDS files (inputs, clones, mutations, cnvs, showcnvs, purities).")
            show("updatetext1")
            show("updatetext2")
            return()
        }
        # All 6 files loaded
        inputs <- NULL
        for(i in 1:numfiles){
            switch (infile$name[i],
                    "clones.RDS" =  clones <<- readRDS(infile$datapath[i]),
                    "mutations.RDS" = mutations <<- readRDS(infile$datapath[i]),
                    "cnvs.RDS" =  cnvs <<- readRDS(infile$datapath[i]),
                    "showcnvs.RDS" = showcnvs <<- readRDS(infile$datapath[i]),
                    "purities.RDS" =  purities <<- readRDS(infile$datapath[i]),
                    "inputs.RDS" = inputs <- i,
                    # 6 files chosen, but not every file is correct
                    {  output$updatetext1 <- renderText(paste("Unknown file:", infile$name[i],". Please choose the following 6 RDS files: inputs, clones, mutations, cnvs, showcnvs, purities."))
                    output$updatetext2 <- renderText(paste("Unknown file:", infile$name[i],". Please choose the following 6 RDS files: inputs, clones, mutations, cnvs, showcnvs, purities."))
                    show("updatetext1")
                    show("updatetext2")
                    return()
                    }
            )
        }
        savedInputs <- readRDS(infile$datapath[inputs])
        inputvalues <- unlist(savedInputs)
        inputIDs <- names(inputvalues)
        # restore all setting values
        for (i in 1:length(inputvalues)) { 
            session$sendInputMessage(inputIDs[i],  list(value=inputvalues[[i]]) )
        }
        update_outputs() #update all table and plot outputs
    })
    
    
    # Update Outputs (Visualization)####
    # All intern variables have to be uptaded before. This function only visualizes all outputs.
    update_outputs <- function(pclones = clones, ppurities = purities, pcnvs = cnvs, pmutations = mutations, pshowcnvs = showcnvs){ 
        output$updatetext1 <- renderText("")
        output$updatetext2 <- renderText("")
        
        # Clonetable
        draw_colortable(pclones)
        show("clonetable")
        
        
        # Mutationtable
        if(is.null(pshowcnvs)){ # no CNVs selected
            allmutations <<- rbind(pmutations, pshowcnvs) 
        }else{ # CNVs present
            allmutations <<- rbind(pmutations,setNames(pshowcnvs, names(pmutations)))
        }
        rownames(allmutations) <- 1:nrow(allmutations)
        output$mutationtable <-renderDataTable(allmutations, options = list(
            pageLength = 10,
            autoWidth = TRUE,
            columnDefs = list(list(width = '10px', targets = c(1:5)))
        ))
        show("mutationtable")
        
        # Puritytable
        output$uipuritytitle <- renderUI({
            h4("Purity Table")
        })
        if(length(ppurities)>1){
            if(ncol(ppurities)>ncol(pclones)-2){#Puritytable exploited, shorten
                ppurities <<- ppurities[,1:(ncol(pclones)-2)]
                purities <<- ppurities[,1:(ncol(pclones)-2)]
            } 
            cols <- 1:ncol(purities)
        }else{
            cols <- 1
        }
        t1 <- ppurities
        output$puritytable <- renderRHandsontable({
            rhandsontable(data.frame(t1), rowHeaderWidth = 75, rowHeaders = NULL) %>%
                hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
                hot_cols(format = "0", halign = "htCenter") %>%
                hot_validate_numeric(cols = cols, min = 1, max = 100) %>%
                htmlwidgets::onRender("
                   function(el, x) {
                    var hot = this.hot
                    $('a[data-value=\"Data\"').click(function() {
                      setTimeout(function() {hot.render();}, 0);
                     })
          }")
        })
        
        show("puritytable")
        outputOptions(output, "puritytable", suspendWhenHidden = FALSE) # Needed to create table without changing into Data tab
        
        # CNVtable
        if(is.null(pcnvs)){ # no CNVs present, provide an empty table
            pcnvs <- data.frame(ID = NA, Chr = NA, Start = NA, End = NA, Type = NA,  Overlap = FALSE, "SNVs" = 0, Function = NA, in_Clone = 0)
            pcnvs <- pcnvs[-1,]
        }
        output$cnvtable <- renderRHandsontable({
            rhandsontable(data.frame(pcnvs)) %>%
                hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
                hot_cols(format = "0", halign = "htCenter", allowInvalid = FALSE, readOnly = TRUE) %>%
                hot_col(col = "Type", type = "autocomplete", source = c("Duplication", "Deletion", "LOH"), strict = TRUE,readOnly = FALSE)%>%
                hot_col(col = "Function", type = "autocomplete", source = c("CNV first","SNV first (passive)","SNV first (active)", "Parallel"), strict = TRUE, readOnly = FALSE)%>%
                hot_col(col = c("Chr","End","Start","in_Clone","Overlap","SNVs"), readOnly = FALSE) %>%
                hot_validate_numeric(cols = 2, min = 1, max = 23) %>%
                hot_validate_numeric(cols = 3, min = 1) %>%
                hot_validate_numeric(cols = 4, min = 1) %>%
                hot_col(col = "ID", readOnly = TRUE) %>%
                htmlwidgets::onRender("
                   function(el, x) {
                    var hot = this.hot
                    $('a[data-value=\"Data\"').click(function() {
                      setTimeout(function() {hot.render();}, 0);
                     })
          }")
        })
        
        show("addcnvbtn")
        show("deletecnvvalue")
        show("deletecnvbtn")
        show("uicnvtitle")
        show("cnvtable")
        
        outputOptions(output, "cnvtable", suspendWhenHidden = FALSE) # Needed to create table without changing into Data tab
        
        #Fishplot
        output$fishPlot <- renderPlot({
            draw_fishplot(pclones) # helper function
        })
        show("fishPlot")
        
        # Add all buttons to the screen (Not seen before first simulation)
        output$uiaddclone <- renderUI({
            actionButton("addclonebtn","Add Clone")
        })
        output$uiaddsample <- renderUI({
            actionButton("addsamplebtn","Add Sample")
        })
        output$uideletesamplevalue <- renderUI({
            numericInput("deletesamplevalue", "Choose sample to delete", value = 1, min = 1, width = 500)
        })
        output$uideletesample <- renderUI({
            actionButton("deletesamplebtn", "Delete sample")
        })
        output$uideleteclonevalue <- renderUI({
            numericInput("deleteclonevalue", "Choose clone to delete", value = 1, min = 1, width = 500)
        })
        output$uideleteclone <- renderUI({
            actionButton("deleteclonebtn", "Delete clone")
        })
        output$uicnvtitle <- renderUI({
            h4("Table of Copy Number Variants")
        })
        output$uiexport2 <- renderUI({
            downloadButton("exportbtn1","Export Data")
        })
        
        show("updatebtn1")
        show("updatebtn2")
        show("uiaddclone")
        show("uiaddsample")
        show("uideletesamplevalue")
        show("uideleteclonevalue")
        show("uideleteclone")
        show("uideletesample")
    }
    
    # Draw Colortable ####
    draw_colortable <-function(table){
        colors <- grDevices::rainbow(nrow(table),s = 0.7)
        hot <- rhandsontable(table, ColorCode = colors, rowHeaderWidth = 75) %>% 
            hot_cols(renderer = color_renderer, format = "00.00", max = 100, min = 0, halign = "htCenter", allowInvalid = FALSE) %>%
            hot_table(highlightCol = TRUE, highlightRow = TRUE,  contextMenu = FALSE ) %>%
            hot_col("Parent", format = "0", halign = "htCenter") 
        output$clonetable <- renderRHandsontable({hot})
        show("clonetable")
    }
    
    # Draw Fishplot ####
    draw_fishplot <- function(pclones){
        samplecount = ncol(pclones)-2
        clonecount = nrow(pclones)
        fractable = as.matrix(pclones[,-(1:2)]) #P arent+Mutations cutted
        if(samplecount == 1){# Trick to make 1 timepoint possible (Doubling it)
            timepoints <- 1:2
            fractable <- cbind(fractable,fractable)
            vlab = c("t_1","")
            vlines = 1
        }else{
            timepoints <- 1:samplecount
            vlab = paste0("t_", 1:samplecount)
            vlines = timepoints
        }
        parents <- as.vector(pclones[,1])
        fish = createFishObject(fractable, parents, timepoints=timepoints, fix.missing.clones = TRUE)
        fish = layoutClones(fish)
        fish = setCol(fish, grDevices::rainbow(clonecount,s = 0.7))
        fishPlot(fish,shape="spline", vlines = vlines, col.vline = "black", cex.title=0.7, 
                 pad.left=0.7, bg.type = "gradient", bg.col = "white", vlab=vlab)
    }
    
    # Data Export ####
    output$exportbtn1  <- downloadHandler(
        filename = function() {paste("Simulation-", Sys.time(), ".zip", sep="")}, # For unique downloadfiles
        content = function(file) {
            
            # EXCELFILE with all data
            wb <- createWorkbook()
            addWorksheet(wb, "Settings")
            addWorksheet(wb, "Plot")
            addWorksheet(wb, "Phylogeny")
            addWorksheet(wb, "Mutations")
            addWorksheet(wb, "SNVs")
            addWorksheet(wb, "CNVs")
            addWorksheet(wb, "Purity")
            
            # Settings 
            # All selected values are compared to theis actual value. Both are printed.
            if(is.null(input$cnvcheckbox)){ # handling checkbox input
                selvalues = c(input$samplevalue, input$evolutionvalue, input$clonecount, input$mutationcount, input$cnvcount, NA,
                              input$coveragevalue, input$purityvalue, input$detectionvalue, input$distancevalue)
                
            } else {
                cnvvalues <- input$cnvcheckbox
                cnvvalues <- paste(cnvvalues, collapse = ",")
                selvalues = c(input$samplevalue, input$evolutionvalue, input$clonecount, input$mutationcount, input$cnvcount, cnvvalues,
                              input$coveragevalue, input$purityvalue, input$detectionvalue, input$distancevalue)
            }
            if(is.null(cnvs)){
                actvalues = c(ncol(clones)-2, "-", nrow(clones), nrow(allmutations), 0, "-", "-", 
                              "(see purity tab)", "-", "-")
            } else {
                actvalues = c(ncol(clones)-2, "-", nrow(clones), nrow(allmutations), nrow(cnvs), "-", "-", 
                              "(see purity tab)", "-", "-")
            }
            settings <- data.frame(
                Setting = c("Samples", "Evolution", "Clones", "Mutations", "CNVs", "Explicit", "Coverage", "Purity", "Detection", "Distance"),
                Selected = selvalues,
                Actual = actvalues
            )
            writeData(wb, "Settings", settings)
            writeData(wb, "Settings", paste("Warning: These Settings could have been overwritten by your manual changes on tables (see differences in Actual)."), startCol = 1, startRow = 15)
            writeData(wb, "Settings", paste("- means that parameter is not trackable."), startCol = 1, startRow = 16)
            setColWidths(wb, "Settings", c(1,2), widths = 10)
            setColWidths(wb, "Settings", 3, widths = 14)
            addStyle(wb, "Settings",style = createStyle(halign = "center"), rows = 1:11, cols = 1:3, gridExpand = TRUE)
            
            # Fishplot
            fishplot <- draw_fishplot(clones)
            insertPlot(wb, "Plot", width = 10, height = 5)
            dev.off()
            
            # Phylogeny
            writeData(wb, "Phylogeny", clones)
            writeData(wb, "Mutations", allmutations)
            writeData(wb, "SNVs", mutations)
            writeData(wb, "CNVs", cnvs)
            writeData(wb, "Purity", purities)
            
            # Files are saved in zip and in local directory, therefore I created directory Output to seperate from code files
            saveWorkbook(wb, file = "Output/Data.xlsx", overwrite = TRUE)
            
            # Single tables for easy further processing
            write.table(clones, file = "Output/clones.tsv", row.names = FALSE, sep ="\t")
            write.table(cnvs, file = "Output/cnv.tsv", row.names = FALSE, sep ="\t")
            write.table(mutations, file = "Output/snv.tsv", row.names = FALSE, sep ="\t")
            write.table(purities, file = "Output/purity.tsv", row.names= FALSE, sep="\t")
            
            if (!file.exists("Output/RDS")){
                dir.create(file.path("Output/RDS"))
            }
            
            # Save current session in RDS directory for restoring possibility
            saveRDS( reactiveValuesToList(input) , file = 'Output/RDS/inputs.RDS')
            saveRDS( clones, "Output/RDS/clones.RDS")
            saveRDS( mutations, "Output/RDS/mutations.RDS")
            saveRDS( cnvs, "Output/RDS/cnvs.RDS")
            saveRDS( purities, "Output/RDS/purities.RDS")
            saveRDS( showcnvs, "Output/RDS/showcnvs.RDS")
            
            fs = c("Data.xlsx", "Output/Data.xlsx", "Output/clones.tsv", "Output/cnv.tsv", "Output/purity.tsv", "Output/snv.tsv","Output/RDS/inputs.RDS", 
                   "Output/RDS/clones.RDS", "Output/RDS/mutations.RDS",
                   "Output/RDS/cnvs.RDS", "Output/RDS/purities.RDS", "Output/RDS/showcnvs.RDS")
            
            zip(zipfile = file, files = fs)
        },
        contentType = "appclication/zip"
    )
    
}
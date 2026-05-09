##before server stuff
##Runs dependencies required to run the application

library(shiny)
library(DT)

#here initalize data using readData.R

source("readData.R") #read the data in using this function

#read in the entire database for displaying and test data
allCompsDisplay <- read.csv("All_comps.csv", encoding = "latin1")
testData <- read.csv("test_data.csv", encoding = "latin1")

##server

server <- function(input, output) {
  
  #set max upload size to 100 Mb to prevent crashing
  options(shiny.maxRequestSize = 100*1024^2)
  
  if (file.exists("count.txt")) {
    saved_count <- as.integer(readLines("count.txt", warn = FALSE))
    if (is.na(saved_count)) saved_count <- 0
  } else {
    saved_count <- 0
  }
  
  if (!file.exists("visit_count.txt")) {
    writeLines("0", "visit_count.txt")
  }

    isolate({
      count <- as.integer(readLines("visit_count.txt", warn = FALSE))
      if (is.na(count)) count <- 0
      count <- count + 1
      writeLines(as.character(count), "visit_count.txt")
    })
    
    # Optional: display to user
    output$visit_count <- renderText({
      paste("This app has been accessed", readLines("visit_count.txt"), "times.")
    })
  
  #initalizes reactive (public) variables
  rvals <- reactiveValues(
    upload = NULL,
    expt.masses = NULL,
    exptRTs = NULL,
    results = NULL,
    results_appended = NULL,
    results_appended_1 = NULL,
    upload2.custom = NULL,
    custom.expt.masses = NULL,
    custom.results = NULL,
    custom.results.appended = NULL,
    custom.results.appended_1 = NULL,
    userdata = NULL,
    name_adduct_mass = NULL,
    v3 = NULL,
    name_and_PRT = NULL,
    custom_RTS = NULL,
    custom.col.names = NULL,
    download.custom.appended = NULL,
    download.noncustom = NULL,
    sample.names = NULL,
    input_vector = NULL,
    timesUsed = NULL,
    searchCount = saved_count
  )
  
    #this is for the input of the file
    observeEvent(input$file1, {
      req(input$file1)
      
      # Read uploaded CSV
      uploaded.data <- read.csv(input$file1$datapath, header = TRUE)
      
      # Save entire uploaded dataset
      rvals$upload <- uploaded.data
      
      # Pull column names
      header.names <- colnames(uploaded.data)
      
      # Extract sample names (columns 3 to N)
      if (length(header.names) > 2) {
        rvals$sample.names <- header.names[3:length(header.names)]
      } else {
        rvals$sample.names <- "N/A"
      }
      
      # Extract m/z values from column 1
      rvals$expt.masses <- uploaded.data[[1]]
      rvals$exptRTs <- uploaded.data[[2]]
      
      # Save selected dataset input (e.g., "1", "2", "3")
      rvals$input_vector <- input$dataset
    })
    
    #select what databases to use based on the user inputs
    db_map <- list("1" = pubchem, "2" = WOS)
    combined_db <- reactive({
      req(input$WOS_pubChem)
      do.call(rbind, db_map[input$WOS_pubChem])
    })
    
    #select what RTs to use based on user inputs
    rt_map <- list("1" = PRTs_pubchem, "2" = PRTs_WOS)
    combined_rt <- reactive({
      req(input$WOS_pubChem)
      do.call(rbind, rt_map[input$WOS_pubChem])
    })
    
    adduct_group_map <- list(
      "1" = c("Monoisotopic"),
      "2" = c("M.plus.H"),
      "3" = c("M.minus.H"),
      "4" = c("2M.plus.H", "2M.plus.NH4", "M.plus.K", "M.plus.Na", "M.plus.NH4")
    )
    
    filtered_adducts <- reactive({
      req(input$dataset, combined_db())
      selected_adducts <- unlist(adduct_group_map[input$dataset])
      combined <- combined_db()
      combined[combined$Adduct %in% selected_adducts, ]
    })
    
    matched_data <- reactive({
      req(
        rvals$expt.masses,
        rvals$exptRTs,
        rvals$upload,
        filtered_adducts(),
        input$tol,
        input$tolMode,
        combined_rt()
      )
      
      expt_data <- rvals$upload
      expt_mz <- rvals$expt.masses
      expt_rt <- rvals$exptRTs
      db <- filtered_adducts()
      rt_table <- combined_rt()
      result <- db[0, ]
      
      for (i in seq_along(expt_mz)) {
        mz <- expt_mz[i]
        rt <- as.numeric(expt_rt[i])  # ensure numeric
        
        delta <- if (input$tolMode == 2) {
          (input$tol / 1e6) * mz  # ppm
        } else {
          input$tol  # dalton
        }
        
        match_rows <- db[db$mz >= (mz - delta) & db$mz <= (mz + delta), ]
        if (nrow(match_rows) == 0) next
        
        match_rows$Expt_mz <- mz
        match_rows$Expt_RT <- rt
        match_rows$ppm_error <- 1e6 * abs(match_rows$mz - mz) / mz
        
        # Safe merge on compound name
        col1_match <- colnames(match_rows)[1]
        col1_rt <- colnames(rt_table)[1]
        colnames(match_rows)[1] <- "merge_key"
        colnames(rt_table)[1] <- "merge_key"
        
        merged_rows <- merge(match_rows, rt_table, by = "merge_key", all.x = TRUE)
        
        colnames(merged_rows)[colnames(merged_rows) == "merge_key"] <- col1_match
        colnames(rt_table)[1] <- col1_rt
        
        if (nrow(merged_rows) == 0) next
        
        # Sample intensities
        if (ncol(expt_data) > 2) {
          user_row <- expt_data[i, 3:ncol(expt_data), drop = FALSE]
          rep_count <- nrow(merged_rows)
          
          if (rep_count > 0) {
            user_intensity_df <- user_row[rep(1, rep_count), , drop = FALSE]
            if (nrow(user_intensity_df) == rep_count) {
              merged_rows <- cbind(merged_rows, user_intensity_df)
            }
          }
        }
        
        result <- rbind(result, merged_rows)
        
      }
        
      return(result)
        
    })
    
    observeEvent(matched_data(), {
      rvals$searchCount <- isolate(rvals$searchCount + 1)
      
      log_entry <- paste0(Sys.time(), " | Search #", rvals$searchCount, "\n")
      cat(log_entry, file = "search_log.txt", append = TRUE)
      
      #save the updated count
      writeLines(as.character(rvals$searchCount), "count.txt")
    })
    
    output$search_count <- renderText({
      paste("Tool has been used", rvals$searchCount, "times.")
    })
    
    #for making table
    output$matches <- renderDataTable({
      data <- matched_data()
      selected <- data[, 1:7, drop = FALSE]
      
      # If column "ppm_error" is one of those first 7 columns
      if ("ppm_error" %in% colnames(selected)) {
        selected$ppm_error <- round(selected$ppm_error, 5)
      }
      
      # Reorder
      new_order <- c(1, 2, 3, 4, 6, 5, 7)
      
      selected <- selected[, new_order]
      
      #need to compute RT error
      percent_rt_match = 100*(pmin((selected[,6])/selected[,7],(selected[,7]/selected[,6]))) 
      RT_match = abs(round(percent_rt_match, 2))
      
      selected <- cbind(selected, RT_match)
      
      colnames(selected) <- c(
        "Compound Name", "Adduct", "Database m/z", "Experimental m/z",
        "ppm Error", "Expt RT", "Predicted RT", "% RT Match"
      )
      
      selected
    })

    #for download
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0("matched_data_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv")
      },
      content = function(file) {
        data <- matched_data()
        if (nrow(data) == 0) {
          writeLines("No matched results", con = file)
        } else {
          #reorder first 7 columns, move column 6 to position 4
          if (ncol(data) >= 7) {
            reordered_cols <- c(1, 2, 3, 4, 6, 5, 7)
            first7 <- data[, reordered_cols]
            remaining <- data[, -(1:7), drop = FALSE]
            percent_rt_match = 100*(pmin((data[,5])/data[,7],(data[,7]/data[,5])))
            RT_match = abs(round(percent_rt_match, 2))
            data <- cbind(first7, RT_match, remaining)
          }
          write.csv(data, file, row.names = FALSE)
        }
      }
    )
    
    #To view all compounds as a table
    allCompsDisplay <- read.csv("All_comps.csv", encoding = "latin1")
    
    # Render the full table
    output$allCompsDisplay <- renderDataTable({
      allCompsDisplay
    })
    
    #for downloading all aminoacids
    output$downloadAAs <- downloadHandler(
      filename = function() {
        "all_aminoacids.csv"
      },
      content = function(file) {
        write.csv(allCompsDisplay, file, row.names = FALSE)
      }
    )
    
    #for downloading test data
    output$downloadTestData <- downloadHandler(
      filename = function() {
        "test_data.csv"
      },
      content = function(file) {
        write.csv(testData, file, row.names = FALSE)
      }
    )
    
    observeEvent(input$dataset_input,{
      req(input$dataset_input)
      rvals$userdata <- read.csv(input$dataset_input$datapath, header = TRUE,)
      userdata <- rvals$userdata
      colnames(userdata) <- c("Compound.Name", "Adduct", "Mass", "RTP")
      
      rvals$name_adduct_mass <- userdata[,1:3]
      
      custom.unique.compounds <- length(unique(userdata[,1]))
      rvals$name_and_PRT <- userdata[1:custom.unique.compounds,]
      rvals$name_and_PRT <- rvals$name_and_PRT[,c(1,4)]
      
    })
    
    observeEvent(input$file2,{
      req(input$file2)
      #reads in the .csv file and saves to global envrionment, set header = TRUE
      rvals$upload2.custom <- read.csv(input$file2$datapath, header = TRUE,)
      
      uploaded.data <- rvals$upload2.custom
      custom.col.names <- colnames(uploaded.data[,3:ncol(uploaded.data)])
      rvals$custom.expt.masses <- uploaded.data[,1]
      
      name_adduct_mass <- rvals$name_adduct_mass
      v3 <- name_adduct_mass[,3]
      
      if (length(input$dataset_input) > 0){
        
        #create an empty vector for markers to go into
        custom.marker.vect <- c()
        
        for (i in 1:length(rvals$custom.expt.masses)){
          
          if (input$tolMode2 == 1){
            
            marker.custom <- which(rvals$custom.expt.masses[i] >= v3-input$tol2 & rvals$custom.expt.masses[i] <= v3+input$tol2)
            
            if (length(marker.custom[i]) > 0){
              
              custom.marker.vect[[i]] <- marker.custom
              
            }
          } else if (input$tolMode2 == 2){
            
            #vector for absolute ppm error to go in (will all be positive)
            ppmAbsolute2 <- c()
            ppmTol2 <- c()
            
            #calculate the error for each entry in the database
            ppmAbsolute2 <- v3*(1+input$tol2/10^6) 
            ppmTol2 <-  ppmAbsolute2 - v3
            
            marker.custom <- which(rvals$custom.expt.masses[i] >= v3-ppmTol2 & rvals$custom.expt.masses[i] <= v3+ppmTol2)
            
            if (length(marker.custom[i]) > 0){
              
              custom.marker.vect[[i]] <- marker.custom
              
            }
            
          }
        }
        #empty vector for position of experimental matches to go
        custom.exp.match <- c()
        for (j in 1:length(custom.marker.vect)){
          
          if (length(unlist(custom.marker.vect[j])) > 0){
            
            custom.no.match <- length(unlist(custom.marker.vect[j]))
            custom.vvv <- (rep(j, length(unlist(custom.marker.vect[j]))))
            custom.exp.match[[j]] <- custom.vvv
            
          }
        }
        
        # gives the matches from the data base as a vector
        custom.location.db <- as.vector(unlist(custom.marker.vect))
        
        # prints out everything from db where there is a match
        custom.db.hits <- name_adduct_mass[custom.location.db,]
        
        # gives the postion of experimental hits from uploaded data
        custom.pos.exp <- unlist(custom.exp.match)
        custom.results <- cbind(custom.db.hits, rvals$upload2.custom[custom.pos.exp,])
        colnames(custom.results) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT")
        custom.print_results <- custom.results[,1:5]
        rvals$custom.results <- custom.print_results
        rvals$custom.results.appended <- custom.results
      }
      
    })
    
    #####Custom database search, copy pasted straight from HormonomoicsDB 
    output$contents_shell <- renderDataTable({
      
      validate(
        #eliminates the error message, so that the error message is much friendlier
        need(input$file2 != "", label = "A .csv file with m/z and RT values")
      )
      
      #psuedo SQL left join
      custom.result.rt <- merge(x = rvals$custom.results.appended[,1:5], y = rvals$name_and_PRT,
                                by.x = "Compound Name", by.y = "Compound.Name", all.x = TRUE)
      
      #rename cols
      colnames(custom.result.rt) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT", "Predicted RT")
      
      #computes difference between expt rt and predicted rt
      #custom.delta.rt <- abs(custom.result.rt[,5] - custom.result.rt[,6]) #changed to below
      custom.percent.delta.rt <- 100*(pmin((custom.result.rt[,6])/custom.result.rt[,7],(custom.result.rt[,7]/custom.result.rt[,6]))) #added 20250823
      
      #converts to %
      #custom.percent.delta.rt <- 100-(custom.delta.rt/custom.result.rt[,5])*100 #changed to above
      
      #binds the percent diff to the results
      custom.results.rt.delta.final <- cbind(custom.result.rt, custom.percent.delta.rt)
      
      #changes column names
      colnames(custom.results.rt.delta.final) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                                   "Experimental m/z", "RT", "Predicted RT",
                                                   "Percent Match RT", rvals$custom.col.names)
      
      #creats msrt for each hit
      custom.results.rt.delta.final.mzrt <- paste0(custom.results.rt.delta.final[,4], "_", custom.results.rt.delta.final[,5])
      
      #calculate ppm error and add it to the output data for download.
      ppmErrorVect <- ((custom.results.rt.delta.final[,3] - custom.results.rt.delta.final[,4])/custom.results.rt.delta.final[,3])*10^6
      
      #bind the ppmError vect to the other data
      custom.results.for.download <- cbind(custom.results.rt.delta.final, ppmErrorVect)
      
      #binds mz_rt to rest of results
      custom.results.for.download <- cbind(custom.results.for.download, custom.results.rt.delta.final.mzrt)
      
      #column names
      colnames(custom.results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                                 "Experimental m/z", "RT", "Predicted RT",
                                                 "Percent Match RT","ppmError", "mzrt")
      
      shell.results.sorted <- custom.results.rt.delta.final[order(custom.results.rt.delta.final[,5]),]
      
      #takes intensities and brings into local envrionment
      custom.experimental.intensities <- rvals$custom.results.appended
      
      #creates mz_rt for each hit
      custom.experimental.intensities.mzrt <- paste0(custom.experimental.intensities[,4],
                                                     "_", custom.experimental.intensities[,5])
      
      #binds back together
      custom.experimental.intensities <- cbind(custom.experimental.intensities.mzrt,
                                               custom.experimental.intensities[,6:ncol(custom.experimental.intensities)])
      
      #removes duplicates
      custom.experimental.intensities <- unique(custom.experimental.intensities)
      
      #renaming the mz_rt column
      colnames(custom.experimental.intensities) <- c("mzrt1")
      
      #inner join
      custom.download.results <- merge(x = custom.results.for.download, y = custom.experimental.intensities,
                                       by.x = "mzrt", by.y = "mzrt1", all.x = TRUE)
      
      #DF for downloading
      custom.download.results <- data.frame(custom.download.results)
      
      #send to global envrio for download through GUI
      rvals$download.custom.appended <- custom.download.results
      
      #compute ppm mass error (don't think people will want to sort by ppm?)
      ppmError2 <- ((shell.results.sorted[,3] - shell.results.sorted[,4])/shell.results.sorted[,4])*10^6
      
      shell.results.sorted <- cbind(shell.results.sorted, ppmError2)
      
      colnames(shell.results.sorted) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                          "RT", "Predicted RT", "Percent Match RT", "ppm Mass Difference")
      
      # Increment the counter
      rvals$searchCount <- isolate(rvals$searchCount + 1)
      
      # Log the event
      log_entry <- paste0(Sys.time(), " | Custom Search #", rvals$searchCount, "\n")
      cat(log_entry, file = "search_log.txt", append = TRUE)
      
      # Save updated count to file
      writeLines(as.character(rvals$searchCount), "count.txt")
      
      #output final results into GUI for user to see
      shell.results.sorted[,1:8]
      
    },
    digits = 5 #always displays 4 decimal points
    )
    
    #download custom data
    output$downloadData_shell <- downloadHandler(
      filename = function(){
        #gives a unique name each time
        paste("aminoacidDB_custom_search_results_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file){
        #writes it to a .csv, reactive so it changes all the time
        write.csv(rvals$download.custom.appended, file)
      }
    )
    
}
  
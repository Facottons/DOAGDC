#' Separate patients in groups
#'
#' \code{groups_identification_cox} is a function designed to separate patients in
#' groups, based in clinical data, and coxHR.
#'
#' @param dataType
#' @param Width
#' @param Height
#' @param Res
#' @param Unit
#' @param image_format
#' @param saveData
#' @param env
#' @param tumor
#' @param dataBase
#' @param workDir
#' @param Name
#' @inheritParams download_gdc
#' @inheritParams concatenate_files
#' @inheritParams groups_identification_mclust
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' groups_identification_coxHR("HIF3A", "gene", tumor = "UCS", dataBase = "legacy")
#' }
groups_identification_coxHR <- function(Name,
                                      dataType,
                                      Width = 3000,
                                      Height = 2000,
                                      Res = 400,
                                      Unit = "px",
                                      image_format = "png",
                                      saveData = TRUE, env,
                                      tumor, dataBase,
                                      workDir) {

    # Library load
    # (xlsx tools httr jsonlite data.table XML ggplot2 survival survminer ggthemes grid yarrr reshape extrafont scales cgdsr)
    #  (plyr png)


    dataType <- gsub(" ", "_", dataType)
    Name <- gsub("-", "_", Name)

    PATH <- file.path(workDir, "DOAGDC", toupper(tumor), "Analyses", tolower(Name))

    dir.create(path = file.path(workDir, "DOAGDC", toupper(tumor), "Analyses"),
               showWarnings = FALSE)

    dir.create(path = file.path(workDir, "DOAGDC", toupper(tumor), "Analyses", tolower(Name)),
               showWarnings = FALSE)

    dir.create(file.path(PATH,
                         paste0("/survival_Results_", toupper(Name))),
               showWarnings = FALSE)
    DIR <- paste0(PATH, "/survival_Results_", toupper(Name))
    dir.create(file.path(DIR, "coxHR"), showWarnings = FALSE) #1-cutoff_finder
    dir.create(file.path(DIR, "kaplan_maier"), showWarnings = FALSE) #2-splitting_n_clinical

    # Download clinical
    if (tolower(dataBase) == "legacy") {
        download_gdc(dataType = "clinical_supplement", tumor = tumor, dataBase = "legacy", workDir = workDir)
        folder_name <- "clinical_supplement_data"
    } else if (tolower(dataBase) == "gdc") {
        download_gdc(dataType = "clinical_supplement", tumor = tumor, dataBase = "gdc", workDir = workDir)
        folder_name <- "gdc_clinical_supplement_data"
    }

    # local functions ####
    # geomSeries <- function(base, max) {
    #     base^(0:floor(log(max, base)))
    # }

    # part1 ####

    # PART B: CLINICAL EVALUATION
    # Add available files

    MANIFEST <- data.frame(read.table(file=file.path(workDir, "DOAGDC", toupper(tumor),
                                                     folder_name, "manifest.sdrf"),
                                      stringsAsFactors = FALSE, header=TRUE, sep="\t"))

    AvailableFiles <- MANIFEST$file_name

    # Create table to receive the relevant clinical points (for all tumors)
    clinicalData_selected <- data.frame(matrix(ncol=0, nrow=length(AvailableFiles)))

    # Add empty values to the columns
    clinicalData_selected[, c("bcr_patient_barcode",
                             "gender",
                             "race",
                             "ethnicity",
                             "days_to_birth",
                             "age_at_initial_pathologic_diagnosis",
                             "tissue_source_site",
                             "icd_o_3_site",
                             "icd_o_3_histology",
                             "vital_status_initial",
                             "days_to_last_followup_initial",
                             "days_to_death_initial",
                             "person_neoplasm_cancer_status_initial",
                             "system_version",
                             "pathologic_stage",
                             "pathologic_T",
                             "pathologic_N",
                             "pathologic_M",
                             "new_tumor_event_after_initial_treatment_initial",
                             "days_to_new_tumor_event_after_initial_treatment_initial",
                             "new_neoplasm_event_type_initial",
                             "new_neoplasm_event_occurrence_anatomic_site_initial",
                             "form_completion_YYYMMDD_initial",
                             "followups",
                             "form_completion_YY_followup",
                             "form_completion_MM_followup",
                             "form_completion_DD_followup",
                             "person_neoplasm_cancer_status_followup",
                             "vital_status_followup",
                             "days_to_last_followup_followup",
                             "new_tumor_event_after_initial_treatment_followup",
                             "new_tumor_event_after_initial_treatment_number_followup",
                             "days_to_new_tumor_event_after_initial_treatment_followup",
                             "new_neoplasm_event_type_followup",
                             "new_neoplasm_event_occurrence_anatomic_site_followup",
                             "days_to_death_followup"
    )] <- "unknown"


    # Loop to fill the general info data
    for(fileNow in 1:length(AvailableFiles)){

        #
        # fileNow <- match("nationwidechildrens.org_clinical.TCGA-BH-A18L.xml", AvailableFiles)
        #Parse xml_data
        # data <- XML::xmlParse(paste("./", tumor, "/", "clinical_xml", "/Raw/", AvailableFiles[fileNow],sep=""))
        data <- XML::xmlParse(file.path(workDir, "DOAGDC", toupper(tumor),
                                        folder_name, AvailableFiles[fileNow]))
        xml_data <- XML::xmlToList(data)


        #### data-to-collect
        # Patient barcode
        clinicalData_selected$bcr_patient_barcode[fileNow] <- xml_data$patient$bcr_patient_barcode$text
        # Primary site (sanity check)
        # clinicalData_selected$tumor_tissue_site[fileNow] <- xml_data$patient$tumor_tissue_site$text

        # Gender
        clinicalData_selected$gender[fileNow] <- xml_data$patient$gender$text

        # Race
        clinicalData_selected$race[fileNow] <- xml_data$patient$race[[1]][[1]]

        #Etnicidade
        clinicalData_selected$ethnicity[fileNow] <- xml_data$patient$ethnicity[[1]]

        # Days to birth
        clinicalData_selected$days_to_birth[fileNow]  <- xml_data$patient$days_to_birth[[1]]

        # age_at_initial_pathologic_diagnosis
        # clinicalData_selected$age_at_initial_pathologic_diagnosis[fileNow] <- xml_data$patient$age_at_initial_pathologic_diagnosis[[1]]

        # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes
        # Place where tissuewere collected
        clinicalData_selected$tissue_source_site[fileNow] <- xml_data$patient$tissue_source_site[[1]]

        #
        # exact icd location site
        # https://training.seer.cancer.gov/breast/abstract-code-stage/codes.html
        clinicalData_selected$icd_o_3_site[fileNow] <- xml_data$patient$icd_o_3_site[[1]]

        # Histological name
        # clinicalData_selected$histological_type[fileNow] <- xml_data$patient$histological_type[[1]]

        # Histological definition
        #http://codes.iarc.fr/codegroup/2
        clinicalData_selected$icd_o_3_histology[fileNow] <- xml_data$patient$icd_o_3_histology[[1]]

        # Vital status at first report
        clinicalData_selected$vital_status_initial[fileNow] <- xml_data$patient$vital_status[[1]]

        # days_to_last_followup_initial
        clinicalData_selected$days_to_last_followup_initial[fileNow] <- xml_data$patient$days_to_last_followup[[1]]

        # days_to_death_initial
        clinicalData_selected$days_to_death_initial[fileNow] <- xml_data$patient$days_to_death[[1]]

        # actual neoplasm status
        if(is.null(xml_data$patient$person_neoplasm_cancer_status[[1]])){
        } else {
            clinicalData_selected$person_neoplasm_cancer_status_initial[fileNow] <- xml_data$patient$person_neoplasm_cancer_status[[1]]
        }


        #Stage check
        if(is.list(xml_data$patient$stage_event)){
            # Staging system
            clinicalData_selected$system_version[fileNow] <- xml_data$patient$stage_event$system_version[[1]]

            # Staging system
            clinicalData_selected$pathologic_stage[fileNow] <- xml_data$patient$stage_event$pathologic_stage[[1]]

            # pathologic_T
            clinicalData_selected$pathologic_T[fileNow] <- xml_data$patient$stage_event$tnm_categories$pathologic_categories$pathologic_T[[1]]

            # pathologic_N
            clinicalData_selected$pathologic_N[fileNow] <- xml_data$patient$stage_event$tnm_categories$pathologic_categories$pathologic_N[[1]]

            # pathologic_M
            clinicalData_selected$pathologic_M[fileNow] <- xml_data$patient$stage_event$tnm_categories$pathologic_categories$pathologic_M[[1]]

        } else {
            # dont do anything
        }


        # (NTE)
        if(is.list(xml_data$patient$new_tumor_events$new_tumor_event)){

            # New tumors events (NTE)
            clinicalData_selected$new_tumor_event_after_initial_treatment_initial[fileNow] <- xml_data$patient$new_tumor_events$new_tumor_event_after_initial_treatment[[1]]

            if(is.null(xml_data$patient$new_tumor_events$new_tumor_event$days_to_new_tumor_event_after_initial_treatment[[1]])){
                clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_initial[fileNow] <- "unknown"
            } else {
                #  days to New tumors events (NTE)
                clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_initial[fileNow] <- xml_data$patient$new_tumor_events$new_tumor_event$days_to_new_tumor_event_after_initial_treatment[[1]]

            }

            if(is.null(xml_data$patient$new_tumor_events$new_tumor_event$new_neoplasm_event_type[[1]])){
                clinicalData_selected$new_neoplasm_event_type_initial[fileNow] <- "unknown"
            } else {
                if(length(xml_data$patient$new_tumor_events$new_tumor_event$new_neoplasm_event_type[[1]])==1){
                    # new_tumor_event_type ***
                    clinicalData_selected$new_neoplasm_event_type_initial[fileNow] <- xml_data$patient$new_tumor_events$new_tumor_event$new_neoplasm_event_type[[1]]

                } else {}

            }

            if(is.null(xml_data$patient$new_tumor_events$new_tumor_event$new_neoplasm_event_occurrence_anatomic_site[[1]])){
                clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_initial[fileNow] <- "unknown"
            } else {
                # new_tumor_event_site ***
                clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_initial[fileNow] <- xml_data$patient$new_tumor_events$new_tumor_event$new_neoplasm_event_occurrence_anatomic_site[[1]]
            }


        }


        # FORM COMPLETION
        if(is.list(xml_data$patient$day_of_form_completion)){
            if(length((xml_data$patient$day_of_form_completion)) == 1){
                clinicalData_selected$form_completion_YYYMMDD_initial[fileNow] <- paste(xml_data$patient$year_of_form_completion$text,
                                                                                        xml_data$patient$month_of_form_completion$text,
                                                                                        xml_data$patient$day_of_form_completion$text,
                                                                                        sep="_")
            } else {

            }

        } else {
            if(length((xml_data$patient$month_of_form_completion)) == 1){
                clinicalData_selected$form_completion_YYYMMDD_initial[fileNow] <- paste(xml_data$patient$year_of_form_completion$text,
                                                                                        xml_data$patient$month_of_form_completion$text,
                                                                                        "unknown",
                                                                                        sep="_")
            } else {}

        }

        # IF THERE ARE FOLLOW UPS
        #
        clinicalData_selected$followups[fileNow] <- 0

        if(is.list(xml_data$patient$follow_ups)){

            #Collect number of followups
            clinicalData_selected$followups[fileNow] <- length(xml_data$patient$follow_ups)

            # Collect year
            clinicalData_selected$form_completion_YY_followup[fileNow] <-
                paste(sapply(1:length(xml_data$patient$follow_ups),
                             function(w){
                                 xml_data$patient$follow_ups[[w]]$year_of_form_completion$text
                             }), collapse="_")

            # Collect month
            clinicalData_selected$form_completion_MM_followup[fileNow] <-
                paste(sapply(1:length(xml_data$patient$follow_ups),
                             function(w){
                                 xml_data$patient$follow_ups[[w]]$month_of_form_completion$text
                             }), collapse="_")

            # Collect day
            clinicalData_selected$form_completion_DD_followup[fileNow] <-
                paste(sapply(1:length(xml_data$patient$follow_ups),
                             function(w){
                                 xml_data$patient$follow_ups[[w]]$day_of_form_completion[[1]]
                             }), collapse="_")
            #
            clinicalData_selected$form_completion_DD_followup[fileNow] <- gsub("form_completion_day", "unknown", clinicalData_selected$form_completion_DD_followup[fileNow])


            # New tumors events (NTE)
            clinicalData_selected$person_neoplasm_cancer_status_followup[fileNow] <- "_"
            clinicalData_selected$vital_status_followup[fileNow] <- "_"
            clinicalData_selected$days_to_last_followup_followup[fileNow] <- "_"
            clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow] <- "_"
            clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup[fileNow] <- "_"
            clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[fileNow] <- "_"
            clinicalData_selected$new_neoplasm_event_type_followup[fileNow] <- "_"
            clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup[fileNow] <- "_"
            clinicalData_selected$days_to_death_followup[fileNow] <- "_"

            ##### GO FOR FOLLOWUPS
            for(followupnow in 1:length(xml_data$patient$follow_ups)){

                #
                clinicalData_selected$person_neoplasm_cancer_status_followup[fileNow] <- paste(clinicalData_selected$person_neoplasm_cancer_status_followup[fileNow],
                                                                                               xml_data$patient$follow_ups[[followupnow]]$person_neoplasm_cancer_status[[1]],
                                                                                               sep="_")

                #
                clinicalData_selected$vital_status_followup[fileNow] <-
                    paste(clinicalData_selected$vital_status_followup[fileNow],
                          xml_data$patient$follow_ups[[followupnow]]$vital_status[[1]],
                          sep="_")

                #
                clinicalData_selected$days_to_last_followup_followup[fileNow] <-
                    paste(clinicalData_selected$days_to_last_followup_followup[fileNow],
                          xml_data$patient$follow_ups[[followupnow]]$days_to_last_followup[[1]],
                          sep="_")

                #
                clinicalData_selected$days_to_death_followup[fileNow] <-
                    paste(clinicalData_selected$days_to_death_followup[fileNow],
                          xml_data$patient$follow_ups[[followupnow]]$days_to_death[[1]],
                          sep="_")


                #
                if(xml_data$patient$follow_ups[[followupnow]]$.attrs[1] == "4.0"){
                    if(xml_data$patient$follow_ups[[followupnow]]$new_tumor_events[[1]][[1]] == "YES"){

                        # Add yes
                        clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow],
                                  xml_data$patient$follow_ups[[followupnow]]$new_tumor_events[[1]][[1]],
                                  sep="_")

                        # Add NTE number
                        clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup[fileNow],
                                  (length(xml_data$patient$follow_ups[[followupnow]]$new_tumor_events)-1),
                                  sep="_")

                        #
                        for(newtumoreventsnow in 2:length(xml_data$patient$follow_ups[[followupnow]]$new_tumor_events)){
                            #
                            clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[fileNow] <-
                                paste(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[fileNow],
                                      xml_data$patient$follow_ups[[followupnow]]$new_tumor_events[[newtumoreventsnow]]$days_to_new_tumor_event_after_initial_treatment[[1]],
                                      sep="_")

                            clinicalData_selected$new_neoplasm_event_type_followup[fileNow] <-
                                paste(clinicalData_selected$new_neoplasm_event_type_followup[fileNow],
                                      xml_data$patient$follow_ups[[followupnow]]$new_tumor_events[[newtumoreventsnow]]$new_neoplasm_event_type[[1]],
                                      sep="_")

                            clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup[fileNow] <-
                                paste(clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup[fileNow],
                                      xml_data$patient$follow_ups[[followupnow]]$new_tumor_events[[newtumoreventsnow]]$new_neoplasm_event_occurrence_anatomic_site[[1]],
                                      sep="_")
                        }
                    } else {
                        # Add NO as NTE
                        clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow],
                                  xml_data$patient$follow_ups[[followupnow]]$new_tumor_events[[1]][[1]],
                                  sep="_")
                    }
                }

                # Version check 2.1
                if(xml_data$patient$follow_ups[[followupnow]]$.attrs[1] == "2.1"){
                    if(xml_data$patient$follow_ups[[followupnow]]$new_tumor_event_after_initial_treatment[[1]] == "YES"){

                        # Add yes
                        clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow],
                                  xml_data$patient$follow_ups[[followupnow]]$new_tumor_event_after_initial_treatment[[1]],
                                  sep="_")

                        # Add NTE number
                        clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup[fileNow],
                                  1,
                                  sep="_")

                        # day to new tumor event
                        clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[fileNow] <-
                            paste(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[fileNow],
                                  xml_data$patient$follow_ups[[followupnow]]$days_to_new_tumor_event_after_initial_treatment[[1]],
                                  sep="_")

                        # NTE type
                        clinicalData_selected$new_neoplasm_event_type_followup[fileNow] <-
                            paste(clinicalData_selected$new_neoplasm_event_type_followup[fileNow],
                                  xml_data$patient$follow_ups[[followupnow]]$new_neoplasm_event_type[[1]],
                                  sep="_")

                        # NTE site
                        clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup[fileNow] <-
                            paste(clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup[fileNow],
                                  xml_data$patient$follow_ups[[followupnow]]$new_neoplasm_event_occurrence_anatomic_site[[1]],
                                  sep="_")

                    } else {
                        # Add NO as NTE
                        clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow],
                                  xml_data$patient$follow_ups[[followupnow]]$new_tumor_event_after_initial_treatment[[1]],
                                  sep="_")

                    }
                }


                # Version check 1.5
                if(xml_data$patient$follow_ups[[followupnow]]$.attrs[1] == "1.5"){

                    # Add yes
                    if(is.numeric(xml_data$patient$follow_ups[[followupnow]]$days_to_new_tumor_event_after_initial_treatment[[1]])){
                        clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow],
                                  "YES",
                                  sep="_")
                    } else {
                        clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow] <-
                            paste(clinicalData_selected$new_tumor_event_after_initial_treatment_followup[fileNow],
                                  "NO",
                                  sep="_")
                    }

                    # Add NTE number
                    clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup[fileNow] <-
                        paste(clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup[fileNow],
                              "unknown",
                              sep="_")

                    # day to new tumor event
                    clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[fileNow] <-
                        paste(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[fileNow],
                              xml_data$patient$follow_ups[[followupnow]]$days_to_new_tumor_event_after_initial_treatment[[1]],
                              sep="_")

                    # NTE type
                    clinicalData_selected$new_neoplasm_event_type_followup[fileNow] <-
                        paste(clinicalData_selected$new_neoplasm_event_type_followup[fileNow],
                              "unknown",
                              sep="_")

                    # NTE site
                    clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup[fileNow] <-
                        paste(clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup[fileNow],
                              "unknown",
                              sep="_")


                }
            }
        }
    }

    # ADD IMPORTANT COLUMNS
    # for overall survival, relapse free and distant metastasis free survival

    # Table fixing
    clinicalData_selected$days_to_death_initial <- gsub("day", "unknown", clinicalData_selected$days_to_death_initial)
    clinicalData_selected$days_to_death_initial <- gsub("false", "unknown", clinicalData_selected$days_to_death_initial)

    clinicalData_selected$person_neoplasm_cancer_status_followup <- gsub("__", "", clinicalData_selected$person_neoplasm_cancer_status_followup)
    clinicalData_selected$person_neoplasm_cancer_status_followup <- gsub("tumor_status", "unknown", clinicalData_selected$person_neoplasm_cancer_status_followup)

    clinicalData_selected$vital_status_followup <- gsub("__", "", clinicalData_selected$vital_status_followup)

    clinicalData_selected$days_to_last_followup_followup <- gsub("__", "", clinicalData_selected$days_to_last_followup_followup)
    clinicalData_selected$days_to_last_followup_followup <- gsub("day", "unknown", clinicalData_selected$days_to_last_followup_followup)
    clinicalData_selected$days_to_last_followup_followup <- gsub("false", "unknown", clinicalData_selected$days_to_last_followup_followup)

    clinicalData_selected$new_tumor_event_after_initial_treatment_followup <- gsub("__", "", clinicalData_selected$new_tumor_event_after_initial_treatment_followup)
    clinicalData_selected$new_tumor_event_after_initial_treatment_followup <- gsub("new_tumor_event_dx_indicator", "unknown", clinicalData_selected$new_tumor_event_after_initial_treatment_followup)

    clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup <- gsub("__", "", clinicalData_selected$new_tumor_event_after_initial_treatment_number_followup)

    clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup <- gsub("__", "", clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup)
    clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup <- gsub("false", "unknown", clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup)

    clinicalData_selected$new_neoplasm_event_type_followup <- gsub("__", "", clinicalData_selected$new_neoplasm_event_type_followup)
    clinicalData_selected$new_neoplasm_event_type_followup <- gsub("new_tumor_event_type", "unknown", clinicalData_selected$new_neoplasm_event_type_followup)

    #
    clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup <- gsub("__", "", clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup)
    clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup <- gsub("new_tumor_event_site", "unknown", clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup)
    clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup <- gsub("Other, specify", "Other", clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup)

    #
    clinicalData_selected$days_to_death_followup <- gsub("__", "", clinicalData_selected$days_to_death_followup)
    clinicalData_selected$days_to_death_followup <- gsub("false", "unknown", clinicalData_selected$days_to_death_followup)
    clinicalData_selected$days_to_death_followup <- gsub("day", "unknown", clinicalData_selected$days_to_death_followup)

    # Save (problema1)
    # rownames(clinicalData_selected) <- clinicalData_selected$bcr_patient_barcode

    #
    # Some random list error
    clinicalData_selected$new_neoplasm_event_type_initial <- unlist(clinicalData_selected$new_neoplasm_event_type_initial)

    # OVERALL SURVIVAL
    clinicalData_selected$final_vital_status <- "NA"
    clinicalData_selected$final_vital_status_times <- "NA"

    for(lineNow in 1:dim(clinicalData_selected)[1]){
        # If there are no followup
        if(clinicalData_selected[lineNow,"followups"] == 0){
            # put vital status
            clinicalData_selected$final_vital_status[lineNow] <- clinicalData_selected[lineNow,"vital_status_initial"]
            #
            if(clinicalData_selected$final_vital_status[lineNow] == "Dead"){
                # put dead time
                clinicalData_selected$final_vital_status_times[lineNow] <- clinicalData_selected[lineNow,"days_to_death_initial"]
            } else {
                # put alive until now
                clinicalData_selected$final_vital_status_times[lineNow] <- clinicalData_selected[lineNow,"days_to_last_followup_initial"]
            }
            # if there are follow ups

        } else {

            # check if dead
            if(("Dead" %in% unlist(strsplit(clinicalData_selected[lineNow,"vital_status_followup"], "_")))){
                # Put dead
                clinicalData_selected$final_vital_status[lineNow] <- "Dead"
                # collect times
                clinicalData_selected$final_vital_status_times[lineNow] <- unlist(strsplit(clinicalData_selected[lineNow,"days_to_death_followup"], "_"))[match("Dead", unlist(strsplit(clinicalData_selected[lineNow,"vital_status_followup"], "_")))]
            } else {
                # alive
                clinicalData_selected$final_vital_status[lineNow] <- "Alive"

                # collect times
                # check if there are any number
                if(grepl("\\d", clinicalData_selected[lineNow,"days_to_last_followup_followup"])){
                    timenow <- unlist(strsplit(clinicalData_selected[lineNow,"days_to_last_followup_followup"], "_"))[max(which(unlist(strsplit(clinicalData_selected[lineNow,"vital_status_followup"], "_")) == "Alive"))]
                } else {
                    timenow <- "unknown"
                }

                # Check to access if timenow is unknown
                if(timenow %in% c("unknown", NA)){
                    clinicalData_selected$final_vital_status_times[lineNow] <- "unknown"
                } else {
                    # check first if original is unknown
                    if(clinicalData_selected[lineNow,"days_to_last_followup_initial"] == "unknown"){
                        clinicalData_selected$final_vital_status_times[lineNow] <- "unknown"
                    } else {

                        if(is.na(as.numeric(clinicalData_selected[lineNow,"days_to_last_followup_initial"]))){
                            clinicalData_selected$final_vital_status_times[lineNow] <- timenow
                        } else {
                            if(as.numeric(timenow) > as.numeric(clinicalData_selected[lineNow,"days_to_last_followup_initial"])){
                                clinicalData_selected$final_vital_status_times[lineNow] <- timenow
                            } else {
                                clinicalData_selected$final_vital_status_times[lineNow] <- clinicalData_selected[lineNow,"days_to_last_followup_initial"]
                            }
                        }


                    }}

            }
        }
    }


    # FINAL RFS (RELAPSE FREE SURVIVAL)
    clinicalData_selected$final_rfs_status <- "NA"
    clinicalData_selected$final_rfs_status_times <- "NA"

    for(lineNow in 1:dim(clinicalData_selected)[1]){

        # If there are no followup
        if(clinicalData_selected[lineNow,"followups"] == 0){

            # Check if dead on begin
            if(clinicalData_selected[lineNow,"vital_status_initial"] == "Dead"){
                # check if it had tumor
                if(clinicalData_selected[lineNow,"person_neoplasm_cancer_status_initial"] == "WITH TUMOR"){
                    #ok
                    clinicalData_selected$final_rfs_status[lineNow] <- "relapse"
                    clinicalData_selected$final_rfs_status_times[lineNow] <- clinicalData_selected[lineNow,"days_to_new_tumor_event_after_initial_treatment_initial"]
                } else {
                    # Discard
                    clinicalData_selected$final_rfs_status[lineNow] <- "dead_other_cause_before_relapse"
                    clinicalData_selected$final_rfs_status_times[lineNow] <- "dead_other_cause_before_relapse"
                }


            } else {

                # No follow up but not dead at begin
                if(clinicalData_selected$person_neoplasm_cancer_status_initial[lineNow] == "TUMOR FREE"){
                    clinicalData_selected$final_rfs_status[lineNow] <- "censored"
                    clinicalData_selected$final_rfs_status_times[lineNow] <- clinicalData_selected$days_to_last_followup_initial[lineNow]
                } else {
                    # if is unknown
                    if(clinicalData_selected$person_neoplasm_cancer_status_initial[lineNow] == "unknown"){
                        clinicalData_selected$final_rfs_status[lineNow] <- "unknown"
                        clinicalData_selected$final_rfs_status_times[lineNow] <- "unknown"
                    } else {
                        # If initial NTE is present
                        if(clinicalData_selected$new_tumor_event_after_initial_treatment_initial[lineNow] == "YES"){
                            clinicalData_selected$final_rfs_status[lineNow] <- "relapse"
                            clinicalData_selected$final_rfs_status_times[lineNow] <- clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_initial[lineNow]
                        } else {
                            clinicalData_selected$final_rfs_status[lineNow] <- "censored"
                            clinicalData_selected$final_rfs_status_times[lineNow] <- clinicalData_selected$days_to_last_followup_initial[lineNow]
                        }
                    }

                }
            }
            # If there are any followup
        } else {

            # NTE vector creation
            nte_status_vector <-
                c(clinicalData_selected$new_tumor_event_after_initial_treatment_initial[lineNow],
                  unlist(strsplit(clinicalData_selected$new_tumor_event_after_initial_treatment_followup[lineNow], "_")))

            tumor_status_vector <-
                c(clinicalData_selected$person_neoplasm_cancer_status_initial[lineNow],
                  unlist(strsplit(clinicalData_selected$person_neoplasm_cancer_status_followup[lineNow], "_")))

            life_status_vector <-
                c(clinicalData_selected$vital_status_initial[lineNow],
                  unlist(strsplit(clinicalData_selected$vital_status_followup[lineNow], "_")))

            nte_timing_vector <-
                c(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_initial[lineNow],
                  unlist(strsplit(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[lineNow], "_")))

            followup_timing_vector <-
                c(clinicalData_selected$days_to_last_followup_initial[lineNow],
                  unlist(strsplit(clinicalData_selected$days_to_last_followup_followup[lineNow], "_")))

            dead_timing_vector <-
                c(clinicalData_selected$days_to_death_initial[lineNow],
                  unlist(strsplit(clinicalData_selected$days_to_death_followup[lineNow], "_")))

            # If there are any NTE
            if("YES" %in% nte_status_vector){
                # add relapse
                clinicalData_selected$final_rfs_status[lineNow] <- "relapse"
                # Check for times
                # Collect the smaller one
                clinicalData_selected$final_rfs_status_times[lineNow] <-
                    min(as.numeric(nte_timing_vector),
                        as.numeric(dead_timing_vector),
                        as.numeric(followup_timing_vector[match("YES", nte_status_vector)]), na.rm=TRUE)
            } else {

                # If there are anything in nte_timing_vector
                if(TRUE %in% (as.numeric(nte_timing_vector)>0)){

                    clinicalData_selected$final_rfs_status[lineNow] <- "relapse"
                    # Add this date
                    clinicalData_selected$final_rfs_status_times[lineNow] <- min(as.numeric(nte_timing_vector), na.rm=TRUE)

                    #If there are any WITH TUMOR
                    #if("WITH TUMOR" %in% tumor_status_vector[2:length(tumor_status_vector)]){
                    #  clinicalData_selected$final_rfs_status[lineNow] <- "relapse"
                    # Check for times
                    # Collect the smaller one
                    #  clinicalData_selected$final_rfs_status_times[lineNow] <-
                    #    min(as.numeric(nte_timing_vector),
                    #        as.numeric(dead_timing_vector),
                    #        as.numeric(followup_timing_vector[match("WITH TUMOR", tumor_status_vector)]), na.rm=TRUE)
                } else {
                    if("TUMOR FREE" %in% tumor_status_vector){
                        clinicalData_selected$final_rfs_status[lineNow] <- "censored"
                        # collect the maximum position with TUMOR FREE
                        clinicalData_selected$final_rfs_status_times[lineNow] <- followup_timing_vector[max(which(match(tumor_status_vector, "TUMOR FREE") == 1))]
                    } else {
                        if(TRUE %in% (followup_timing_vector > -100)){
                            clinicalData_selected$final_rfs_status[lineNow] <- "censored"
                            # collect the maximum position with values
                            clinicalData_selected$final_rfs_status_times[lineNow] <- followup_timing_vector[max(which(followup_timing_vector > 0))]
                        } else {
                            if("unknown" %in% tumor_status_vector){
                                clinicalData_selected$final_rfs_status[lineNow] <- "unknown"
                                # Check for times
                                clinicalData_selected$final_rfs_status_times[lineNow] <- "unknown"
                            } else {
                                clinicalData_selected$final_rfs_status[lineNow] <- "censored"
                                clinicalData_selected$final_rfs_status_times[lineNow] <-
                                    max(as.numeric(dead_timing_vector),
                                        as.numeric(followup_timing_vector), na.rm=TRUE)
                            }
                        }
                    }
                }
            }
        }
    }

    # DMFS DISTANT METASTASIS FREE SURVIVAL
    clinicalData_selected$final_dmfs_status <- "unknown"
    clinicalData_selected$final_dmfs_status_times <- "unknown"

    for(lineNow in 1:dim(clinicalData_selected)[1]){

        # If there are no relapse, no metastasis, just copy
        if(clinicalData_selected$final_rfs_status[lineNow] == "censored"){
            # copy
            clinicalData_selected$final_dmfs_status[lineNow] <- clinicalData_selected$final_rfs_status[lineNow]
            #
            clinicalData_selected$final_dmfs_status_times[lineNow] <- clinicalData_selected$final_rfs_status_times[lineNow]

            # if there was relapse, check if it was metastatic
        } else {
            if(clinicalData_selected$final_rfs_status[lineNow] == "relapse"){

                # If there are metastasis in first form, just copy
                if(clinicalData_selected$new_neoplasm_event_type_initial[lineNow] == "Distant Metastasis"){

                    # copy
                    clinicalData_selected$final_dmfs_status[lineNow] <- "metastasis"
                    #
                    clinicalData_selected$final_dmfs_status_times[lineNow] <- clinicalData_selected$final_rfs_status_times[lineNow]

                    # Check for further metastasis
                } else {}

                # Create some vectors
                nte_type_vector <-
                    c(clinicalData_selected$new_neoplasm_event_type_initial[lineNow],
                      unlist(strsplit(clinicalData_selected$new_neoplasm_event_type_followup[lineNow], "_")))

                nte_timing_vector <-
                    c(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_initial[lineNow],
                      unlist(strsplit(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[lineNow], "_")))


                clinicalData_selected$new_neoplasm_event_type_followup
                clinicalData_selected$new_neoplasm_event_occurrence_anatomic_site_followup

                # There are any metastasis on vector
                if("Distant Metastasis" %in% nte_type_vector){
                    # copy
                    clinicalData_selected$final_dmfs_status[lineNow] <- "metastasis"
                    #
                    clinicalData_selected$final_dmfs_status_times[lineNow] <- nte_timing_vector[match("Distant Metastasis", nte_type_vector)]

                    # if there are NO metastasis in vector, put the highest available date
                    # even with other kind of relapses
                } else {

                    dead_timing_vector <-
                        c(clinicalData_selected$days_to_death_initial[lineNow],
                          unlist(strsplit(clinicalData_selected$days_to_death_followup[lineNow], "_")))

                    followup_timing_vector <-
                        c(clinicalData_selected$days_to_last_followup_initial[lineNow],
                          unlist(strsplit(clinicalData_selected$days_to_last_followup_followup[lineNow], "_")))

                    nte_timing_vector <-
                        c(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_initial[lineNow],
                          unlist(strsplit(clinicalData_selected$days_to_new_tumor_event_after_initial_treatment_followup[lineNow], "_")))

                    # copy
                    clinicalData_selected$final_dmfs_status[lineNow] <- "censored"
                    #
                    clinicalData_selected$final_dmfs_status_times[lineNow] <- max(as.numeric(c(dead_timing_vector, followup_timing_vector, nte_timing_vector)), na.rm=TRUE)
                }
            }
        }
    }

    # Save table
    write.csv(clinicalData_selected, file.path(DIR, paste0(tumor, "_clinical_xml", ".csv")))

    # part2 ####

    # Table File position

    if(missing(env)){
        # envir_link <- paste(toupper(tumor), toupper(dataBase), gsub(" ", "_", tolower(dataType)), "tumor_data", sep = "_")
        stop(message("Please, insert the Environment name as 'env' argument!"))
    } else {
        envir_link <- deparse(substitute(env))
    }

    string_vars <- list(envir_link = get(envir_link))

    assign("clinical_xml", clinicalData_selected, envir = get(envir_link))

    #get expression values
    if (tolower(dataType) == "methylation") {
        # all <- string_vars[["envir_link"]]$gene_tumor_normalized
        tmp <- paste0("methylation_tumor_filtered_selected_", Name)
    } else {
        tmp <- paste0(tolower(dataType), "_tumor_normalized_selected_", toupper(Name))
    }

    framesList <- eval(parse(text=paste0('string_vars[["envir_link"]]$', paste0(tmp))))

    framesList <- as.data.frame(framesList)

    # duplicate fix
    # tmp <- unname(sapply(rownames(PRAD_LEGACY_isoform_tumor_data$isoform_tumor_normalized_selected_Name), function(w){
    #     paste(unlist(strsplit(w, "-"))[1:3], collapse="-")
    # }))
    # tmp2 <- cbind(PRAD_LEGACY_isoform_tumor_data$isoform_tumor_normalized_selected_Name, tmp)
    # tmp2$tmp <- as.character(unlist(tmp2$tmp))
    # tmp3 <- as.data.frame(matrix(nrow = length(unique(tmp2$tmp)), ncol = 2))
    # count <- 0
    # for (patients in unique(tmp2$tmp)){
    #     count <- count + 1
    #     selector <- patients == tmp2$tmp
    #     tmp3[count, ] <- tmp2[selector, ][1, ]
    # }
    # rownames(tmp3) <- tmp3[, 2]
    # colnames(tmp3)[1] <- Name
    # framesList <- tmp3[, -2, drop = FALSE]
    # remove(tmp2, tmp3)

    # clinical
    clinicalList <- clinicalData_selected
    rownames(clinicalList) <- clinicalList$bcr_patient_barcode


    # PATIENT SPLITTING METRICS

    # create empty table for filling
    # best_z_table <- matrix(ncol=1, nrow=length(framesList))
    # colnames(best_z_table) <- Name
    # rownames(best_z_table) <- names(framesList)
    # best_z_table <- as.data.frame(best_z_table)

    best_z_table <- NULL

    # add patient code
    framesList$PatientCode <- rownames(framesList)

    # collect ordering
    matchingOrdering <- match(framesList[,"PatientCode"], rownames(clinicalList))

    # Add clinical values to the frames list
    framesList <- cbind(framesList,
                                   clinicalList[matchingOrdering, c("race", "gender", "pathologic_stage", "pathologic_T", "pathologic_N", "pathologic_M",
                                                                               "final_vital_status", "final_vital_status_times", "final_rfs_status", "final_rfs_status_times",
                                                                               "final_dmfs_status", "final_dmfs_status_times")])
    #
    remove(matchingOrdering)

    # Loop for each gene
    # Collect geneNow name


    # # create z-score just for tumor samples
    # TissueTypeSimple <- unname(sapply(rownames(framesList), function(w){
    #     paste(unlist(strsplit(w, "-"))[4], collapse="-")}))
    # TissueTypeSimple <- unname(sapply(TissueTypeSimple, function(w){
    #     paste(unlist(strsplit(w, ""))[1:2], collapse="")}))
    # framesList[,"TissueTypeSimple"] <- TissueTypeSimple
    # tumorValues_RSEM_UQ <- framesList[TissueTypeSimple == "01", Name]
    #
    # # Calculate z-scores
    # tumorValueszscore <- (tumorValues_RSEM_UQ - mean(tumorValues_RSEM_UQ)) / sd(tumorValues_RSEM_UQ)
    #
    # # Save values to table
    # framesList[, "log2p1"] <- log2(framesList[, Name]+1)
    # framesList[, "zscore"] <- NA
    # framesList[framesList$TissueTypeSimple == "01", "zscore"] <- tumorValueszscore
    #
    # # Collect gene with
    # TableNow <- framesList[TissueTypeSimple == "01", ]


    TissueTypeSimple <- unname(sapply(rownames(framesList), function(w){
        paste(unlist(strsplit(w, "-"))[4], collapse="-")}))
    TissueTypeSimple <- unname(sapply(TissueTypeSimple, function(w){
        paste(unlist(strsplit(w, ""))[1:2], collapse="")}))
    framesList[,"TissueTypeSimple"] <- TissueTypeSimple

    # Calculate z-scores
    tumorValueszscore <- as.numeric(scale(framesList[, Name]))

    # Save values to table
    framesList[, "log2p1"] <- log2(framesList[, Name]+1)
    framesList[, "zscore"] <- NA
    framesList[, "zscore"] <- tumorValueszscore

    # Collect gene with
    TableNow <- framesList

    # Sort by expression level
    TableNow <- TableNow[order(TableNow[, "zscore"]), ]

    ### Remove missing data
    TableNow <- TableNow[(!TableNow$final_rfs_status == "unknown"),]
    TableNow <- TableNow[(!TableNow$final_rfs_status == "dead_other_cause_before_relapse"),]
    TableNow <- TableNow[(!TableNow$final_rfs_status == "vital_status"),]
    TableNow <- TableNow[!is.na(TableNow$final_rfs_status),]
    TableNow <- TableNow[!is.nan(TableNow$final_rfs_status),]
    TableNow <- TableNow[(!TableNow$final_rfs_status_times == "unknown"),]
    TableNow <- TableNow[(!TableNow$final_rfs_status_times == "dead_other_cause_before_relapse"),]
    TableNow <- TableNow[(!TableNow$final_rfs_status_times == "day"),]
    TableNow <- TableNow[!is.na(TableNow$final_rfs_status_times),]
    TableNow <- TableNow[!is.nan(TableNow$final_rfs_status_times),]
    TableNow$final_rfs_status_times <- as.numeric(TableNow$final_rfs_status_times)

    # remove patients without relapse information
    TableNow$final_rfs_status <- as.integer(gsub("relapse",1, gsub("censored", 0, TableNow$final_rfs_status)))


    #
    # Start optimization
    # as Cutoff Finder: A Comprehensive and Straightforward Web Application
    #            Enabling Rapid Biomarker Cutoff Optimization
    #

    # define optimization ranges using at least 4% of the patients
    # or 10
    nmin <- round(0.04*dim(TableNow)[1])
    if(nmin < 10){
        nmin <- 10
    }
    nmax <- dim(TableNow)[1] - nmin
    patient_cuts_sequence <- nmin:nmax

    # create matrix to receive the results
    cutoff_optimization <- matrix(nrow=length(patient_cuts_sequence), ncol=6)
    colnames(cutoff_optimization) <- c("patientNo", "zscorecut", "HR", "HR_min95", "HR_max95", "logrank_p")
    cutoff_optimization <- as.data.frame(cutoff_optimization)
    cutoff_optimization$patientNo <- patient_cuts_sequence

    for(patientcutnow in patient_cuts_sequence){
        if(TableNow[, 1][patientcutnow] == 0){
            # check if it isnt the last zero
            if(TableNow[, 1][patientcutnow+1] == 0){
                # CoxHR
                cutoff_optimization[patientcutnow,"HR"] <- Inf

                # CoxHR min .95
                cutoff_optimization[patientcutnow,"HRmin95"] <- Inf

                # CoxHR max .95
                cutoff_optimization[patientcutnow,"HRmax95"] <- Inf

                # logrank pvalue
                cutoff_optimization[patientcutnow,"logrank_p"] <- 1

                # if is the last zero, do normally
            } else {
                # Create the classification vector for now
                classification_vector <- c(rep("low", patientcutnow), rep("high", dim(TableNow)[1]-patientcutnow))
                classification_vector <- factor(classification_vector, levels=c("low", "high"))

                # create the model summary
                coxmodel <- suppressWarnings(summary(survival::coxph(survival::Surv(final_rfs_status_times, final_rfs_status) ~ classification_vector,
                                                     data = TableNow)))

                # collect line now for table
                lineNow <- match(patientcutnow, cutoff_optimization[,"patientNo"])

                cutoff_optimization[lineNow,"zscorecut"] <- TableNow[patientcutnow, "zscore"]

                # CoxHR
                cutoff_optimization[lineNow,"HR"] <- coxmodel$conf.int[1]

                # CoxHR min .95
                cutoff_optimization[lineNow,"HR_min95"] <- coxmodel$conf.int[3]

                # CoxHR max .95
                cutoff_optimization[lineNow,"HR_max95"] <- coxmodel$conf.int[4]

                # logrank pvalue
                cutoff_optimization[lineNow,"logrank_p"] <- coxmodel$logtest[3]
            }
        } else {

            # Create the classification vector for now
            classification_vector <- c(rep("low", patientcutnow), rep("high", dim(TableNow)[1]-patientcutnow))
            classification_vector <- factor(classification_vector, levels=c("low", "high"))

            # create the model summary
            coxmodel <- suppressWarnings(summary(survival::coxph(survival::Surv(final_rfs_status_times, final_rfs_status) ~ classification_vector,
                                      data = TableNow)))

            # collect line now for table
            lineNow <- match(patientcutnow, cutoff_optimization[,"patientNo"])

            cutoff_optimization[lineNow,"zscorecut"] <- TableNow[patientcutnow, "zscore"]

            # CoxHR
            cutoff_optimization[lineNow,"HR"] <- coxmodel$conf.int[1]
            # CoxHR min .95
            cutoff_optimization[lineNow,"HR_min95"] <- coxmodel$conf.int[3]
            # CoxHR max .95
            cutoff_optimization[lineNow,"HR_max95"] <- coxmodel$conf.int[4]
            # logrank pvalue
            cutoff_optimization[lineNow,"logrank_p"] <- coxmodel$logtest[3]
        }
    }

    # remove infinite rows
    # pick infinite rows
    removeRows <- -sort(unique(c(which(is.infinite(cutoff_optimization$HR)),
                                 which(is.infinite(cutoff_optimization$HR_min95)),
                                 which(is.infinite(cutoff_optimization$HR_max95)),
                                 which(is.na(cutoff_optimization$logrank_p)))))

    cutoff_optimization <- cutoff_optimization[!is.na(cutoff_optimization$HR_min95), ]

    # save optim table
    write.csv(cutoff_optimization,
              file=file.path(DIR, "coxHR", paste0("Optimization_", Name, ".csv")))

    if (nrow(cutoff_optimization) == 0) {
        stop(message("\nLow expression values. Please, try another expression data!\n"))
    }

    # Collect best z score
    best_z_table <- suppressWarnings(cutoff_optimization$zscorecut[which(min(cutoff_optimization$logrank_p,
                                                                             na.rm=TRUE) == cutoff_optimization$logrank_p)][1])

    if(nrow(cutoff_optimization[removeRows, ]) > 0){
        # plot
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "coxHR",
                                     paste0("with_confidence_interval_", Name, ".png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "coxHR",
                                     paste0("with_confidence_interval_", Name, ".svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(5,5,3,5))
        plot(x=cutoff_optimization$zscorecut,
             y=cutoff_optimization$HR, type="l",
             ylim=c(min(cutoff_optimization$HR_min95, na.rm = TRUE),
                    max(cutoff_optimization$HR_max95, na.rm = TRUE)),
             xlim=c(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                    max(cutoff_optimization$zscorecut, na.rm = TRUE)),
             lwd=2, lty=1, xlab="Expression z-score", ylab="Cox HR RFS", axes=FALSE,
             col=rgb(8,88,158, 240, max=255))

        axis(side = 1, lwd = 2, cex.axis=1.5, cex=1.5,
             at=round(seq(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                          max(cutoff_optimization$zscorecut, na.rm = TRUE), by = 0.1), 2))
        axis(side = 2, lwd = 2, cex.axis=1.5, cex=1.5,
             col=rgb(8,88,158, 255, max=255), las = 1)

        lines(x=cutoff_optimization$zscorecut,
              y=cutoff_optimization$HR_min95,
              lwd=1, lty=3, col=rgb(19, 153, 24, 200, max=255))
        lines(x=cutoff_optimization$zscorecut,
              y=cutoff_optimization$HR_max95,
              lwd=1, lty=3, col=rgb(19, 153, 24, 200, max=255))

        # abline(v=best_z_table, lwd=2, lty=3, col=rgb(0,0,0, 255, max=255))
        text(best_z_table, max(cutoff_optimization$HR, na.rm = TRUE)+.12, paste0("z-score = ", round(best_z_table, 2)), cex=0.8)
        points(best_z_table, max(cutoff_optimization$HR, na.rm = TRUE))

        points(TableNow[,"zscore"],
               rep(min(cutoff_optimization$HR_min95, na.rm = TRUE), length(TableNow[,"zscore"])),
               pch="|", col=rgb(37,37,37, 90, max=255))

        par(new = T)

        plot(x=cutoff_optimization$zscorecut,
             y=cutoff_optimization$logrank_p, type="l",
             ylim=c(max(cutoff_optimization$logrank_p, na.rm = TRUE),
                    min(cutoff_optimization$logrank_p, na.rm = TRUE)-0.1),
             xlim=c(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                    max(cutoff_optimization$zscorecut, na.rm = TRUE)),
             lwd=2, lty=1, xlab=NA, ylab=NA,
             axes=FALSE, cex=1.5, col=rgb(252,78,42, 150, max=255))
        axis(side = 4, lwd = 2, cex.axis=1.5, cex=1.5, col=rgb(252,78,42, 255, max=255), las =1)
        mtext(side = 4, line = 3.8, "p logrank", cex=1.2)

        min_p <- min(cutoff_optimization$logrank_p , na.rm = TRUE)

        # abline(h = min_p, lwd=2, lty=3, col=rgb(0,0,0, 255, max=255))
        text(best_z_table, min_p-.01, paste0("p logrank = ", round(min_p, 4)), cex=0.8)
        points(best_z_table, min_p)

        # par(xpd = TRUE)
        # text(x = par("usr")[2]+0.2, y = mean(par("usr")[3:4]), "p logrank", srt = 270, cex=1.3)
        dev.off()

        # plot HR log2
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "coxHR", paste0("Log2_with_confidence_interval_",
                                                          Name, ".png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "coxHR", paste0("Log2_with_confidence_interval_",
                                                          Name, ".svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        par(mar = c(5,5,3,5))
        plot(x=cutoff_optimization$zscorecut,
             y=log2(cutoff_optimization$HR), type="l",
             ylim=c(min(log2(cutoff_optimization$HR_min95), na.rm = TRUE),
                    max(log2(cutoff_optimization$HR_max95), na.rm = TRUE)),
             xlim=c(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                    max(cutoff_optimization$zscorecut, na.rm = TRUE)),
             lwd=2, lty=1, xlab="Expression z-score", ylab=expression('Log'[2]*'(Cox HR RFS)'), axes=F,
             col=rgb(8,88,158, 240, max=255))

        axis(side = 1, lwd = 2, cex.axis=1.5, cex=1.5,
             at=round(seq(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                          max(cutoff_optimization$zscorecut, na.rm = TRUE), by = 0.1), 2))
        axis(side = 2, lwd = 2, cex.axis=1.5, cex=1.5,
             col=rgb(8,88,158, 255, max=255), las =1)

        lines(x=cutoff_optimization$zscorecut,
              y=log2(cutoff_optimization$HR_min95),
              lwd=1, lty=3, col=rgb(19, 153, 24, 200, max=255))
        lines(x=cutoff_optimization$zscorecut,
              y=log2(cutoff_optimization$HR_max95),
              lwd=1, lty=3, col=rgb(19, 153, 24, 200, max=255))

        text(best_z_table, max(log2(cutoff_optimization$HR), na.rm = TRUE)+.105,
             paste0("z-score = ", round(best_z_table, 2)), cex=0.8)
        points(best_z_table, max(log2(cutoff_optimization$HR), na.rm = TRUE))

        points(TableNow[,"zscore"],
               rep(min(log2(cutoff_optimization$HR_min95), na.rm = TRUE), length(TableNow[,"zscore"])),
               pch="|", col=rgb(37,37,37, 90, max=255))
        par(new = T)

        plot(x=cutoff_optimization$zscorecut,
             y=cutoff_optimization$logrank_p, type="l",
             ylim=c(max(cutoff_optimization$logrank_p, na.rm = TRUE),
                    min(cutoff_optimization$logrank_p, na.rm = TRUE)),
             xlim=c(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                    max(cutoff_optimization$zscorecut, na.rm = TRUE)),
             lwd=2, lty=1, xlab=NA, ylab=NA,
             axes=F, cex=1.5, col=rgb(252,78,42, 150, max=255))
        axis(side = 4, lwd = 2, cex.axis=1.5, cex=1.5, col=rgb(252,78,42, 255, max=255), las =1)
        mtext(side = 4, line = 3.8, "p logrank", cex=1.2)

        min_p <- min(cutoff_optimization$logrank_p , na.rm = TRUE)

        text(best_z_table, min_p-.015, paste0("p logrank = ", round(min_p, 4)), cex=0.8)
        points(best_z_table, min_p)

        dev.off()

    }

    # hide outliers
    # plot
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "coxHR",
                                 paste0("hide_confidence_interval", Name, ".png")),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "coxHR",
                                 paste0("hide_confidence_interval", Name, ".svg")),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    par(mar = c(5,5,3,5))
    plot(x=cutoff_optimization$zscorecut,
         y=cutoff_optimization$HR, type="l",
         lwd=2, lty=1, xlab="Expression z-score", ylab="Cox HR RFS", axes=FALSE,
         col=rgb(8,88,158, alpha = 100, max=255))

    axis(side = 1, lwd = 2, cex.axis=1.5, cex=1.5,
         at=round(seq(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                      max(cutoff_optimization$zscorecut, na.rm = TRUE), by = 0.1), 2))
    axis(side = 2, lwd = 2, cex.axis=1.5, cex=1.5, col=rgb(8,88,158, 255, max=255), las = 1)
    # las = 1, at = round(seq(min(cutoff_optimization$HR, na.rm = TRUE),
    #                   max(cutoff_optimization$HR, na.rm = TRUE), by = 0.2), 1))

    text(best_z_table, max(cutoff_optimization$HR, na.rm = TRUE)+.12,
         paste0("z-score = ", round(best_z_table, 2)), cex=0.8)
    points(best_z_table, max(cutoff_optimization$HR, na.rm = TRUE))

    points(cutoff_optimization$zscorecut,
           rep(min(cutoff_optimization$HR, na.rm = TRUE), length(cutoff_optimization$zscorecut)),
           pch="|", col=rgb(37,37,37, 90, max=255))

    par(new = T)
    plot(x=cutoff_optimization$zscorecut,
         y=cutoff_optimization$logrank_p, type="l",
         ylim=c(max(cutoff_optimization$logrank_p, na.rm = TRUE), min(cutoff_optimization$logrank_p, na.rm = TRUE)),
         lwd=3, lty=1, xlab=NA, ylab=NA,
         axes=F, cex=1.5, col=rgb(252,78,42, alpha = 100, max=255))
    axis(side = 4, lwd = 2, cex.axis=1.5, cex=1.5, col=rgb(252,78,42, 255, max=255), las = 1)
    mtext(side = 4, line = 3.8, "p logrank", cex=1.2)

    min_p <- min(cutoff_optimization$logrank_p , na.rm = TRUE)
    # abline(h = min_p, lwd=2, lty=3, col=rgb(0,0,0, 255, max=255))
    text(best_z_table, min_p-.01, paste0("p logrank = ", round(min_p, 4)), cex=0.8)
    points(best_z_table, min_p)

    dev.off()

    # plot HR log2
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "coxHR",
                                 paste0("Log2_hide_confidence_interval_", Name, ".png")),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "coxHR",
                                 paste0("Log2_hide_confidence_interval_", Name, ".svg")),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    par(mar = c(5,5,3,5))
    plot(x=cutoff_optimization$zscorecut,
         y=log2(cutoff_optimization$HR), type="l",
         lwd=2, lty=1, xlab="Expression z-score", ylab=expression('Log'[2]*'(Cox HR RFS)'), axes=F,
         col=rgb(8,88,158, alpha = 100, max=255))

    axis(side = 1, lwd = 2, cex.axis=1.5, cex=1.5,
         at=round(seq(min(cutoff_optimization$zscorecut, na.rm = TRUE),
                      max(cutoff_optimization$zscorecut, na.rm = TRUE), by = 0.1), 2))
    axis(side = 2, lwd = 2, cex.axis=1.5, cex=1.5, col=rgb(8,88,158, 255, max=255), las = 1)
    # las = 1, at = round(seq(min(log2(cutoff_optimization$HR), na.rm = TRUE),
    #                         max(log2(cutoff_optimization$HR), na.rm = TRUE), by = 0.2), 1))

    text(best_z_table, max(log2(cutoff_optimization$HR), na.rm = TRUE)+.105,
         paste0("z-score = ", round(best_z_table, 2)), cex=0.8)
    points(best_z_table, max(log2(cutoff_optimization$HR), na.rm = TRUE))

    points(cutoff_optimization$zscorecut,
           rep(min(log2(cutoff_optimization$HR), na.rm = TRUE), length(cutoff_optimization$zscorecut)),
           pch="|", col=rgb(37,37,37, 90, max=255))

    par(new = T)

    plot(x=cutoff_optimization$zscorecut,
         y=cutoff_optimization$logrank_p, type="l",
         ylim=c(max(cutoff_optimization$logrank_p, na.rm = TRUE), min(cutoff_optimization$logrank_p, na.rm = TRUE)),
         lwd=3, lty=1, xlab=NA, ylab=NA,
         axes=F, cex=1.5, col=rgb(252,78,42, alpha = 100, max=255))
    axis(side = 4, lwd = 2, cex.axis=1.5, cex=1.5, col=rgb(252,78,42, 255, max=255), las = 1)
    mtext(side = 4, line = 3.8, "p logrank", cex=1.2)

    min_p <- min(cutoff_optimization$logrank_p , na.rm = TRUE)
    # abline(h = min_p, lwd=2, lty=3, col=rgb(0,0,0, 255, max=255))
    text(best_z_table, min_p-.01, paste0("p logrank = ", round(min_p, 4)), cex=0.8)
    points(best_z_table, min_p)

    dev.off()

    # Save best z table
    # save optim table
    write.csv(best_z_table, file=file.path(DIR, "coxHR",
                                           paste0("Best_z_Table", ".csv")))

    # part3 ####
    # KAPLAN MEYER ACCORDINGLY TO LEVEL CUTOFF

    # add unknown to everyone
    framesList[, "classification"] <- "unknown"
    # add low expressing tag
    framesList[which(framesList[,"zscore"] <= best_z_table), "classification"] <- "low"
    # add high expressing tag
    framesList[which(framesList[,"zscore"] > best_z_table), "classification"] <- "high"
    # create full classification
    framesList[,"full_classification"] <- framesList[,"classification"]

    # copy lines to move
    linestomove <- which(framesList[,"classification"] == "unknown")

    framesList[linestomove, "full_classification"] <- framesList[linestomove,"TissueTypeSimple"]
    # factor it
    framesList[, "full_classification"] <- factor(framesList[, "full_classification"],
                                                  levels = c("low", "high"))

    # remove used variable
    remove(linestomove)


    #Plotting
    ### Comparison normal - low - high levels by log2 RSEM UQ level

    # define the formula now
    formulaNow <- formula(framesList[, "log2p1"] ~ framesList[, "full_classification"])

    if (nrow(framesList) > 0) {
        # Pirate plotting - Basic
        # plot HR log2
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "kaplan_maier",
                                     paste0(Name, "_PiratePlot_log2RSEM_all_conditions.png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "kaplan_maier",
                                     paste0(Name, "_PiratePlot_log2RSEM_all_conditions.svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        #bottom, left, top and right
        par(mar=c(3,4.8,4,2))
        yarrr::pirateplot(formula = formulaNow,
                          data=framesList, gl.lwd = 0,
                          ylab = expression('Log'[2]*'(Expression + 1)'),
                          # main = paste(Name, "_", names(framesList), sep=""),
                          inf.method="iqr", cex.lab=1.2, cex.axis=1.8,
                          inf.disp="rect", jitter.val = 0.08, theme=2, cex.names=1.5,
                          pal="pony", avg.line.fun = median, point.cex = 1.3,
                          inf.f.col=rgb(200,200,200, 255, max=255), xlab = "",
                          inf.b.col=rgb(120,120,120, 255, max=255))
        dev.off()
    }

    # ANOVA followed by tukey
    sink(file=file.path(DIR, "kaplan_maier", paste0(Name, "_PiratePlot_log2RSEM",
                                                    "_all_conditions_anova_tukey.txt")))

    fit <- aov(formulaNow,
               data=framesList)

    cat("********************\n")
    cat(paste("A - ANOVA","\n", sep=""))
    print(summary(fit))

    cat("********************\n")
    cat(paste("B - Shapiro Test for Normality","\n", sep=""))
    cat("\n")
    cat(paste("*B1 - all values test","\n", sep=""))
    print(shapiro.test(framesList[,"log2p1"]))
    cat(paste("*B2 - fit residuals test","\n", sep=""))
    print(shapiro.test(residuals(fit)))

    cat("\n")
    cat("********************\n")
    cat(paste("C - Tukey Honest Significant Differences post-test","\n", sep=""))
    print(TukeyHSD(fit))

    sink()

    # Kruskal-wallis followd by Wilcoxon test (for nono-parametric assumption)
    sink(file=file.path(DIR, "kaplan_maier", paste0(Name, "_PiratePlot_log2RSEM",
                                                    "_all_conditions_Kruskal_Wilcoxon.txt")))
    cat("********************\n")
    cat(paste("A - Kruskal-Wallis test","\n", sep=""))

    print(kruskal.test(data=framesList,
                       formulaNow))

    cat("********************\n")
    cat(paste("B - Non-adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
    print(pairwise.wilcox.test(x=framesList[,"log2p1"],
                               g=framesList[,"full_classification"],
                               p.adjust.method = "none",
                               paired = FALSE,
                               exact = TRUE,
                               correct = TRUE,
                               alternative="two.sided"))

    cat("\n")
    cat("********************\n")
    cat(paste("C - Benjamini FDR adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
    print(pairwise.wilcox.test(x=framesList[,"log2p1"],
                               g=framesList[,"full_classification"],
                               p.adjust.method = "fdr",
                               paired = FALSE,
                               exact = TRUE,
                               correct = TRUE,
                               alternative="two.sided"))

    sink()

    # Clean formula now
    remove(formulaNow)


    # Create table for patient numbers
    PatientAmountTable <- matrix(ncol=4, nrow=8)
    colnames(PatientAmountTable) <- c("type", "classification", "amount", "percent")
    PatientAmountTable <- data.frame(PatientAmountTable)

    # use previously defined markers

    # Copy PatientAmountTable
    PatientAmountTable_now <- PatientAmountTable

    # OVERALL SURVIVAL ANALYSIS
    # Kaplan meyer using best scenario
    KaplanTableNow <- framesList[,c(Name,
                                               "classification",
                                               "final_vital_status_times",
                                               "final_vital_status",
                                               "final_rfs_status",
                                               "final_rfs_status_times",
                                               "final_dmfs_status",
                                               "final_dmfs_status_times")]

    # Kaplan
    colnames(KaplanTableNow)[2] <- "classification"

    # Collect just rows with classification
    KaplanTableNow <- KaplanTableNow[c(which(KaplanTableNow[,"classification"] == "low"),
                                       which(KaplanTableNow[,"classification"] == "high")),]


    # Factorize classification
    KaplanTableNow$classification_factorized <- factor(KaplanTableNow$classification, levels=c("low", "high"))

    ############# COLLECT DATA

    # Fill table with patient amount - ALL
    rowsNow <- c(1,2)
    PatientAmountTable_now[rowsNow,"type"] <- "All"
    PatientAmountTable_now[rowsNow,"classification"] <- c("low", "high")
    PatientAmountTable_now[rowsNow[1],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["low"])
    PatientAmountTable_now[rowsNow[2],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["high"])
    PatientAmountTable_now[rowsNow[1],"percent"] <- PatientAmountTable_now[rowsNow[1], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])
    PatientAmountTable_now[rowsNow[2],"percent"] <- PatientAmountTable_now[rowsNow[2], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])


    # Retrieve non-NA entries for final vital status
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_vital_status == "unknown"),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_vital_status == "vital_status"),]
    KaplanTableNow <- KaplanTableNow[!is.na(KaplanTableNow$final_vital_status),]
    KaplanTableNow <- KaplanTableNow[!is.nan(KaplanTableNow$final_vital_status),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_vital_status_times == "unknown"),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_vital_status_times == "day"),]
    KaplanTableNow <- KaplanTableNow[!is.na(KaplanTableNow$final_vital_status_times),]
    KaplanTableNow <- KaplanTableNow[!is.nan(KaplanTableNow$final_vital_status_times),]
    KaplanTableNow$final_vital_status_times <- as.numeric(KaplanTableNow$final_vital_status_times)

    # Overall survival
    KaplanTableNow$final_vital_status <- as.integer(gsub("Dead", 1, gsub("Alive", 0, KaplanTableNow$final_vital_status)))


    # Fill table with patient amount - OSurv
    rowsNow <- c(3,4)
    PatientAmountTable_now[rowsNow,"type"] <- "OSurv"
    PatientAmountTable_now[rowsNow,"classification"] <- c("low", "high")
    PatientAmountTable_now[rowsNow[1],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["low"])
    PatientAmountTable_now[rowsNow[2],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["high"])

    PatientAmountTable_now[rowsNow[1],"percent"] <- PatientAmountTable_now[rowsNow[1], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])
    PatientAmountTable_now[rowsNow[2],"percent"] <- PatientAmountTable_now[rowsNow[2], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])


    # Perform the fitting
    remove(fit)
    fit <- survival::survfit(survival::Surv(final_vital_status_times, final_vital_status) ~ classification_factorized,
                   data = KaplanTableNow)

    # Report p-value
    sink(file=file.path(DIR, "kaplan_maier",
                        paste0(Name, "Kaplan_overall_survival", ".txt")))
    print(summary(survival::coxph(survival::Surv(final_vital_status_times, final_vital_status) ~ classification_factorized,
                        data = KaplanTableNow)))
    sink()


    # Plot data
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name, "Kaplan_overall_survival.png")),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name, "Kaplan_overall_survival.svg")),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    par(mar=c(3,3,3,6))
    print(survminer::ggsurvplot(fit, data = KaplanTableNow, risk.table = TRUE,
                     linetype=c(1,1), conf.int = TRUE,
                     #palette = "grey",
                     palette =c("#4575b4", "#d73027"),
                     ggtheme = ggplot2::theme_bw(),
                     size = 2,
                     #xscale="d_m",
                     font.main = 30,
                     font.x =  26,
                     font.y = 26,
                     font.tickslab = 24))
    dev.off()

    # Create 5yr status
    KaplanTableNow$yr5_status <- "unknown"
    KaplanTableNow$yr5_status[KaplanTableNow$final_vital_status_times >= (365*5)] <- "alive_5yr"
    KaplanTableNow$yr5_status[(KaplanTableNow$final_vital_status_times < (365*5) & KaplanTableNow$final_vital_status == 1)] <- "dead_5yr"

    KaplanTableNow$yr5_status <- as.factor(KaplanTableNow$yr5_status)

    # define the formula now
    formulaNow <- formula(log2(get(Name)+1) ~ yr5_status)

    message("Pirate plotting \n")

    if (nrow(framesList) > 0) {
        # Pirate plotting - Basic
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "kaplan_maier",
                                     paste0(Name,"PiratePlot_log2RSEM__5yr_overallsurvival.png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "kaplan_maier",
                                     paste0(Name,"PiratePlot_log2RSEM__5yr_overallsurvival.svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        #bottom, left, top and right
        par(mar=c(3,4.8,4,2))
        yarrr::pirateplot(formula = formulaNow,
                          data=KaplanTableNow,
                          ylab = expression('Log'[2]*'(Expression + 1)'),
                          xlab = "",
                          # main = paste(Name, "_", names(framesList), sep=""),
                          inf.method="iqr", cex.lab=1.2, cex.axis=1.8,
                          inf.disp="rect", jitter.val = 0.08, theme=2, cex.names=1.5,
                          pal="pony", avg.line.fun = median, point.cex = 1.3,
                          inf.f.col=rgb(200,200,200, 255, max=255),
                          inf.b.col=rgb(120,120,120, 255, max=255))
        #sortx = "alphabetical")
        dev.off()
    }


    # ANOVA followed by tukey
    sink(file=file.path(DIR, "kaplan_maier",
                        paste0(Name,"PiratePlot_log2RSEM__5yr_overallsurvival_anova_tukey.txt")))

    fit <- aov(formulaNow,
               data=KaplanTableNow)

    cat("********************\n")
    cat(paste("A - ANOVA","\n", sep=""))
    print(summary(fit))


    cat("********************\n")
    cat(paste("B - Shapiro Test for Normality","\n", sep=""))
    cat("\n")
    cat(paste("*B1 - all values test","\n", sep=""))
    print(shapiro.test(log2(KaplanTableNow[,Name]+1)))
    cat(paste("*B2 - fit residuals test","\n", sep=""))
    print(shapiro.test(residuals(fit)))

    cat("\n")
    cat("********************\n")
    cat(paste("C - Tukey Honest Significant Differences post-test","\n", sep=""))
    print(TukeyHSD(fit))

    sink()

    # Kruskal-wallis followd by Wilcoxon test (for nono-parametric assumption)
    sink(file=file.path(DIR, "kaplan_maier",
                        paste0(Name,"PiratePlot_log2RSEM__5yr_overallsurvival_kruskal_wilcoxon.txt")))
    cat("********************\n")
    cat(paste("A - Kruskal-Wallis test","\n", sep=""))

    print(kruskal.test(data=KaplanTableNow,
                       formulaNow))

    cat("********************\n")
    cat(paste("B - Non-adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
    print(pairwise.wilcox.test(x=log2(KaplanTableNow[,Name]+1),
                               g=KaplanTableNow[,"yr5_status"],
                               p.adjust.method = "none",
                               paired = FALSE,
                               exact = TRUE,
                               correct = TRUE,
                               alternative="two.sided"))

    cat("\n")
    cat("********************\n")
    cat(paste("C - Benjamini FDR adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
    print(pairwise.wilcox.test(x=log2(KaplanTableNow[,Name]+1),
                               g=KaplanTableNow[,"yr5_status"],
                               p.adjust.method = "fdr",
                               paired = FALSE,
                               exact = TRUE,
                               correct = TRUE,
                               alternative="two.sided"))

    sink()

    # Clean formula now
    remove(formulaNow)

    # relapse
    # Kaplan meyer using best scenario
    KaplanTableNow <- framesList[,c(Name,
                                               "classification",
                                               "final_vital_status_times",
                                               "final_vital_status",
                                               "final_rfs_status",
                                               "final_rfs_status_times",
                                               "final_dmfs_status",
                                               "final_dmfs_status_times")]

    # Kaplan
    colnames(KaplanTableNow)[2] <- "classification"

    # Collect just rows with classification
    KaplanTableNow <- KaplanTableNow[c(which(KaplanTableNow[,"classification"] == "low"),
                                       which(KaplanTableNow[,"classification"] == "high")),]


    # Factorize classification
    KaplanTableNow$classification_factorized <- factor(KaplanTableNow$classification, levels=c("low", "high"))

    # Retrieve non-NA entries for final vital status
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_rfs_status == "unknown"),]
    KaplanTableNow <- KaplanTableNow[!is.na(KaplanTableNow$final_rfs_status),]
    KaplanTableNow <- KaplanTableNow[!is.nan(KaplanTableNow$final_rfs_status),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_rfs_status_times == "unknown"),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_rfs_status_times == "dead_other_cause_before_relapse"),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_rfs_status_times == "day"),]
    KaplanTableNow <- KaplanTableNow[!is.na(KaplanTableNow$final_rfs_status_times),]
    KaplanTableNow <- KaplanTableNow[!is.nan(KaplanTableNow$final_rfs_status_times),]
    KaplanTableNow <- KaplanTableNow[!is.infinite(KaplanTableNow$final_rfs_status_times),]
    KaplanTableNow$final_rfs_status_times <- as.numeric(KaplanTableNow$final_rfs_status_times)

    # Overall survival
    KaplanTableNow$final_rfs_status <- as.integer(gsub("relapse", 1, gsub("censored", 0, KaplanTableNow$final_rfs_status)))

    # Fill table with patient amount - Relapse
    rowsNow <- c(5,6)
    PatientAmountTable_now[rowsNow,"type"] <- "Relapse"
    PatientAmountTable_now[rowsNow,"classification"] <- c("low", "high")
    PatientAmountTable_now[rowsNow[1],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["low"])
    PatientAmountTable_now[rowsNow[2],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["high"])

    PatientAmountTable_now[rowsNow[1],"percent"] <- PatientAmountTable_now[rowsNow[1], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])
    PatientAmountTable_now[rowsNow[2],"percent"] <- PatientAmountTable_now[rowsNow[2], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])


    # Report p-value
    sink(file=file.path(DIR, "kaplan_maier",
                        paste0(Name,"Kaplan_rfs_survival", ".txt")))
    print(summary(survival::coxph(survival::Surv(final_rfs_status_times, final_rfs_status) ~ classification_factorized,
                        data = KaplanTableNow)))
    sink()


    # Perform the fitting
    remove(fit)
    fit <- survival::survfit(survival::Surv(final_rfs_status_times, final_rfs_status) ~ classification_factorized,
                   data = KaplanTableNow)

    # Plot data
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name,"Kaplan_rfs_survival", ".png")),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name,"Kaplan_rfs_survival", ".svg")),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    par(mar=c(3,3,3,6))
    print(survminer::ggsurvplot(fit, data = KaplanTableNow, risk.table = TRUE,
                     linetype=c(1,1), conf.int = TRUE,
                     #palette = "grey",
                     palette =c("#4575b4", "#d73027"),
                     ggtheme = ggplot2::theme_bw(),
                     size = 2,
                     #xscale="d_m",
                     font.main = 30,
                     font.x =  26,
                     font.y = 26,
                     font.tickslab = 24))
    dev.off()

    # Create 5yr status
    KaplanTableNow$yr5_status <- "unknown"
    KaplanTableNow$yr5_status[KaplanTableNow$final_rfs_status_times >= (365*5)] <- "nonrelapse_5yr"
    KaplanTableNow$yr5_status[(KaplanTableNow$final_rfs_status_times < (365*5) & KaplanTableNow$final_rfs_status == 1)] <- "relapse_5yr"

    #
    KaplanTableNow$yr5_status <- as.factor(KaplanTableNow$yr5_status)

    # define the formula now
    formulaNow <- formula(log2(get(Name)+1) ~ yr5_status)

    message("Pirate plotting \n")

    if (nrow(framesList) > 0) {
        # Pirate plotting - Basic
        if (tolower(image_format) == "png") {
            png(filename = file.path(DIR, "kaplan_maier",
                                     paste0(Name,"PiratePlot_log2RSEM__5yr_rfs.png")),
                width = Width, height = Height, res = Res, units = Unit)
        } else if (tolower(image_format) == "svg") {
            svg(filename = file.path(DIR, "kaplan_maier",
                                     paste0(Name,"PiratePlot_log2RSEM__5yr_rfs.svg")),
                width = Width, height = Height, onefile = TRUE)
        } else {
            stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
        }
        #bottom, left, top and right
        par(mar=c(3,4.8,4,2))
        yarrr::pirateplot(formula = formulaNow,
                          data=KaplanTableNow,
                          xlab = "",
                          ylab = expression('Log'[2]*'(Expression + 1)'),
                          # main = paste(Name, "_", names(framesList), sep=""),
                          inf.method="iqr", cex.lab=1.2, cex.axis=1.8,
                          inf.disp="rect", jitter.val = 0.08, theme=2, cex.names=1.5,
                          pal="pony", avg.line.fun = median,
                          inf.f.col=rgb(200,200,200, 255, max=255),
                          inf.b.col=rgb(120,120,120, 255, max=255))
        #sortx = "alphabetical")
        dev.off()
    }

    # ANOVA followed by tukey
    sink(file=file.path(DIR, "kaplan_maier",
                        paste0(Name,"PiratePlot_log2RSEM__5yr_rfs_anova_tukey.txt")))

    fit <- aov(formulaNow,
               data=KaplanTableNow)

    cat("********************\n")
    cat(paste("A - ANOVA","\n", sep=""))
    print(summary(fit))


    cat("********************\n")
    cat(paste("B - Shapiro Test for Normality","\n", sep=""))
    cat("\n")
    cat(paste("*B1 - all values test","\n", sep=""))
    print(shapiro.test(log2(KaplanTableNow[,Name]+1)))
    cat(paste("*B2 - fit residuals test","\n", sep=""))
    print(shapiro.test(residuals(fit)))

    cat("\n")
    cat("********************\n")
    cat(paste("C - Tukey Honest Significant Differences post-test","\n", sep=""))
    print(TukeyHSD(fit))

    sink()

    # Kruskal-wallis followd by Wilcoxon test (for nono-parametric assumption)
    sink(file=file.path(DIR, "kaplan_maier",
                        paste0(Name,"_PiratePlot_log2RSEM__5yr_rfs_kruskal_wilcoxon.txt")))
    cat("********************\n")
    cat(paste("A - Kruskal-Wallis test","\n", sep=""))

    print(kruskal.test(data=KaplanTableNow,
                       formulaNow))

    cat("********************\n")
    cat(paste("B - Non-adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
    print(pairwise.wilcox.test(x=log2(KaplanTableNow[,Name]+1),
                               g=KaplanTableNow[,"yr5_status"],
                               p.adjust.method = "none",
                               paired = FALSE,
                               exact = TRUE,
                               correct = TRUE,
                               alternative="two.sided"))

    cat("\n")
    cat("********************\n")
    cat(paste("C - Benjamini FDR adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
    print(pairwise.wilcox.test(x=log2(KaplanTableNow[,Name]+1),
                               g=KaplanTableNow[,"yr5_status"],
                               p.adjust.method = "fdr",
                               paired = FALSE,
                               exact = TRUE,
                               correct = TRUE,
                               alternative="two.sided"))

    sink()

    # Clean formula now
    remove(formulaNow)


    # dmfs
    # Kaplan meyer using best scenario
    KaplanTableNow <- framesList[,c(Name,
                                               "classification",
                                               "final_vital_status_times",
                                               "final_vital_status",
                                               "final_rfs_status",
                                               "final_rfs_status_times",
                                               "final_dmfs_status",
                                               "final_dmfs_status_times")]

    # Kaplan
    colnames(KaplanTableNow)[2] <- "classification"

    # Collect just rows with classification
    KaplanTableNow <- KaplanTableNow[c(which(KaplanTableNow[,"classification"] == "low"),
                                       which(KaplanTableNow[,"classification"] == "high")),]

    # Factorize classification
    KaplanTableNow$classification_factorized <- factor(KaplanTableNow$classification, levels=c("low", "high"))

    # Retrieve non-NA entries for final vital status
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_dmfs_status == "unknown"),]
    KaplanTableNow <- KaplanTableNow[!is.na(KaplanTableNow$final_dmfs_status),]
    KaplanTableNow <- KaplanTableNow[!is.nan(KaplanTableNow$final_dmfs_status),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_dmfs_status_times == "unknown"),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_dmfs_status_times == "dead_other_cause_before_relapse"),]
    KaplanTableNow <- KaplanTableNow[(!KaplanTableNow$final_dmfs_status_times == "day"),]
    KaplanTableNow <- KaplanTableNow[!is.na(KaplanTableNow$final_dmfs_status_times),]
    KaplanTableNow <- KaplanTableNow[!is.nan(KaplanTableNow$final_dmfs_status_times),]
    KaplanTableNow <- KaplanTableNow[!is.infinite(KaplanTableNow$final_dmfs_status_times),]
    KaplanTableNow$final_dmfs_status_times <- as.numeric(KaplanTableNow$final_dmfs_status_times)

    # Overall survival
    KaplanTableNow$final_dmfs_status <- as.integer(gsub("metastasis", 1, gsub("censored", 0, KaplanTableNow$final_dmfs_status)))


    # Fill table with patient amount - DMFS
    rowsNow <- c(7,8)
    PatientAmountTable_now[rowsNow,"type"] <- "DMFS"
    PatientAmountTable_now[rowsNow,"classification"] <- c("low", "high")
    PatientAmountTable_now[rowsNow[1],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["low"])
    PatientAmountTable_now[rowsNow[2],"amount"] <- as.numeric(table(KaplanTableNow[,"classification"])["high"])

    PatientAmountTable_now[rowsNow[1],"percent"] <- PatientAmountTable_now[rowsNow[1], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])
    PatientAmountTable_now[rowsNow[2],"percent"] <- PatientAmountTable_now[rowsNow[2], "amount"] / (PatientAmountTable_now[rowsNow[1], "amount"] + PatientAmountTable_now[rowsNow[2], "amount"])


    # Skip fitting

    # Perform the fitting
    remove(fit)
    fit <- survival::survfit(survival::Surv(final_dmfs_status_times, final_dmfs_status) ~ classification_factorized,
                   data = KaplanTableNow)

    #
    if(sum(KaplanTableNow$final_dmfs_status) > 0){
        # Report p-value
        sink(file=file.path(DIR, "kaplan_maier",
                            paste0(Name,"_Kaplan_dmfs_survival", ".txt")))
        print(summary(survival::coxph(survival::Surv(final_dmfs_status_times, final_dmfs_status) ~ classification_factorized,
                            data = KaplanTableNow)))
        sink()
    }

    # Plot data
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name,"Kaplan_dmfs_survival", ".png")),
            width = Width, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name,"Kaplan_dmfs_survival", ".svg")),
            width = Width, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    par(mar=c(3,3,3,6))
    print(survminer::ggsurvplot(fit, data = KaplanTableNow, risk.table = TRUE,
                     linetype=c(1,1), conf.int = TRUE,
                     #palette = "grey",
                     palette =c("#4575b4", "#d73027"),
                     ggtheme = ggplot2::theme_bw(),
                     size = 2,
                     #xscale="d_m",
                     font.main = 30,
                     font.x =  26,
                     font.y = 26,
                     font.tickslab = 24))
    dev.off()

    # Create 5yr status
    KaplanTableNow$yr5_status <- "unknown"
    KaplanTableNow$yr5_status[KaplanTableNow$final_dmfs_status_times >= (365*5)] <- "nonmetas_5yr"
    KaplanTableNow$yr5_status[(KaplanTableNow$final_dmfs_status_times < (365*5) & KaplanTableNow$final_dmfs_status == 1)] <- "metas_5yr"

    #
    KaplanTableNow$yr5_status <- as.factor(KaplanTableNow$yr5_status)

    # define the formula now
    formulaNow <- formula(log2(get(Name)+1) ~ yr5_status)

    # if (nrow(framesList) > 0) {
    #     # Pirate plotting - Basic
    #     if (tolower(image_format) == "png") {
    #         png(filename = file.path(DIR, "kaplan_maier",
    #                                  paste0(Name,"_PiratePlot_log2RSEM__5yr_dmfs.png")),
    #             width = Width, height = Height, res = Res, units = Unit)
    #     } else if (tolower(image_format) == "svg") {
    #         svg(filename = file.path(DIR, "kaplan_maier",
    #                                  paste0(Name,"_PiratePlot_log2RSEM__5yr_dmfs.svg")),
    #             width = Width, height = Height, onefile = TRUE)
    #     } else {
    #         stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    #     }
    #     #bottom, left, top and right
    #     par(mar=c(3,4.8,4,2))
    #     yarrr::pirateplot(formula = formulaNow,
    #                       data=KaplanTableNow,
    #                       xlab = "",
    #                       ylab = expression('Log'[2]*'(Expression + 1)'),
    #                       # main = paste(Name, "_", names(framesList), sep=""),
    #                       inf.method="iqr", cex.lab=1.2, cex.axis=1.8,
    #                       inf.disp="rect", jitter.val = 0.08, theme=2, cex.names=1.5,
    #                       pal="pony", avg.line.fun = median,
    #                       inf.f.col=rgb(200,200,200, 255, max=255),
    #                       inf.b.col=rgb(120,120,120, 255, max=255))
    #     #sortx = "alphabetical")
    #     dev.off()
    # }

    # check to skip one conditions
    if(length(unique(KaplanTableNow$yr5_status))>=2){
        # ANOVA followed by tukey
        sink(file=file.path(DIR, "kaplan_maier",
                            paste0(Name,"_PiratePlot_log2RSEM__5yr_dmfs_anova_tukey.txt")))

        fit <- aov(formulaNow,
                   data=KaplanTableNow)

        cat("********************\n")
        cat(paste("A - ANOVA","\n", sep=""))
        print(summary(fit))


        cat("********************\n")
        cat(paste("B - Shapiro Test for Normality","\n", sep=""))
        cat("\n")
        cat(paste("*B1 - all values test","\n", sep=""))
        print(shapiro.test(log2(KaplanTableNow[,Name]+1)))
        cat(paste("*B2 - fit residuals test","\n", sep=""))
        print(shapiro.test(residuals(fit)))

        cat("\n")
        cat("********************\n")
        cat(paste("C - Tukey Honest Significant Differences post-test","\n", sep=""))
        print(TukeyHSD(fit))

        sink()

        # if("yes"=="no"){
        # Kruskal-wallis followd by Wilcoxon test (for nono-parametric assumption)
        sink(file=file.path(DIR, "kaplan_maier",
                            paste0(Name,"_PiratePlot_log2RSEM__5yr_dmfs_kruskal_wilcoxon.txt")))
        cat("********************\n")
        cat(paste("A - Kruskal-Wallis test","\n", sep=""))

        print(kruskal.test(data=KaplanTableNow,
                           formulaNow))

        cat("********************\n")
        cat(paste("B - Non-adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
        print(pairwise.wilcox.test(x=log2(KaplanTableNow[,Name]+1),
                                   g=KaplanTableNow[,"yr5_status"],
                                   p.adjust.method = "none",
                                   paired = FALSE,
                                   exact = TRUE,
                                   correct = TRUE,
                                   alternative="two.sided"))

        cat("\n")
        cat("********************\n")
        cat(paste("C - Benjamini FDR adjusted Wilcoxon (Mann-Whitney U)","\n", sep=""))
        print(pairwise.wilcox.test(x=log2(KaplanTableNow[,Name]+1),
                                   g=KaplanTableNow[,"yr5_status"],
                                   p.adjust.method = "fdr",
                                   paired = FALSE,
                                   exact = TRUE,
                                   correct = TRUE,
                                   alternative="two.sided"))

        sink()
    }


    # Clean formula now
    remove(formulaNow)
    remove(KaplanTableNow)

    ### patient amount plot
    #Fix factor
    PatientAmountTable_now$type <- factor(PatientAmountTable_now$type, levels=c("All", "OSurv", "Relapse", "DMFS"))
    PatientAmountTable_now$classification <- factor(PatientAmountTable_now$classification, levels=c("high", "low"))

    # Add position for plotting
    PatientAmountTable_now$pos <- 0
    for (type in levels(PatientAmountTable_now$type)) {
        lines <- PatientAmountTable_now$type == type
        numbers <- PatientAmountTable_now$amount[lines]
        PatientAmountTable_now$pos[lines] <- cumsum(numbers) - (0.5 * numbers)
    }

    # Add percents
    PatientAmountTable_now$percent <- (PatientAmountTable_now$percent * 100)
    PatientAmountTable_now$percent <- round(PatientAmountTable_now$percent, 0)


    # Pirate plotting - Basic
    if (tolower(image_format) == "png") {
        png(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name,"_PatientAmount.png")),
            width = Width/2, height = Height, res = Res, units = Unit)
    } else if (tolower(image_format) == "svg") {
        svg(filename = file.path(DIR, "kaplan_maier",
                                 paste0(Name,"_PatientAmount.svg")),
            width = Width/2, height = Height, onefile = TRUE)
    } else {
        stop(message("Please, Insert a valid image_format! ('png' or 'svg')"))
    }
    # Start ggplot2 for patient amount
    p4 <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_bar(ggplot2::aes(y = amount, x = type, fill = classification), data = PatientAmountTable_now,
                 stat="identity") +
        ggplot2::geom_text(data=PatientAmountTable_now, ggplot2::aes(x = type, y = pos, label = paste0(percent, "%")),
                  size=3.6) +
        ggplot2::scale_fill_manual(values=c(rgb(215,48,39,220,max=255), rgb(69,117,200,180,max=255))) +
        ggplot2::theme(legend.position="bottom", legend.direction="horizontal",
              legend.title = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.text.x=ggplot2::element_text(colour = "black", size=14),
              axis.text.y=ggplot2::element_text(colour = "black", size=13)
        )
    print(p4)
    dev.off()


    # Save table with patient amounts
    write.csv(PatientAmountTable_now, file.path(DIR, "kaplan_maier", paste0("PatientAmount_", Name, ".csv")))


    # CLINICAL CONTINGENCY TABLE

    #create empty table
    ClinicoPathologicalTable <- matrix()

    # pathologic_stage
    clinicalNameNow <- "pathologic_stage"
    # Collect data
    clinicoNow <- framesList[,c("full_classification",
                                           clinicalNameNow)]

    # summary stage
    clinicoNow$summary_simple <- "remove"

    "I" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC")]
    "II" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")]
    "III" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")]
    "IV" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC")]

    # Collect just tumors low high
    clinicoNow <- clinicoNow[clinicoNow[,"full_classification"] %in% c("low", "high"),]

    # remove non classified lines
    clinicoNow <- clinicoNow[!clinicoNow$summary_simple == "remove",]

    # Drop unused levels
    clinicoNow[,"full_classification"] <-
        droplevels(clinicoNow[,"full_classification"])

    #
    # calculate N sizes
    n_now <- table(clinicoNow[,"full_classification"])
    n_prop_now <- round(prop.table(n_now)*100, 1)


    # check
    if(dim(clinicoNow)[1] > 2){

        # Create table
        contingencyNow <- table(clinicoNow$summary_simple,
                                clinicoNow[,"full_classification"])
        contingencyNow_percent <- round(prop.table(contingencyNow)*100, 1)

        # join information in order
        # as total followed by percents
        ClinicoPathologicalTable_now <- rbind(c(clinicalNameNow, "", "", "", round(chisq.test(contingencyNow)$p.value, 4)),
                                              c("size", n_now[1], n_prop_now[1], n_now[2], n_prop_now[2]),
                                              cbind(rownames(contingencyNow),
                                                    unname(contingencyNow[,1,drop=TRUE]),
                                                    unname(contingencyNow_percent[,1,drop=TRUE]),
                                                    unname(contingencyNow[,2,drop=TRUE]),
                                                    unname(contingencyNow_percent[,2,drop=TRUE]))
        )

        colnames(ClinicoPathologicalTable_now) <- c("type", colnames(contingencyNow)[1], colnames(contingencyNow_percent)[1], colnames(contingencyNow)[2], colnames(contingencyNow_percent)[2])

        # Add or join
        #ClinicoPathologicalTable <- rbind(ClinicoPathologicalTable, ClinicoPathologicalTable_now)
        ClinicoPathologicalTable <- ClinicoPathologicalTable_now


    } else {
        ClinicoPathologicalTable <- matrix(ncol=5, nrow=0)
        colnames(ClinicoPathologicalTable) <- c("type", "high", "high", "low", "low")
    }


    # pathologic_T
    clinicalNameNow <- "pathologic_T"
    # Collect data
    clinicoNow <- framesList[,c("full_classification",
                                           clinicalNameNow)]

    # summary stage
    clinicoNow$summary_simple <- "remove"

    "T1" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("T1", "T1a", "T1b", "T1c", "T1d")]
    "T2" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("T2", "T2a", "T2b", "T2c", "T2d")]
    "T3" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("T3", "T3a", "T3b", "T3c", "T3d")]
    "T4" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("T4", "T4a", "T4b", "T4c", "T4d")]

    # Collect just tumors low high
    clinicoNow <- clinicoNow[clinicoNow[,"full_classification"] %in% c("low", "high"),]

    # remove non classified lines
    clinicoNow <- clinicoNow[!clinicoNow$summary_simple == "remove",]

    # Drop unused levels
    clinicoNow[,"full_classification"] <-
        droplevels(clinicoNow[,"full_classification"])

    #
    #
    # calculate N sizes
    n_now <- table(clinicoNow[,"full_classification"])
    n_prop_now <- round(prop.table(n_now)*100, 1)


    # check
    if(dim(clinicoNow)[1] > 2){

        # Create table
        contingencyNow <- table(clinicoNow$summary_simple,
                                clinicoNow[,"full_classification"])
        contingencyNow_percent <- round(prop.table(contingencyNow)*100, 1)

        # join information in order
        # as total followed by percents
        ClinicoPathologicalTable_now <- rbind(c(clinicalNameNow, "", "", "", round(chisq.test(contingencyNow)$p.value, 4)),
                                              c("size", n_now[1], n_prop_now[1], n_now[2], n_prop_now[2]),
                                              cbind(rownames(contingencyNow),
                                                    unname(contingencyNow[,1,drop=TRUE]),
                                                    unname(contingencyNow_percent[,1,drop=TRUE]),
                                                    unname(contingencyNow[,2,drop=TRUE]),
                                                    unname(contingencyNow_percent[,2,drop=TRUE]))
        )

        colnames(ClinicoPathologicalTable_now) <- c("type", colnames(contingencyNow)[1], colnames(contingencyNow_percent)[1], colnames(contingencyNow)[2], colnames(contingencyNow_percent)[2])

        # Add or join
        ClinicoPathologicalTable <- rbind(ClinicoPathologicalTable, ClinicoPathologicalTable_now)
    }

    # pathologic_N
    clinicalNameNow <- "pathologic_N"
    # Collect data
    clinicoNow <- framesList[,c("full_classification",
                                           clinicalNameNow)]

    # summary stage
    clinicoNow$summary_simple <- "remove"

    "N0" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("N0", "N0 (i-)", "N0 (i+)", "N0 (mol+)")]
    "N1" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("N1", "N1a", "N1b", "N1c", "N1d", "N1mi")]
    "N2" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("N2", "N2a", "N2b", "N2c", "N2d")]
    "N3" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("N3", "N3a", "N3b", "N3c", "N3d")]

    # Collect just tumors low high
    clinicoNow <- clinicoNow[clinicoNow[,"full_classification"] %in% c("low", "high"),]

    # remove non classified lines
    clinicoNow <- clinicoNow[!clinicoNow$summary_simple == "remove",]

    # Drop unused levels
    clinicoNow[,"full_classification"] <-
        droplevels(clinicoNow[,"full_classification"])

    #
    #
    # calculate N sizes
    n_now <- table(clinicoNow[,"full_classification"])
    n_prop_now <- round(prop.table(n_now)*100, 1)


    # check
    if(dim(clinicoNow)[1] > 2){

        # Create table
        contingencyNow <- table(clinicoNow$summary_simple,
                                clinicoNow[,"full_classification"])
        contingencyNow_percent <- round(prop.table(contingencyNow)*100, 1)

        # join information in order
        # as total followed by percents
        ClinicoPathologicalTable_now <- rbind(c(clinicalNameNow, "", "", "", round(chisq.test(contingencyNow)$p.value, 4)),
                                              c("size", n_now[1], n_prop_now[1], n_now[2], n_prop_now[2]),
                                              cbind(rownames(contingencyNow),
                                                    unname(contingencyNow[,1,drop=TRUE]),
                                                    unname(contingencyNow_percent[,1,drop=TRUE]),
                                                    unname(contingencyNow[,2,drop=TRUE]),
                                                    unname(contingencyNow_percent[,2,drop=TRUE]))
        )

        colnames(ClinicoPathologicalTable_now) <- c("type", colnames(contingencyNow)[1], colnames(contingencyNow_percent)[1], colnames(contingencyNow)[2], colnames(contingencyNow_percent)[2])

        # Add or join
        ClinicoPathologicalTable <- rbind(ClinicoPathologicalTable, ClinicoPathologicalTable_now)


    }

    # pathologic_M
    clinicalNameNow <- "pathologic_M"
    # Collect data
    clinicoNow <- framesList[,c("full_classification",
                                           clinicalNameNow)]

    # summary stage
    clinicoNow$summary_simple <- "remove"

    "M0" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("cM0 (i+)", "M0")]
    "M1" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("M1")]

    # Collect just tumors low high
    clinicoNow <- clinicoNow[clinicoNow[,"full_classification"] %in% c("low", "high"),]

    # remove non classified lines
    clinicoNow <- clinicoNow[!clinicoNow$summary_simple == "remove",]

    # Drop unused levels
    clinicoNow[,"full_classification"] <-
        droplevels(clinicoNow[,"full_classification"])


    # calculate N sizes
    n_now <- table(clinicoNow[,"full_classification"])
    n_prop_now <- round(prop.table(n_now)*100, 1)


    # check
    if(dim(clinicoNow)[1] > 2){

        # Create table
        contingencyNow <- table(clinicoNow$summary_simple,
                                clinicoNow[,"full_classification"])
        contingencyNow_percent <- round(prop.table(contingencyNow)*100, 1)

        # join information in order
        # as total followed by percents
        ClinicoPathologicalTable_now <- rbind(c(clinicalNameNow, "", "", "", round(chisq.test(contingencyNow)$p.value, 4)),
                                              c("size", n_now[1], n_prop_now[1], n_now[2], n_prop_now[2]),
                                              cbind(rownames(contingencyNow),
                                                    unname(contingencyNow[,1,drop=TRUE]),
                                                    unname(contingencyNow_percent[,1,drop=TRUE]),
                                                    unname(contingencyNow[,2,drop=TRUE]),
                                                    unname(contingencyNow_percent[,2,drop=TRUE]))
        )

        colnames(ClinicoPathologicalTable_now) <- c("type", colnames(contingencyNow)[1], colnames(contingencyNow_percent)[1], colnames(contingencyNow)[2], colnames(contingencyNow_percent)[2])

        # Add or join
        ClinicoPathologicalTable <- rbind(ClinicoPathologicalTable, ClinicoPathologicalTable_now)
    }

    # gender
    clinicalNameNow <- "gender"
    # Collect data
    clinicoNow <- framesList[,c("full_classification",
                                           clinicalNameNow)]

    # summary stage
    clinicoNow$summary_simple <- "remove"

    "Male" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("MALE")]
    "Female" -> clinicoNow$summary_simple[clinicoNow[,clinicalNameNow] %in% c("FEMALE")]

    # Collect just tumors low high
    clinicoNow <- clinicoNow[clinicoNow[,"full_classification"] %in% c("low", "high"),]

    # remove non classified lines
    clinicoNow <- clinicoNow[!clinicoNow$summary_simple == "remove",]

    # Drop unused levels
    clinicoNow[,"full_classification"] <-
        droplevels(clinicoNow[,"full_classification"])

    #
    #
    # calculate N sizes
    n_now <- table(clinicoNow[,"full_classification"])
    n_prop_now <- round(prop.table(n_now)*100, 1)


    # check
    if(dim(clinicoNow)[1] > 2){

        # Create table
        contingencyNow <- table(clinicoNow$summary_simple,
                                clinicoNow[,"full_classification"])
        contingencyNow_percent <- round(prop.table(contingencyNow)*100, 1)

        # join information in order
        # as total followed by percents
        ClinicoPathologicalTable_now <- rbind(c(clinicalNameNow, "", "", "", round(chisq.test(contingencyNow)$p.value, 4)),
                                              c("size", n_now[1], n_prop_now[1], n_now[2], n_prop_now[2]),
                                              cbind(rownames(contingencyNow),
                                                    unname(contingencyNow[,1,drop=TRUE]),
                                                    unname(contingencyNow_percent[,1,drop=TRUE]),
                                                    unname(contingencyNow[,2,drop=TRUE]),
                                                    unname(contingencyNow_percent[,2,drop=TRUE]))
        )

        colnames(ClinicoPathologicalTable_now) <- c("type", colnames(contingencyNow)[1], colnames(contingencyNow_percent)[1], colnames(contingencyNow)[2], colnames(contingencyNow_percent)[2])

        # Add or join
        ClinicoPathologicalTable <- rbind(ClinicoPathologicalTable, ClinicoPathologicalTable_now)
    }

    # Save clinicopathological table
    write.table(ClinicoPathologicalTable,
                file = file.path(DIR, "kaplan_maier", paste0(Name,"_Clinical_table.txt")), sep="\t")


    clinical_groups <- framesList[, c("full_classification",  "classification")]
    rownames(clinical_groups) <- framesList[, "PatientCode"]

    assign("clinical_groups", clinical_groups, envir = get(envir_link))

    write.csv(clinical_groups, file = file.path(DIR, "clinical_groups"))

    gc()
    message("Done!\n")

}

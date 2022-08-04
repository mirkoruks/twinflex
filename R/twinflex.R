###############################################################################
#' Estimates different types of twin models
#' @rdname twinflex
#' @description twinflex is a wrapper function for OpenMx. It allows you to estimate complex twin models with just one line of code. Furthermore you can combine different types of twin models (univariate, GxE, with or without covariates, Cholesky vs. ACE-beta) individually.
#' @param acevars a string vector containing the variables whose variance/covariance you want to decompose. For multivariate models, just use several variables. At the moment only the Cholesky parametrization is implemented. The first variable in the acevars vector is considered to be the first in this causal or temporary sense as well; the second variable as the second and so on. So have in mind to check the order of your variables in the multivariate case.
#' @param zyg please recode your zygosity variable as follows: 1 = MZ and 2 = DZ, the variable name must be `zyg`
#' @param data a raw data frame (please just remove any SPSS or Stata labels/notes etc. before).
#' @param sep a separator referring to the element that separates the variable names and the twin specific suffix. For example, the variables in the data frame are called `IQ_1` for twin 1 and `IQ_2` for twin 2. The `acevars` argument then is `acevars = 'IQ'` and the separator is `separator = '_'`
#' @param covvars a string vector with covariates. It is possible to use long-formatted (1 variable per twin pair) or wide-formatted (1 variable per twin) variables. Please make sure that you use the wide format only for variables that have a within-twin-pair variance.
#' @param covariance a logical to define whether the covariates should be modelled as random covariates in the model for the covariance matrix or as fixed covariates in the model for the mean.
#' @param type a string vector to define the "type" of a multivariate model. For a Cholesky model (default) use `type = 'Chol'`. For the "ACE-beta" model, use `type = 'aceb'`
#' @param moderation a string vector indicating a moderated path coefficient. Notation is comparable to MPLUS. The effects of covariates cannot be moderated. For a univariate twin model the argument `moderation = Y BY Mod` with `Y` as the outcome and `Mod` as the moderator refers to a moderation of the effects of the ACE components on `Y` by `Mod`. For a multivariate model there are two options which can be combined: a) the argument `moderation = Y BY Mod` refers to a moderation of the effects of the ACE components unique to `Y` on `Y` by `Mod`; b) the argument `moderation = X -> Y BY Mod` refers to a moderation of the effects of the ACE components of `X` on `Y` (and the phenotypic effect of `X` on `Y` when `type = 'aceb'`) by `Mod`.
#' @param ModCov a string vector indicating whether main effect of moderator should be estimated for all variables (`ModCov = "All"`) or for only the variables for which the incoming paths are moderated (`ModCov = "DV"`). Default is `ModCov = "DV"`.
#' @param ordinal a string vector defining the ordinal variables (binary or ordinal)
#' @param TryHard a logical to define whether TryHard should be used (see OpenMx documentation for details). Recommended for complex models.
#' @param Tries a numeric indicating the number of attempts to run the model in addition to the first (see OpenMx documentation for details). Default is 10.
#' @param exh a logical indicating model fitting stops making additional attempts once it reaches an acceptable solution (when `exh = FALSE` as the default).
#' @param Optimizer a string vector indicating the optimizer. Default is `Optimizer = 'SLSQP'` (see OpenMx documentation for details).
#' @param lboundACE a logical indicating whether to set lower bounds for the estimation of the univariate ACE effects.
#' @param dzA a numerical vector indicating the genetic correlation for DZ twins (changeable to adjust for assortative mating). Default is `dzA = 0.5`
#' @param dzC a numerical vector indicating the shared environmental correlation for DZ twins. Default is `dzC = 1`. Change to `dzC = 0.25` to estimate an ADE model.

#' @import OpenMx tidyverse dplyr
#' @importFrom umx umxThresholdMatrix
#' @importFrom stats na.omit pnorm qnorm var
#' @export
# Begin function
twinflex <- function(acevars = NULL, zyg = "zyg", sep = "", data = NULL, covvars = NULL, covariance = FALSE, type = "Chol", moderation = NULL, ModCov = "DV", ordinal = NULL, TryHard = FALSE, Tries = 10, exh = TRUE, Optimizer = "SLSQP", lboundACE = TRUE, dzA = 0.5, dzC = 1) {

    ############################################
    #-------------- Preparations --------------#
    ############################################
    if (is.null(covvars) & covariance == TRUE) {
        stop("If you don't provide any covariates, please set covariance = FALSE")
    }
    ###############################
    # Select acevars
    ###############################
    acevars1 <-    paste0(acevars,sep,"1") # ACE vars twin 1
    acevars2 <-    paste0(acevars,sep,"2") # ACE vars twin 2
    vars_wide <- function(vars, sep, num) {
        result <- paste(vars, rep(num, each = length(vars)), sep = sep)
        return(result)
    }
    acevars_wide <- vars_wide(vars = acevars, sep = sep, num = c(1,2))
    check_acevars <- acevars_wide %in% colnames(data)
    if (FALSE %in% check_acevars) {
        stop(paste0(c("The following variable(s) provided in the 'acevars' argument cannot be found in the data frame provided in the 'data' argument:", paste(acevars2[!check_acevars], collapse = ", ")), collapse=" "))
    }


    ###############################
    # Select covvars
    ###############################
    if (!is.null(covvars)) {
        check_covvars_long <- covvars %in% colnames(data)
        covvars_long <- covvars[check_covvars_long]
        if (FALSE %in% check_covvars_long) {
            covvars_wide <- vars_wide(vars = covvars[!check_covvars_long], sep = sep, num = c(1,2))
            check_covvars_wide <- covvars_wide %in% colnames(data)
            if (FALSE %in% check_covvars_wide) {
                covvars_error_wide <- covvars_wide[!check_covvars_wide][1:(length(covvars_wide[!check_covvars_wide])/2)]
                covvars_error <- gsub(pattern = paste0(sep,"1"), replacement = "",x = covvars_error_wide)
                stop(paste0(c("The following variable(s) provided in the 'covvars' argument cannot be found in the data frame provided in the 'data' argument:", paste(covvars_error, collapse = ", ")), collapse=" "))
            }
        } else {
            covvars_wide <- NULL
        }
        covvars_all <- c(covvars_long, covvars_wide)
        covvars1 <- covvars_wide[1:(length(covvars_wide)/2)]
        covvars2 <- covvars_wide[(1+(length(covvars_wide))/2):length(covvars_wide)]
    } else {
        covvars_all <- NULL
    }
    ###############################
    # select moderators
    ###############################
    # Separate acevars and moderators
    if (!is.null(moderation)) {
        moderatorlist <- list()
        modraw <- unlist(strsplit(moderation, split = "BY"))
        modraw <- gsub(pattern = " ", replacement = "", modraw)
        modraw
        index <- seq(from = 2, to = length(modraw), by = 2)
        start <- 1
        for (i in 1:length(index)) {
            end <- index[i]
            nmod <- 1
            mod_UV <- NA
            modsplit <- modraw[start:end]
            moderator <- modsplit[2]
            if (grepl(pattern = "+", x = moderator)) {
                moderator <- unlist(strsplit(moderator, split = "\\+"))
                nmod <- length(moderator)
            }
            mod_acevars <- modsplit[1]
            mod_AV <- mod_acevars
            if (grepl(pattern = "->", x = mod_acevars)) {
                mod_mult <- TRUE
                mod_AV <- unlist(strsplit(mod_acevars, split = "->"))[2]
                mod_UV <- unlist(strsplit(mod_acevars, split = "->"))[1]
                if (which(acevars == mod_AV) < which(acevars == mod_UV)) {
                    stop("Check your moderations! The order of the variables does not correspond to the order given in the 'acevars' argument!")
                }
                moderatorlist[[i]] <- list(AV = which(acevars == mod_AV), UV = which(acevars == mod_UV), mod_mult = mod_mult, nmod = nmod, moderators = moderator)
            } else {
                mod_mult <- FALSE
                mod_UV <- NA
                moderatorlist[[i]] <- list(AV = which(acevars == mod_AV), UV = NA, mod_mult = mod_mult, nmod = nmod, moderators = moderator)
            }
            start <- end+1
        }

        # Check if acevars exist
        mod_acevars <- NULL
        for (i in 1:length(moderatorlist)) {
            mod_acevars <- c(mod_acevars, acevars[moderatorlist[[i]][["AV"]]])
            if (moderatorlist[[i]][["mod_mult"]] == TRUE) {
                mod_acevars <- c(mod_acevars, acevars[moderatorlist[[i]][["UV"]]])
            }
        }
        mod_acevars <- unique(mod_acevars)
        mod_acevars_wide <- vars_wide(vars = mod_acevars, sep = sep, num = c(1,2))
        check_mod_acevars_wide <- mod_acevars_wide %in% colnames(data)
        if (FALSE %in% check_mod_acevars_wide) {
            stop(paste0(c("The following variable(s) provided in the 'moderation' argument cannot be found in the data frame provided in the 'data' argument:", paste(mod_acevars[!check_mod_acevars_wide], collapse = ", ")), collapse=" "))
        }

        # Check if moderators exist
        mod_modvars <- NULL
        for (i in 1:length(moderatorlist)) {
            mod_modvars <- c(mod_modvars, moderatorlist[[i]][["moderators"]])
        }
        mod_modvars <- unique(mod_modvars)
        check_mod_modvars_long <- mod_modvars %in% colnames(data)
        if (TRUE %in% check_mod_modvars_long) {
            mod_modvars_long <- mod_modvars[check_mod_modvars_long]
        } else {
            mod_modvars_long <- NULL
        }
        if (FALSE %in% check_mod_modvars_long) {
            mod_modvars_wide <- vars_wide(vars = mod_modvars[!check_mod_modvars_long], sep = sep, num = c(1,2))
            check_mod_modvars_wide <- mod_modvars_wide %in% colnames(data)
            if (FALSE %in% check_mod_modvars_wide) {
                mod_modvars_error_wide <- mod_modvars_wide[!check_mod_modvars_wide][1:(length(mod_modvars_wide[!check_mod_modvars_wide])/2)]
                mod_modvars_error <- gsub(pattern = paste0(sep,"1"), replacement = "",x = mod_modvars_error_wide)
                stop(paste0(c("The following variable(s) provided in the 'moderation' argument cannot be found in the data frame provided in the 'data' argument:", paste(mod_modvars_error, collapse = ", ")), collapse=" "))
            }
        } else {
            mod_modvars_wide <- NULL
        }

        moderatorlist2 <- list()
        uniquemod <- NULL
        for (i in 1:length(moderatorlist)) {
            uniquemod <- c(uniquemod, moderatorlist[[i]][["moderators"]])
        }
        uniquemod <- unique(uniquemod)
        if (length(uniquemod) < 5) {
            uniquemod <- c(uniquemod,rep(NA,(5-length(uniquemod))))
        } else if (length(uniquemod) > 5) {
            stop("Max. 5 different moderators are allowed")
        }
        for (i in 1:length(uniquemod)) {
            avvars <- NULL
            uvvars <- NULL
            if (!is.na(uniquemod[i])) {
                for (j in 1:length(moderatorlist)) {
                    if (uniquemod[i] %in% moderatorlist[[j]][["moderators"]]) {
                        avvars <- c(avvars, moderatorlist[[j]][["AV"]])
                        uvvars <- c(uvvars, moderatorlist[[j]][["UV"]])
                    }
                }
            }
            if (uniquemod[i] %in% mod_modvars_long) {
                Long <- TRUE

                Def1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = paste0("data.",uniquemod[i]), name = paste0("defmod",i,"t1"))
                Def2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = paste0("data.",uniquemod[i]), name = paste0("defmod",i,"t2"))
            } else if (paste0(uniquemod[i],sep,"1") %in% mod_modvars_wide) {
                Long <- FALSE
                Def1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = paste0("data.",paste0(uniquemod[i],sep,"1")), name = paste0("defmod",i,"t1"))
                Def2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = paste0("data.",paste0(uniquemod[i],sep,"2")), name = paste0("defmod",i,"t2"))
            } else if (is.na(uniquemod[i])) {
                avvars <- NA
                uvvars <- NA
                Long <- NA
                Def1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = NA, name = paste0("defmod",i,"t1"))
                Def2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = NA, name = paste0("defmod",i,"t2"))
            }
            moderatorlist2[[i]] <- list(AV = avvars, UV = uvvars, Long = Long, Def1 = Def1, Def2 = Def2)
        }
    } else {
        mod_modvars_long <- NULL
        mod_modvars_wide <- NULL
    }

    #############################################
    #-------- Store number of variables --------#
    #############################################
    nv <- length(acevars_wide)/2
    ntv <- nv*2

    nACE <- nv*3
    ntACE <- ntv*3

    if (!is.null(covvars)) {
        ncl <- length(covvars_long)
        ncw <- length(covvars_wide)/2
    } else {
        ncl <- 0
        ncw <- 0
    }
    ntcl <- ncl
    ntcw <- ncw*2
    ntc <- ntcl + ntcw


    if (covariance == TRUE) {
        ntotal <- ntv + ntc + ntACE
    } else {
        ntotal <- ntv + ntACE
    }

    #####################################################
    #-------- Check if binary variables present --------#
    #####################################################
    binary <- NULL
    if (!is.null(ordinal)) {
        for (i in 1:length(ordinal)) {
            ordvar <- paste0(ordinal[i],sep,"1")
            ordval <- unique(data[,ordvar][!is.na(data[,ordvar])])
            if (length(ordval) == 2) {
                binary <- c(binary, TRUE)
            } else {
                binary <- NULL
            }
        }
    }
    if (!is.null(binary) & covariance == TRUE) {
        stop("At the moment, the combination 'random covariates' and 'binary outcomes' is still quite experimental.")
    }
    ############################################################
    # Matrix A
    ############################################################

    ##############################
    # Lower and Gap
    ##############################
    if (covariance == FALSE) {
        AMatLower <- mxMatrix(type = "Zero",
                              nrow = ntACE,
                              ncol = ntotal,
                              name = "ALower")
    } else {
        AMatLower <- mxMatrix(type = "Zero",
                              nrow = ntACE + ntc,
                              ncol = ntotal,
                              name = "ALower")
    }
    AMatGap <- mxMatrix(type = "Zero",
                        nrow = nv,
                        ncol = nv,
                        name = "AGap")

    ##############################
    # ACE
    ##############################
    LabA <- matrix(paste0("a", 1:nv, rep(1:nv, each = nv)), nrow = nv, ncol = nv)
    LabA[upper.tri(LabA, diag = FALSE)] <- NA
    FreeA <- !is.na(LabA)

    if (dzC == 0.25) {
        LabC <- matrix(paste0("d", 1:nv, rep(1:nv, each = nv)), nrow = nv, ncol = nv)
    } else {
        LabC <- matrix(paste0("c", 1:nv, rep(1:nv, each = nv)), nrow = nv, ncol = nv)
    }
    LabC[upper.tri(LabC, diag = FALSE)] <- NA
    FreeC <- !is.na(LabC)
    LabE <- matrix(paste0("e", 1:nv, rep(1:nv, each = nv)), nrow = nv, ncol = nv)

    varace <- sqrt(diag(as.matrix(var(data[,acevars1], use = "pairwise"))))
    ValACE <- matrix(rep(varace, each = nv), nrow = nv, ncol = nv, byrow = TRUE)
    diag(ValACE) <- diag(ValACE)*0.3
    ValACE[lower.tri(ValACE, diag = FALSE)] <- ValACE[lower.tri(ValACE, diag = FALSE)]*0.2
    ValACE[upper.tri(ValACE, diag = FALSE)] <- 0

    ValA <- ValACE
    ValC <- ValACE

    LboundACE <- matrix(NA, nrow = nv, ncol = nv)
    if (lboundACE == TRUE) {
        diag(LboundACE) <- 0.001
    }

    if (type == "aceb") {
        LabE[upper.tri(LabE, diag = FALSE)] <- NA
        LabE[lower.tri(LabE, diag = FALSE)] <- NA
        FreeE <- !is.na(LabE)
        ValE <- ValACE
        ValE[lower.tri(ValE, diag = FALSE)] <- 0
    } else {
        LabE[upper.tri(LabE, diag = FALSE)] <- NA
        FreeE <- !is.na(LabE)
        ValE <- ValACE
    }


    matA <- mxMatrix(type = "Full",
                     nrow = nv,
                     ncol = nv,
                     labels = LabA,
                     values = ValA,
                     lbound = LboundACE,
                     free = FreeA,
                     name = "AA")

    matC <- mxMatrix(type = "Full",
                     nrow = nv,
                     ncol = nv,
                     labels = LabC,
                     values = ValC,
                     lbound = LboundACE,
                     free = FreeC,
                     name = "AC")

    matE <- mxMatrix(type = "Full",
                     nrow = nv,
                     ncol = nv,
                     labels = LabE,
                     values = ValE,
                     lbound = LboundACE,
                     free = FreeE,
                     name = "AE")

    if (!is.null(moderation)) {
        ACEmodLabList <- list()
        for (i in 1:length(moderatorlist2)) {
            AmodLab <- matrix(NA, nrow = nv, ncol = nv)
            CmodLab <- matrix(NA, nrow = nv, ncol = nv)
            EmodLab <- matrix(NA, nrow = nv, ncol = nv)
            BmodLab <- matrix(NA, nrow = nv, ncol = nv)
            for (j in 1:length(moderatorlist2[[i]][["AV"]])) {
                AV <- moderatorlist2[[i]][["AV"]][j]

                if (!is.na(moderatorlist2[[i]][["UV"]][j])) {
                    UV <- moderatorlist2[[i]][["UV"]][j]
                } else {
                    UV <- AV
                }
                AmodLab[AV,UV] <- paste0("bm",i,"a",AV,UV)
                CmodLab[AV,UV] <- paste0("bm",i,"c",AV,UV)
                EmodLab[AV,UV] <- paste0("bm",i,"e",AV,UV)
                BmodLab[AV,UV] <- paste0("bm",i,"b",AV,UV)
            }
            if (type == "aceb") {
                EmodLab[lower.tri(EmodLab, diag = FALSE)] <- NA
                diag(BmodLab) <- NA
            } else {
                BmodLab <- NA
            }
            AFree <- !is.na(AmodLab)
            Amod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, free = AFree, labels = AmodLab, name = paste0("AAM",i))
            CFree <- !is.na(CmodLab)
            Cmod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, free = CFree, labels = CmodLab, name = paste0("ACM",i))
            EFree <- !is.na(EmodLab)
            Emod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, free = EFree, labels = EmodLab, name = paste0("AEM",i))
            BFree <- !is.na(BmodLab)
            Bmod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, free = BFree, labels = BmodLab, name = paste0("ABM",i))
            ACEmodLabList[[i]] <- list(A = Amod, C = Cmod, E = Emod, B = Bmod)
        }

        # Save matrices with interaction effects as R objects
        for (i in 1:length(ACEmodLabList)) {
            assign(paste0("AMod",i), ACEmodLabList[[i]][["A"]])
            assign(paste0("CMod",i), ACEmodLabList[[i]][["C"]])
            assign(paste0("EMod",i), ACEmodLabList[[i]][["E"]])
            assign(paste0("BMod",i), ACEmodLabList[[i]][["B"]])
        }
        # Save definition variables with moderators as R objects
        for (i in 1:length(moderatorlist2)) {
            assign(paste0("Mod",i,"Def1"), moderatorlist2[[i]][["Def1"]])
            assign(paste0("Mod",i,"Def2"), moderatorlist2[[i]][["Def2"]])
        }
        aModFull1 <- mxAlgebra(expression = AA + AAM1*defmod1t1 + AAM2*defmod2t1 + AAM3*defmod3t1 + AAM4*defmod4t1 + AAM5*defmod5t1, name = "AMod1")
        aModFull2 <- mxAlgebra(expression = AA + AAM1*defmod1t2 + AAM2*defmod2t2 + AAM3*defmod3t2 + AAM4*defmod4t2 + AAM5*defmod5t2, name = "AMod2")

        cModFull1 <- mxAlgebra(expression = AC + ACM1*defmod1t1 + ACM2*defmod2t1 + ACM3*defmod3t1 + ACM4*defmod4t1 + ACM5*defmod5t1, name = "CMod1")
        cModFull2 <- mxAlgebra(expression = AC + ACM1*defmod1t2 + ACM2*defmod2t2 + ACM3*defmod3t2 + ACM4*defmod4t2 + ACM5*defmod5t2, name = "CMod2")

        eModFull1 <- mxAlgebra(expression = AE + AEM1*defmod1t1 + AEM2*defmod2t1 + AEM3*defmod3t1 + AEM4*defmod4t1 + AEM5*defmod5t1, name = "EMod1")
        eModFull2 <- mxAlgebra(expression = AE + AEM1*defmod1t2 + AEM2*defmod2t2 + AEM3*defmod3t2 + AEM4*defmod4t2 + AEM5*defmod5t2, name = "EMod2")

        bModFull1 <- mxAlgebra(expression = AB + ABM1*defmod1t1 + ABM2*defmod2t1 + ABM3*defmod3t1 + ABM4*defmod4t1 + ABM5*defmod5t1, name = "BMod1")
        bModFull2 <- mxAlgebra(expression = AB + ABM1*defmod1t2 + ABM2*defmod2t2 + ABM3*defmod3t2 + ABM4*defmod4t2 + ABM5*defmod5t2, name = "BMod2")


        matACEFull <- mxAlgebra(expression = rbind(cbind(AMod1, CMod1, EMod1, AGap, AGap, AGap),
                                                   cbind(AGap, AGap, AGap, AMod2, CMod2, EMod2)),
                                name = "ACEF")
        matACEFullNull <- mxAlgebra(expression = rbind(cbind(AA, AC, AE, AGap, AGap, AGap),
                                                       cbind(AGap, AGap, AGap, AA, AC, AE)),
                                    name = "ACEFNull")
    } else {
        matACEFull <- mxAlgebra(expression = rbind(cbind(AA, AC, AE, AGap, AGap, AGap),
                                                   cbind(AGap, AGap, AGap, AA, AC, AE)),
                                name = "ACEF")
        matACEFullNull <- NULL
    }

    ##############################
    # Beta
    ##############################
    LabB <- matrix(paste0("b", 1:nv, rep(1:nv, each = nv)), nrow = nv, ncol = nv)
    LabB[upper.tri(LabB, diag = TRUE)] <- NA
    FreeB <- !is.na(LabB)
    varace <- var(data[,acevars1], use = "pairwise")
    nenner <- diag(varace)^-1
    nenner <- (diag(varace)^-1)
    nenner <- matrix(nenner, nrow = nv, ncol = nv, byrow = TRUE)
    nenner[upper.tri(nenner, diag = TRUE)] <- 0
    zaehler <- varace
    zaehler[upper.tri(zaehler, diag = TRUE)] <- 0
    ValB <- nenner*zaehler
    rownames(ValB) <- NULL
    colnames(ValB) <- NULL
    if (type != "aceb") {
        ValB <- 0
        FreeB <- FALSE
    }
    matB <- mxMatrix(type = "Full",
                     nrow = nv,
                     ncol = nv,
                     labels = LabB,
                     values = ValB,
                     free = FreeB,
                     name = "AB")
    if (!is.null(moderation)) {
        matBFull <- mxAlgebra(expression = rbind(cbind(BMod1,AGap),
                                                 cbind(AGap,BMod2)),
                              name = "ABF")
        matBFullNull <- mxAlgebra(expression = rbind(cbind(AB,AGap),
                                                     cbind(AGap,AB)),
                                  name = "ABFNull")
    } else {
        matBFull <- mxAlgebra(expression = rbind(cbind(AB,AGap),
                                                 cbind(AGap,AB)),
                              name = "ABF")
    }

    ##############################
    # Covariates
    ##############################
    if (!is.null(covvars) & covariance == TRUE) {
        if (ncl > 0) {
            LabCovLong <- paste0("bcovl", 1:nv, rep(1:ncl, each = nv))
            Ycov1 <- as.matrix(subset(na.omit(data), select = acevars1))
            Ycov2 <- as.matrix(subset(na.omit(data), select = acevars2))
            XcovLong <- cbind(1,as.matrix(subset(na.omit(data), select = c(covvars_long))))
            pathCovLongStart1 <- (t(solve(t(XcovLong)%*%XcovLong)%*%t(XcovLong)%*%Ycov1))[,2:ncol(XcovLong)]
            pathCovLongStart2 <- (t(solve(t(XcovLong)%*%XcovLong)%*%t(XcovLong)%*%Ycov2))[,2:ncol(XcovLong)]
            ValCovLong <- (pathCovLongStart1 + pathCovLongStart2) / 2
            colnames(ValCovLong) <- NULL
            rownames(ValCovLong) <- NULL
        }
        if (ncw > 0) {
            # Wide
            LabCovWide <- paste0("bcovw", 1:nv, rep(1:ncw, each = nv))
            Ycov1 <- as.matrix(subset(na.omit(data), select = acevars1))
            Ycov2 <- as.matrix(subset(na.omit(data), select = acevars2))
            XCovWide1 <- cbind(1,as.matrix(subset(na.omit(data), select = c(covvars1))))
            XCovWide2 <- cbind(1,as.matrix(subset(na.omit(data), select = c(covvars2))))
            pathCovWideStart11 <- (t(solve(t(XCovWide1)%*%XCovWide1)%*%t(XCovWide1)%*%Ycov1))[,2:ncol(XCovWide1)]
            pathCovWideStart21 <- (t(solve(t(XCovWide1)%*%XCovWide1)%*%t(XCovWide1)%*%Ycov2))[,2:ncol(XCovWide1)]
            pathCovWideStart12 <- (t(solve(t(XCovWide2)%*%XCovWide2)%*%t(XCovWide2)%*%Ycov1))[,2:ncol(XCovWide2)]
            pathCovWideStart22 <- (t(solve(t(XCovWide2)%*%XCovWide2)%*%t(XCovWide2)%*%Ycov2))[,2:ncol(XCovWide2)]
            ValCovWide <- (pathCovWideStart11 + pathCovWideStart21 + pathCovWideStart12 + pathCovWideStart22) / 4
            colnames(ValCovWide) <- NULL
            rownames(ValCovWide) <- NULL
        }
        if (ncl > 0) {
            AMatCovLong <- mxMatrix(type = "Full",
                                    nrow = nv,
                                    ncol = ncl,
                                    labels = LabCovLong,
                                    values = ValCovLong,
                                    free = TRUE,
                                    name = "ACL")
        }
        if (ncw > 0) {
            AMatCovWide <- mxMatrix(type = "Full",
                                    nrow = nv,
                                    ncol = ncw,
                                    labels = LabCovWide,
                                    values = ValCovWide,
                                    free = TRUE,
                                    name = "ACW")
            AMatCovGap <- mxMatrix(type = "Zero",
                                   nrow = nv,
                                   ncol = ncw,
                                   name = "ACGap")
            if (ncl > 0) {
                AMatCov <- mxAlgebra(expression = rbind(cbind(ACL,ACW,ACGap),
                                                        cbind(ACL,ACGap,ACW)),
                                     name = "ACov")
            } else {
                AMatCov <- mxAlgebra(expression = rbind(cbind(ACW,ACGap),
                                                        cbind(ACGap,ACW)),
                                     name = "ACov")
            }
        } else {
            AMatCov <- mxAlgebra(expression = rbind(ACL,ACL),
                                 name = "ACov")
        }
        AMat <- mxAlgebra(expression = rbind(cbind(ABF, ACov, ACEF),
                                             ALower),
                          name = "A")
        if (!is.null(moderation) & !is.null(binary)) {
            AMatNull <- mxAlgebra(expression = rbind(cbind(ABFNull, ACov, ACEFNull),
                                                     ALower),
                                  name = "ANull")
        }

    } else {
        AMat <- mxAlgebra(expression = rbind(cbind(ABF, ACEF),
                                             ALower),
                          name = "A")


        if (!is.null(moderation) & !is.null(binary)) {
            AMatNull <- mxAlgebra(expression = rbind(cbind(ABFNull, ACEFNull),
                                                     ALower),
                                  name = "ANull")
        }
    }


    ##############################
    # Matrix S
    ##############################

    SMatUpper <- mxMatrix(type = "Zero",
                          nrow = ntv,
                          ncol = ntotal,
                          name = "SUpper")

    SMatACE <- mxMatrix(type = "Diag",
                        nrow = nv,
                        ncol = nv,
                        values = 1,
                        name = "SACE")

    if (covariance == FALSE) {
        SMatLower <- mxMatrix(type = "Zero",
                              nrow = ntACE,
                              ncol = ntv,
                              name = "SLower")
    } else {
        SMatLower <- mxMatrix(type = "Zero",
                              nrow = ntACE,
                              ncol = ntv + ntc,
                              name = "SLower")
    }

    SMatGap <- mxMatrix(type = "Zero",
                        nrow = nv,
                        ncol = nv,
                        name = "SGap")

    SMat1MZ <- mxAlgebra(expression = rbind(cbind(SACE,SGap,SGap,SACE,SGap,SGap),
                                            cbind(SGap,SACE,SGap,SGap,SACE,SGap),
                                            cbind(SGap,SGap,SACE,SGap,SGap,SGap),
                                            cbind(SACE,SGap,SGap,SACE,SGap,SGap),
                                            cbind(SGap,SACE,SGap,SGap,SACE,SGap),
                                            cbind(SGap,SGap,SGap,SGap,SGap,SACE)),
                         name = "SACEMZ")

    matdzA <- mxMatrix(type = "Full", nrow =  1, ncol = 1, free = FALSE, values = dzA, name = "dzA")
    matdzC <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = dzC, name = "dzC")

    SMat1DZ <- mxAlgebra(expression = rbind(cbind(SACE,SGap,SGap,dzA*SACE,SGap,SGap),
                                            cbind(SGap,SACE,SGap,SGap,dzC*SACE,SGap),
                                            cbind(SGap,SGap,SACE,SGap,SGap,SGap),
                                            cbind(dzA*SACE,SGap,SGap,SACE,SGap,SGap),
                                            cbind(SGap,dzC*SACE,SGap,SGap,SACE,SGap),
                                            cbind(SGap,SGap,SGap,SGap,SGap,SACE)),
                         name = "SACEDZ")

    if (!is.null(covvars) & covariance == TRUE) {
        # Labels
        if (ncw > 0) {
            CovSLabCWW <- matrix(paste0("covwcw", 1:ncw, "cw", rep(1:ncw, each = ncw)), nrow = ncw, ncol = ncw)
            CovSLabCWW[upper.tri(CovSLabCWW, diag = FALSE)] <- (CovSLabCWW)[lower.tri(CovSLabCWW, diag = FALSE)]
            diag(CovSLabCWW) <- paste0("varcw",1:ncw)
            CovSLabCWB <- matrix(paste0("covbcw", 1:ncw, "cw", rep(1:ncw, each = ncw)), nrow = ncw, ncol = ncw)
            CovSLabCWB[upper.tri(CovSLabCWB, diag = FALSE)] <- (CovSLabCWB)[lower.tri(CovSLabCWB, diag = FALSE)]
            CovSLabCW <- rbind(cbind(CovSLabCWW, CovSLabCWB),
                               cbind(CovSLabCWB, CovSLabCWW))
        }
        if (ncl > 0) {
            CovSLabCL <- matrix(paste0("covcl", 1:ncl, "cl", rep(1:ncl, each = ncl)), nrow = ncl, ncol = ncl)
            CovSLabCL[upper.tri(CovSLabCL, diag = FALSE)] <- (CovSLabCL)[lower.tri(CovSLabCL, diag = FALSE)]
            diag(CovSLabCL) <- paste0("varcl",1:ncl)
            if (ncw > 0) {
                CovSLabCLW <- matrix(paste0(rep(paste0("covcl",1:ncl), each = ntcw), rep(paste0("cw",1:ncw), 2)), nrow = ntcw, ncol = ntcl)
                CovSLabCov <- rbind(cbind(CovSLabCL, t(CovSLabCLW)),
                                    cbind(CovSLabCLW, CovSLabCW))
            } else {
                CovSLabCov <- CovSLabCL
            }
        } else {
            CovSLabCov <- CovSLabCW
        }
        # Start values
        if (ncw > 0) {
            CovSValCWW <- suppressWarnings(var(data[,covvars1], use = "pairwise"))
            CovSValCWB <- suppressWarnings(var(data[,covvars_wide], use = "pairwise"))[(ncw+1):ntcw,1:ncw]
            CovSValCWB[upper.tri(CovSValCWB, diag = FALSE)] <- (CovSValCWB)[lower.tri(CovSValCWB, diag = FALSE)]
            CovSValCW <- rbind(cbind(CovSValCWW, CovSValCWB),
                               cbind(CovSValCWB, CovSValCWW))
        }
        if (ncl > 0) {
            CovSValCL <- suppressWarnings(var(data[,covvars_long], use = "pairwise"))
            if (ncw > 0) {
                CovSValCLW <- matrix(rep(suppressWarnings(var(data[,c(covvars_long,covvars1)], use = "pairwise"))[(ncl+1):(ncl+ncw),1:ncl],2), byrow = TRUE, nrow = ntcw, ncol = ntcl)
                CovSValCov <- rbind(cbind(CovSValCL, t(CovSValCLW)),
                                    cbind(CovSValCLW, CovSValCW))
                colnames(CovSValCov) <- NULL
                rownames(CovSValCov) <- NULL
            } else {
                CovSValCov <- CovSValCL
            }
        } else {
            CovSValCov <- CovSValCWW
        }

        if (any(is.na(diag(CovSValCov)))) {
            if (any(is.na(CovSValCov[lower.tri(CovSValCov, diag = FALSE)]))) {
                stop("I tried to find some good starting values for the covariance matrix of the covariates. However, for some of your covariates no variances and covariances can be calculated based on the data frame provided in 'data'. This probably results in problems estimating the model parameters. So please check your covariates!")
            } else {
                stop("I tried to find some good starting values for the covariance matrix of the covariates. However, for some of your covariates no variances can be calculated based on the data frame provided in 'data'. This probably results in problems estimating the model parameters. So please check your covariates!")
            }
        }
        CovSLboCov <- matrix(NA, nrow = ntc, ncol = ntc)
        diag(CovSLboCov) <- 0.0001
        SMatCov <- mxMatrix(type = "Full",
                            nrow = ntc,
                            ncol = ntc,
                            labels = CovSLabCov,
                            free = TRUE,
                            values = CovSValCov,
                            lbound = CovSLboCov,
                            name = "SCov")

        SMatCovGap1 <- mxMatrix(type = "Zero",
                                nrow = ntc,
                                ncol = ntv,
                                name = "SCovG1")
        SMatCovGap2 <- mxMatrix(type = "Zero",
                                nrow = ntc,
                                ncol = ntACE,
                                name = "SCovG2")

        SMatMZ <- mxAlgebra(expression = rbind(SUpper,
                                               cbind(SCovG1, SCov, SCovG2),
                                               cbind(SLower,SACEMZ)),
                            name = "SMZ")
        SMatDZ <- mxAlgebra(expression = rbind(SUpper,
                                               cbind(SCovG1, SCov, SCovG2),
                                               cbind(SLower,SACEDZ)),
                            name = "SDZ")
    } else {
        SMatMZ <- mxAlgebra(expression = rbind(SUpper,
                                               cbind(SLower,SACEMZ)),
                            name = "SMZ")
        SMatDZ <- mxAlgebra(expression = rbind(SUpper,
                                               cbind(SLower,SACEDZ)),
                            name = "SDZ")
    }

    #########
    ## Fix Variance of binary var to 1
    #########
    if (!is.null(ordinal) & !is.null(binary)) {
        nbin <- sum(binary)
        nord <- length(binary)-sum(binary)
        filterma <- matrix(as.numeric(c(acevars %in% ordinal[binary],rep(FALSE,nv))), byrow = TRUE, nrow = nbin, ncol = ntv)
        binindex <- which(acevars %in% ordinal[binary])
        for (i in 1:nrow(filterma)) {
            filterma[i,-binindex[i]] <- 0
        }
        if (covariance == TRUE) {
            filterma <- cbind(filterma, matrix(0, nrow = nbin, ncol = ntc))
            filterbin <- mxMatrix(type = "Full", nrow = nbin, ncol = ntv + ntc,
                                  values = filterma,
                                  name = "fbin")
        } else {
            filterbin <- mxMatrix(type = "Full", nrow = nbin, ncol = ntv,
                                  values = filterma,
                                  name = "fbin")
        }

        SMatACEVar1 <- mxAlgebra(expression = rbind(cbind(SACE,SGap,SGap,SACE,SGap,SGap),
                                                    cbind(SGap,SACE,SGap,SGap,SACE,SGap),
                                                    cbind(SGap,SGap,SACE,SGap,SGap,SGap),
                                                    cbind(SACE,SGap,SGap,SACE,SGap,SGap),
                                                    cbind(SGap,SACE,SGap,SGap,SACE,SGap),
                                                    cbind(SGap,SGap,SGap,SGap,SGap,SACE)),
                                 name = "SACEVar1")
        if (covariance == FALSE) {
            SMatVar1 <- mxAlgebra(expression = rbind(SUpper,
                                                     cbind(SLower,SACEVar1)),
                                  name = "SVar1")
        } else {
            SMatVar1 <- mxAlgebra(expression = rbind(SUpper,
                                                     cbind(SCovG1, SCov, SCovG2),
                                                     cbind(SLower,SACEVar1)),
                                  name = "SVar1")
        }

        if (!is.null(moderation)) {
            CovConstr <- mxAlgebra(expression = Filter%*%solve(I-ANull)%*%SVar1%*%t(solve(I-ANull))%*%t(Filter), name = "expCovVar1")
            binarycov <- mxAlgebra(expression = fbin %*%expCovVar1 %*% t(fbin), name = "binCov")
        } else {
            CovConstr <- mxAlgebra(expression = Filter%*%solve(I-A)%*%SVar1%*%t(solve(I-A))%*%t(Filter), name = "expCovVar1")
            binarycov <- mxAlgebra(expression = fbin %*%expCovVar1 %*% t(fbin), name = "binCov")
        }
        one <- mxMatrix(type = "Unit", nrow = nbin, ncol = 1, name = "Unit")
        var1 <- mxConstraint(expression = diag2vec(binCov)==Unit , name = "VConstraint1")
        parsvarbin <- c(SMatACEVar1, SMatVar1, CovConstr, binarycov, one, var1, SMatACE, SMatGap, filterbin)
    }

    #########
    ## MEAN
    #########

    # Vector of latent variables (all set to zero) just need to concatenate them to the manifests to get the dimensions right
    latentmeans <- mxMatrix(type = "Full", nrow = 1, ncol = ntACE, name = "lmeans")


    if (!is.null(moderation)) {
        pathModMLabVec <- NULL
        for (i in 1:length(moderatorlist2)) {
            pathModMLabVec <- c(pathModMLabVec, paste0("bModM", 1:nv, rep(i, each = nv)), rep(NA,ntv), paste0("bModM", 1:nv, rep(i, each = nv)))
        }
        pathModMLab <- matrix(pathModMLabVec, byrow = TRUE, nrow = 10, ncol = ntv)

        pathModMFree <- matrix(FALSE, nrow = 10, ncol = ntv)
        for (i in 1:length(moderatorlist2)) {
            if (!(uniquemod[i] %in% acevars)) {
                if (ModCov == "DV") {
                pathModMFree[((i*2)-1),unique(moderatorlist2[[i]][["AV"]])] <- TRUE
                pathModMFree[((i*2)),unique(moderatorlist2[[i]][["AV"]])+nv] <- TRUE
                }
                if (ModCov == "All" & i <= length(moderator)) {
                pathModMFree[((i*2)-1),1:nv] <- TRUE
                pathModMFree[((i*2)),((1:nv)+nv)] <- TRUE
                }
            }
        }

        pathModM <- mxMatrix(type = "Full", nrow = 10, ncol = ntv, byrow = FALSE,
                             free = pathModMFree,
                             values = 0,
                             labels = pathModMLab,
                             name = "pModM")

        DefModMLab <- NULL
        for (i in 1:length(moderatorlist2)) {
            DefModMLab <- c(DefModMLab,moderatorlist2[[i]][["Def1"]]@labels,moderatorlist2[[i]][["Def2"]]@labels)
        }
        DefModM      <- mxMatrix( type="Full", nrow=1, ncol=10, free=FALSE, labels=DefModMLab, name="dModM" )
        RegModM <- mxAlgebra(expression = dModM %*% pModM, name = "RegModM")
    }

    if (covariance == FALSE & !is.null(covvars)) {
        if (ncl > 0) {
            LabCovMLong <- paste0("bcov", 1:nv, rep(1:(ncl), each = ntv))
            if (ncw > 0) {
                LabCovMWide <- c(paste0("bcov", 1:nv, rep((ncl+1):(ncw+ncl), each = ntv)),paste0("bcov", 1:nv, rep((ncl+1):(ncw+ncl), each = ntv)))
            } else {
                LabCovMWide <- NULL
            }
        } else {
            LabCovMLong <- NULL
            LabCovMWide <- c(paste0("bcov", 1:nv, rep((1):(ncw), each = ntv)),paste0("bcov", 1:nv, rep((1):(ncw), each = ntv)))
        }
        LabCovM <- matrix(c(LabCovMLong, LabCovMWide), byrow = TRUE, nrow = ntc, ncol = ntv)


        # Long
        if (ncl > 0) {
            Ycov1 <- as.matrix(subset(na.omit(data), select = c(acevars1,acevars1)))
            Ycov2 <- as.matrix(subset(na.omit(data), select = c(acevars2,acevars2)))
            XcovLong <- cbind(1,as.matrix(subset(na.omit(data), select = c(covvars_long))))
            pathCovLongStart1 <- (t(solve(t(XcovLong)%*%XcovLong)%*%t(XcovLong)%*%Ycov1))[,2:ncol(XcovLong)]
            pathCovLongStart2 <- (t(solve(t(XcovLong)%*%XcovLong)%*%t(XcovLong)%*%Ycov2))[,2:ncol(XcovLong)]
            ValCovLong <- t((pathCovLongStart1 + pathCovLongStart2) / 2)
            colnames(ValCovLong) <- NULL
            rownames(ValCovLong) <- NULL
        } else {
            ValCovLong <- NULL
        }
        if (ncw > 0) {
            # Wide
            Ycov1 <- as.matrix(subset(na.omit(data), select = c(acevars1,acevars1)))
            Ycov2 <- as.matrix(subset(na.omit(data), select = c(acevars2,acevars2)))
            XCovWide1 <- cbind(1,as.matrix(subset(na.omit(data), select = c(covvars1))))
            XCovWide2 <- cbind(1,as.matrix(subset(na.omit(data), select = c(covvars2))))
            pathCovWideStart11 <- (t(solve(t(XCovWide1)%*%XCovWide1)%*%t(XCovWide1)%*%Ycov1))[,2:ncol(XCovWide1)]
            pathCovWideStart21 <- (t(solve(t(XCovWide1)%*%XCovWide1)%*%t(XCovWide1)%*%Ycov2))[,2:ncol(XCovWide1)]
            pathCovWideStart12 <- (t(solve(t(XCovWide2)%*%XCovWide2)%*%t(XCovWide2)%*%Ycov1))[,2:ncol(XCovWide2)]
            pathCovWideStart22 <- (t(solve(t(XCovWide2)%*%XCovWide2)%*%t(XCovWide2)%*%Ycov2))[,2:ncol(XCovWide2)]
            ValCovWide <- t((pathCovWideStart11 + pathCovWideStart21 + pathCovWideStart12 + pathCovWideStart22) / 4)
            colnames(ValCovWide) <- NULL
            rownames(ValCovWide) <- NULL
        } else {
            ValCovWide <- NULL
        }
        ValCovM <- rbind(ValCovLong, ValCovWide,ValCovWide)


        pathCovM <- mxMatrix(type = "Full", nrow = ntc, ncol = ntv, byrow = FALSE,
                             free = TRUE,
                             values = ValCovM,
                             labels = LabCovM,
                             name = "pCovM") #Dimension = c*ntv

        # Matrix of definition variables for main effect of covariates in mean
        DefCovMLab <- paste0("data",".",covvars_all)
        DefCovM      <- mxMatrix( type="Full", nrow=1, ncol=ntc, free=FALSE, labels=DefCovMLab, name="dMCov" )
        # Matrix effects on means


        if (!is.null(moderation)) {
            effCovM <- mxAlgebra(expression = dMCov%*%pCovM + RegModM, name = "eCovM")
        } else {
            effCovM <- mxAlgebra(expression = dMCov%*%pCovM, name = "eCovM")
        }
        effMeanCovFull <- mxAlgebra(expression = cbind(eCovM,lmeans), name = "effMCovFull")
    } else if (!is.null(moderation)) {
        effCovM <- mxAlgebra(expression = RegModM, name = "eCovM")
        if (covariance == TRUE) {
            covmeans <- mxMatrix(type = "Full", nrow = 1, ncol = ntc, name = "cmeans")
            effMeanCovFull <- mxAlgebra(expression = cbind(eCovM,cmeans, lmeans), name = "effMCovFull")
        } else {
            effMeanCovFull <- mxAlgebra(expression = cbind(eCovM, lmeans), name = "effMCovFull")
        }
    }

    # Unmoderated means
    if (length(acevars) == 1) {
        MeanVal1 <- c(mean(data[,acevars1], na.rm = TRUE))
    } else {
        MeanVal1 <- c(colMeans(data[,acevars1], na.rm = TRUE))
    }

    MeanLab1 <- c(rep(paste0("mu",1:nv), 2))

    if (ncl > 0 & covariance == TRUE) {
        MeanLabCovL <- c(paste0("mucl",1:ncl))

        if (ncl == 1) {
            MeanValCovL <- c(mean(data[,covvars_long], na.rm = TRUE))
        } else {
            MeanValCovL <- c(colMeans(data[,covvars_long], na.rm = TRUE))
        }
    } else {
        MeanLabCovL <- NULL
        MeanValCovL <- NULL
    }
    if (ncw > 0 & covariance == TRUE) {
        MeanLabCovW <- c(rep(paste0("mucw",1:ncw), 2))
        if (length(covvars1 ) == 1) {
            MeanValCovW <- rep((mean(data[,covvars_wide[1:ncw]], na.rm = TRUE) + mean(data[,covvars_wide[(ncw+1):ntcw]], na.rm = TRUE))/2,2)
        } else {
            MeanValCovW <- rep((colMeans(data[,covvars_wide[1:ncw]], na.rm = TRUE) + colMeans(data[,covvars_wide[(ncw+1):ntcw]], na.rm = TRUE))/2,2)
        }
    } else {
        MeanLabCovW <- NULL
        MeanValCovW <- NULL
    }
    MeanLabCov <- c(MeanLabCovL, MeanLabCovW)
    MeanValCov <- c(MeanValCovL, MeanValCovW)
    MeanLabACE <- rep(c(paste0("Amu",1:nv), paste0("Cmu",1:nv), paste0("Emu",1:nv)), 2)
    MeanValACE <- rep(0, ntACE)
    MeanLab <- c(MeanLab1, MeanLabCov, MeanLabACE)
    MeanFree1 <- rep(TRUE, nv)
    if (!is.null(ordinal) & !is.null(binary)) {
        #which(acevars_wide %in% paste0(ordinal[binary],sep,c(1:2)))

        MeanFree1[binindex] <- FALSE
        MeanVal1[binindex] <- 0
    }
    if (!is.null(covvars) & covariance == TRUE) {
        MeanCovFree <- rep(TRUE, ntc)
    } else {
        MeanCovFree <- NULL
    }
    MeanFree <- c(MeanFree1, MeanFree1, MeanCovFree, rep(FALSE, ntACE))
    MeanVal <- c(MeanVal1, MeanVal1, MeanValCov, MeanValACE)

    M0 <- mxMatrix(type = "Full", nrow = ntotal, ncol = 1, free = MeanFree, values = MeanVal, labels = MeanLab, name = "M")

    #########
    ## Threshold Matrix
    #########
    if (!is.null(ordinal)) {
        ordinal_wide <- vars_wide(ordinal, sep, num = c(1,2))
        for (i in 1:length(ordinal_wide)) {
            categories <- sort(unique(data[,ordinal_wide[i]][!is.na(data[,ordinal_wide[i]])]))
            data[,ordinal_wide[i]] <- mxFactor(data[,ordinal_wide[i]], levels = categories)
        }
        thresh <- umxThresholdMatrix(df = data, fullVarNames = ordinal_wide, method = "Mehta", sep = sep)
    }

    #----------#
    # Matrix F
    #----------#
    if (covariance == TRUE) {
        matF1 <- mxMatrix(type = "Diag", nrow = ntv+ntc, ncol = ntv+ntc, values = 1, name = "F1")
        matF2 <- mxMatrix(type = "Full", nrow = ntv+ntc, ncol = ntACE, name = "F2")
    } else {
        matF1 <- mxMatrix(type = "Diag", nrow = ntv, ncol = ntv, values = 1, name = "F1")
        matF2 <- mxMatrix(type = "Full", nrow = ntv, ncol = ntACE, name = "F2")
    }

    matF <- mxAlgebra(expression = cbind(F1,F2), name = "Filter")


    #----------#
    # Matrix I
    #----------#
    matI <- mxMatrix(type = "Diag", nrow = ntotal, ncol = ntotal, values = 1, name = "I")

    #-------------#
    # Expected Means
    #-------------#
    if (!is.null(moderation) | (!is.null(covvars) & covariance == FALSE)) {
        modMean <- mxAlgebra(expression =M+t(effMCovFull), name = "modM")
        mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%modM), name = "expMean")
    } else {
        mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%M), name = "expMean")
    }

    #-------------#
    # Expected Covariance Matrix
    #-------------#
    # MZ
    covMZ <- mxAlgebra(expression = Filter%*%solve(I-A)%*%SMZ%*%t(solve(I-A))%*%t(Filter), name = "expCovMZ")
    # DZ
    covDZ <- mxAlgebra(expression = Filter%*%solve(I-A)%*%SDZ%*%t(solve(I-A))%*%t(Filter), name = "expCovDZ")

    #-------------#
    # Data Objects
    #-------------#
    if (length(unique(data$zyg)) != 2) {
        stop("Zygosity variable must be coded as follows: 1 = MZ, 2 = DZ. Please, recode the zygosity variable.")
    }

    mzData    <- subset(data, zyg==1, unique(c(acevars_wide, covvars_all, mod_modvars_long, mod_modvars_wide)))
    dzData    <- subset(data, zyg==2, unique(c(acevars_wide, covvars_all, mod_modvars_long, mod_modvars_wide)))
    dataMZ    <- mxData( observed=mzData, type="raw" )
    dataDZ    <- mxData( observed=dzData, type="raw" )

    #---------------------#
    # Expectation Objects
    #---------------------#
    selVars <- c(acevars_wide)
    if (covariance == TRUE) {
        selVars <- c(selVars, covvars_all)
    }
    if (is.null(ordinal)) {
        expMZ     <- mxExpectationNormal(covariance="expCovMZ", means="expMean",
                                         dimnames= c(selVars))
        expDZ     <- mxExpectationNormal(covariance="expCovDZ", means="expMean",
                                         dimnames= c(selVars))
    } else {
        expMZ     <- mxExpectationNormal(covariance="expCovMZ", means="expMean", thresholds = "threshMat", threshnames = ordinal_wide,
                                         dimnames= c(selVars))
        expDZ     <- mxExpectationNormal(covariance="expCovDZ", means="expMean", thresholds = "threshMat", threshnames = ordinal_wide,
                                         dimnames= c(selVars))
    }
    #---------------------#
    # Fit Function
    #---------------------#
    fitfun     <- mxFitFunctionML()

    #---------------------#
    # Parameter List
    #---------------------#

    pars <- c(matA, matC, matE, matB, M0, AMatGap, SMatGap, AMatLower, SMatUpper, SMatLower, SMatACE, matF1, matF2, matF, matI, latentmeans, fitfun) # hier alle Matrizen ohne freie Parameter
    def <- c(mean, AMat, matBFull, matACEFull) # hier alle Matrizen, die Definitionsvariablen enthalten
    MZobj <- c(SMat1MZ, SMatMZ, covMZ, dataMZ, expMZ) # MZ spezifische Objekte
    DZobj <- c(SMat1DZ, SMatDZ, covDZ, dataDZ, expDZ, matdzA, matdzC) # DZ spezifische Objekte

    if (!is.null(covvars) & covariance == TRUE) {
        if (ncl > 0) {
            pars <- c(pars, AMatCovLong, SMatCov)
        }
        if (ncw > 0) {
            pars <- c(pars, AMatCovWide, SMatCov, AMatCovGap)
        }
        pars <- c(pars, AMatCov, SMatCovGap1, SMatCovGap2)
        if (!is.null(moderation)) {
            pars <- c(pars, covmeans)
        }
    }

    if (!is.null(ordinal)) {
        pars <- c(pars, thresh)
        if (!is.null(binary)) {
            if (is.null(moderation)) {
                pars <- c(pars,AMat, matBFull, matACEFull)
            } else {
                pars <- c(pars, AMatNull, matBFullNull, matACEFullNull)
            }
        }
    }

    if (!is.null(covvars) & covariance == FALSE) {
        pars <- c(pars, pathCovM)
        def <- c(def, DefCovM, effCovM, effMeanCovFull)
    }

    if (!is.null(moderation)) {
        pars <- c(pars, pathModM, BMod1, BMod2, BMod3, BMod4, BMod5, AMod1, AMod2, AMod3, AMod4, AMod5, CMod1, CMod2, CMod3, CMod4, CMod5, EMod1, EMod2, EMod3, EMod4, EMod5)
        def <- c(def, DefModM, RegModM, effCovM, effMeanCovFull, aModFull1, aModFull2, cModFull1, cModFull2, eModFull1, eModFull2, bModFull1, bModFull2, Mod1Def1, Mod1Def2, Mod2Def1, Mod2Def2, Mod3Def1, Mod3Def2, Mod4Def1, Mod4Def2, Mod5Def1, Mod5Def2)
    }
    if (!is.null(moderation) | (!is.null(covvars) & covariance == FALSE)) {
        def <- c(def, modMean)
    }

    #---------------------#
    # Standardization
    #---------------------#
    if (is.null(covvars) | covariance == FALSE) {
        matIden <- mxMatrix(type = "Iden", nrow = nv, ncol = nv, name = "Iden")
        matF1Var <- mxMatrix(type = "Diag", nrow = nv, ncol = nv, values = 1, name = "F1Var")
        matF2Var <- mxMatrix(type = "Full", nrow = nv, ncol = nv, name = "F2Var")
        matFVar <- mxAlgebra(expression = cbind(F1Var,F2Var), name = "FVar")
        matVar <- mxAlgebra(expression = (vec2diag(FVar%*%diag2vec(expCovMZ))), name = "Var")
        matSD <- mxAlgebra( expression=solve(sqrt(Iden*Var)), name="iSD")
        matSD2 <- mxAlgebra( expression=(sqrt(Iden*Var)), name="iSD2")
        matStdA <- mxAlgebra(expression = iSD %*% AA, name = "stda")
        matStdC <- mxAlgebra(expression = iSD %*% AC, name = "stdc")
        matStdE <- mxAlgebra(expression = iSD %*% AE, name = "stde")
        matStdB <- mxAlgebra(expression = iSD %*% AB %*% iSD2, name = "stdb")
        std <- c(matIden, matF1Var, matF2Var, matFVar, matVar,matSD, matSD2, matStdA, matStdC, matStdE, matStdB)
    } else {
        std <- NULL
    }

    #---------------------#
    # Model Objects
    #---------------------#
    modelMZ   <- mxModel(pars, def, MZobj, fitfun, std, name="MZ")
    modelDZ   <- mxModel(pars, def, DZobj, fitfun, name="DZ")
    multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
    if (is.null(binary)) {
        modelACE  <- mxModel("Full ACE",pars, modelMZ, modelDZ, multi)
    } else {
        parsvarbin <- c(parsvarbin, matF1, matF2, matI)
        modelACE  <- mxModel("Full ACE",pars, modelMZ, modelDZ, multi, parsvarbin)
    }
    #---------------------#
    # Run Model
    #---------------------#
    set.seed(1)
    if (Optimizer == "SLSQP") {
        mxOption(NULL, 'Default optimizer', 'SLSQP')
    } else if (Optimizer == "NPSOL") {
        mxOption(NULL, 'Default optimizer', 'NPSOL')
    } else {
        mxOption(NULL, 'Default optimizer', 'CSOLNP')
    }

    if (TryHard == TRUE) {
        if (!is.null(ordinal)) {
            fitACE    <- mxTryHardOrdinal(modelACE, extraTries = Tries, exhaustive = exh)
        } else {
            fitACE    <- mxTryHard(modelACE, extraTries = Tries, exhaustive = exh)
        }

    } else {
        fitACE    <- mxRun(modelACE)
    }

    return(fitACE)

}

# Output function for twinflex
twinOutput <- function(fit, rounded = FALSE, cf = 0.95, n = 4) {
    outputlist <- list()
    Par <- names(fit@output[["estimate"]])
    Est <- fit@output[["estimate"]]
    SE <- fit@output[["standardErrors"]]
    rownames(SE) <- NULL
    names(Est) <- NULL
    zv <- qnorm(1- (1 - cf) / 2)

    outputUStd <- data.frame("Par" = Par, "Est" = Est,"SE" = SE)
    outputUStd <- outputUStd %>%
        mutate(zval = Est/SE) %>%
        mutate(pval = 2*pnorm(q=abs(zval), lower.tail=FALSE)) %>%
        mutate(lowerCI = Est - zv*SE) %>%
        mutate(upperCI = Est + zv*SE)

    if (rounded == TRUE) {
        outputUStd <- outputUStd %>%
            mutate_if(is.numeric, round, digits = n)
    }

    outputlist[["All parameters (unstandardized)"]] <- outputUStd

    # Estimated parameters
    pathAStd <- c(fit@submodels[["MZ"]]@algebras[["stda"]]$result)
    pathCStd <- c(fit@submodels[["MZ"]]@algebras[["stdc"]]$result)
    pathEStd <- c(fit@submodels[["MZ"]]@algebras[["stde"]]$result)
    pathBStd <- c(fit@submodels[["MZ"]]@algebras[["stdb"]]$result)
    pathoutputStd <- c(pathAStd, pathCStd, pathEStd, pathBStd)

    # Standard Errors
    matStdASE <- c(x = mxSE("MZ.stda", model = fit, silent = TRUE))
    matStdCSE <- c(x = mxSE("MZ.stdc", model = fit, silent = TRUE))
    matStdESE <- c(x = mxSE("MZ.stde", model = fit, silent = TRUE))
    matStdBSE <- c(x = mxSE("MZ.stdb", model = fit, silent = TRUE))
    pathSEoutputStd <- c(matStdASE, matStdCSE, matStdESE, matStdBSE)

    # Labels
    pathALabel <- c(fit@matrices[["AA"]]$labels)
    pathCLabel <- c(fit@matrices[["AC"]]$labels)
    pathELabel <- c(fit@matrices[["AE"]]$labels)
    pathBLabel <- c(fit@matrices[["AB"]]$labels)
    pathLaboutput <- c(pathALabel, pathCLabel, pathELabel, pathBLabel)

    # Free parameter
    pathAFree <- c(fit@matrices[["AA"]]$free)
    pathCFree <- c(fit@matrices[["AC"]]$free)
    pathEFree <- c(fit@matrices[["AE"]]$free)
    pathBFree <- c(fit@matrices[["AB"]]$free)
    pathFreeoutput <- c(pathAFree, pathCFree, pathEFree, pathBFree)

    outputStd <- na.omit(data.frame("Par" = pathLaboutput, "EstStd" = pathoutputStd, "SEStd" = pathSEoutputStd, "Free" = pathFreeoutput))
    outputStd <- outputStd %>%
        mutate(SEStd = ifelse(Free == FALSE, NA, SEStd)) %>%
        select(-Free)

    # z and p values
    outputStd <- outputStd %>%
        mutate(zvalStd = ifelse(!is.na(SEStd), EstStd/SEStd, NA)) %>%
        mutate(pvalStd = ifelse(!is.na(SEStd), 2*pnorm(q=abs(zvalStd), lower.tail=FALSE), NA)) %>%
        mutate(lowerCIStd = EstStd - zv*SEStd) %>%
        mutate(upperCIStd = EstStd + zv*SEStd)

    if (rounded == TRUE) {
        outputStd <- outputStd %>%
            mutate_if(is.numeric, round, digits = n)
    }
    outputlist[["All parameters (standardized)"]] <- outputStd
    return(outputlist)
}

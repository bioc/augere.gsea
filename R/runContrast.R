#' Differential gene set analyses with contrasts
#'
#' Analyze a RNA-seq dataset for differential expression across a collection of gene sets.
#' This requires access to the original values, unlike \code{\link{runPrecomputed}}.
#' 
#' @inheritParams augere.de::runVoom
#' @inheritParams runPrecomputed 
#' @param comparison Character vector of length no greater than 2.
#' For length 2, this specifies the groups to be compared, whereas for length 1, this specifies the covariate to test. 
#' See \code{\link[augere.de]{processSimpleComparisons}} for more details.
#'
#' Alternatively, a named character vector with no more than 2 unique names,
#' see \code{\link[augere.de]{processSimpleComparisons}} for more details.
#' 
#' Unlike \code{\link[augere.de]{runVoom}}, a list of vectors is not accepted. 
#' @param contrast String, function or vector specifying a custom contrast,
#' see \code{?\link[augere.de]{processCustomContrasts}} for more details.
#' 
#' Unlike \code{\link[augere.de]{runVoom}}, a list of contrasts is not accepted.
#' @param methods Character vector specifying the methods to run.
#' This should contain at least one of the following:
#' \itemize{
#' \item \code{"mroast"}, which \code{\link[limma]{mroast}} from the \pkg{limma} package.
#' This is a self-contained gene set test where the null hypothesis is that the genes in the set are not differentially expressed.
#' \item \code{"fry"}, which calls \code{\link[limma]{fry}} from the \pkg{limma} package.
#' This is also a self-contained gene set test that is a fast approximation of \code{mroast}.
#' \item \code{"camera"}, which calls \code{\link[limma]{camera}} from the \pkg{limma} package.
#' This is a competitive gene set test where the null hypothesis is that the genes in the set are not more DE than the genes outside the set.
#' \item \code{"romer"}, which calls \code{\link[limma]{romer}} from the \pkg{limma} package.
#' This is also a competitive test that takes a different approach to accounting for correlations between genes in the same set.
#' }
#'
#' @details
#' Some of the methods involve randomization, so for full reproducibility, users should call \code{\link{set.seed}} before running \code{runContrast}.
#'
#' Note that, even if the user does not call \code{\link{set.seed}},
#' \code{runContrast} will automatically insert \code{\link{set.seed}} statements into the Rmarkdown report 
#' prior to any GSEA functions that involve randomization.
#' Each \code{set.seed} call has a hard-coded seed to ensure that future compilation of the generated report will give the same result.
#' However, different calls to \code{runContrast} will use different (randomly selected) seeds to avoid systematic biases.
#' Thus, if full reproducibility of \code{runContrast} is required, users should set the seed themselves before calling \code{runContrast}.
#'
#' @return
#' A Rmarkdown report named \code{report.Rmd} is written inside \code{output.dir}. 
#' This contains all commmands used to reproduce the analysis. 
#'
#' If \code{dry.run=FALSE}, a list of \link[S4Vectors]{DataFrame}s is returned where each DataFrame contains the enrichment results for a \code{method}.
#' Each row corresponds to a gene set in \code{sets} and is named accordingly.
#' All DataFrames are guaranteed to have (at least) the following fields:
#' \itemize{
#' \item \code{NumGenes}, the number of genes in the set (after removing genes that were not tested).
#' \item \code{Direction}, the net direction of the change in expression within the set (either \code{"up"} or \code{"down"}).
#' \item \code{PValue}, the p-value for enrichment in each gene set.
#' \item \code{FDR}, the Benjamini-Hochberg-adjusted p-value. 
#' }
#' Additional fields may be present for specific methods.
#'
#' If \code{save.results=TRUE}, the results are saved in a \code{results} subdirectory of \code{output.dir}.
#'
#' If \code{dry.run=FALSE}, only the report is created, and \code{NULL} is returned.
#' 
#' @author Aaron Lun
#' @examples
#' x <- augere.de::loadExampleDataset()
#' all.genes <- rownames(x)
#'
#' sets <- list(
#'     A = sample(all.genes, 92),
#'     B = sample(all.genes, 212),
#'     C = sample(all.genes, 12),
#'     D = sample(all.genes, 38),
#'     E = sample(all.genes, 55)
#' )
#'
#' output.dir <- tempfile()
#' res <- runContrast(
#'     x,
#'     sets,
#'     group="dex", 
#'     comparison=c("trt", "untrt"),
#'     output=output.dir
#' )
#'
#' res
#' list.files(output.dir, recursive=TRUE)
#' 
#' @export
#' @import augere.core augere.de limma
runContrast <- function(
    x,
    sets, 
    groups, 
    comparison, 
    covariates = NULL,
    block = NULL, 
    subset.factor = NULL,
    subset.levels = NULL,
    subset.groups = TRUE,
    design = NULL,
    contrast = NULL, 
    dc.block = NULL, 
    robust = TRUE, 
    quality = TRUE,
    trend = FALSE,
    methods = c("fry", "camera"), 
    assay = 1,
    annotation = NULL,
    metadata = NULL,
    output.dir = "contrast", 
    author = NULL,
    dry.run = FALSE, 
    save.results = TRUE,
    suppress.plots = FALSE
) {
    restore.fun <- resetInputCache()
    on.exit(restore.fun(), after=FALSE, add=TRUE)

    if (is.null(author)) {
        author <- Sys.info()[["user"]]
    }

    dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
    fname <- file.path(output.dir, "report.Rmd")

    template <- system.file("templates", "contrast.Rmd", package="augere.gsea", mustWork=TRUE)
    lines <- readLines(template)
    parsed <- parseRmdTemplate(lines)

    parsed[["setup-data"]] <- processInputCommands(x, name="se")
    parsed[["setup-sets"]] <- processInputCommands(sets, name="sets")

    replacements <- list(
        ASSAY = deparseToString(assay),
        AUTHOR = paste(sprintf("  - %s", author), collapse="\n"),
        RAW_AUTHOR = deparseToString(as.list(author))
    )

    # Setting up the design and contrasts.
    if (!is.null(subset.factor)) {
        parsed[["subset-data"]] <- replacePlaceholders(parsed[["subset-data"]], list(SUBSET_FACTOR=deparseToString(subset.factor), SUBSET_LEVELS=deparseToString(subset.levels)))
    } else {
        parsed[["subset-data"]] <- NULL
    }

    if (!is.null(design) && !is.null(contrast)) {
        parsed[["create-design"]] <- processCustomDesignMatrix(design=design, se.name="se")
        contrast.info <- processCustomContrasts(contrast)
        replacements$FILTER_ARGS <- "design=design"
        parsed[["subset-group"]] <- NULL

    } else {
        parsed[["create-design"]] <- processSimpleDesignMatrix(groups=groups, block=block, covariates=covariates, se.name="se")
        contrast.info <- processSimpleComparisons(comparison)

        if (is.null(groups)) {
            replacements$FILTER_ARGS <- "design=design"
            parsed[["subset-group"]] <- NULL

        } else {
            replacements$FILTER_ARGS <- "group=model.data$group."

            group.levels <- NULL
            if (subset.groups) {
                group.levels <- findSubsetGroups(contrast.info)
            }
            if (!is.null(group.levels)) {
                parsed[["subset-group"]] <- replacePlaceholders(
                    parsed[["subset-group"]],
                    list(
                        GROUP_FACTOR=deparseToString(groups),
                        GROUP_LEVELS=deparseToString(group.levels)
                    )
                )
            } else {
                parsed[["subset-group"]] <- NULL
            }
        }
    }

    if (length(contrast.info) != 1) {
        stop("expected exactly one contrast/comparison")
    }
    parsed[["create-contrast"]] <- contrast.info[[1]]$commands
    replacements$CONTRAST_NAME <- deparse(contrast.info[[1]]$title)

    # Now fitting the linear model and EB shrinking.
    voom.args <- ""
    lm.args <- character(0)
    if (!is.null(dc.block)) {
        replacements$DUPCOR_BLOCK <- deparseToString(dc.block)
        voom.args <- c(voom.args, "block=dc.block")
        lm.args <- c("block=dc.block", "correlation=fit$correlation")
    } else {
        parsed[["duplicate-correlation-text"]] <- NULL
        parsed[["duplicate-correlation-block"]] <- NULL
    }

    if (quality) {
        voom.args <- c(voom.args, "sample.weights=TRUE")
    } else {
        parsed[["quality-text"]] <- NULL
    }

    eb.args <- character(0) 
    if (trend) {
        eb.args <- c(eb.args, "trend=TRUE")
    }

    if (robust) {
        eb.args <- c(eb.args, "robust=TRUE")
        replacements$EXTRA_EB_CAPT <- " with outliers marked in red"
    } else {
        replacements$EXTRA_EB_CAPT <- ""
        parsed[["robust-text"]] <- NULL
    }

    replacements$VOOM_ARGS <- paste(voom.args, collapse=", ")
    replacements$EB_OPTS <- paste(c("", eb.args), collapse=", ")

    # Setting up some of the outputs.
    if (!is.null(annotation)) {
        replacements$ANNO_FIELDS <- deparseToString(annotation)
    } else {
        parsed[["common-annotation"]] <- NULL
    }

    merge.metadata <- !is.null(metadata)
    if (merge.metadata) {
        replacements$COMMON_METADATA <- deparseToString(metadata)
    } else {
        parsed[["create-common-metadata"]] <- NULL
    }

    meta.cmds <- processContrastMetadata(contrast.info[[1]])
    meta.cmds[1] <- paste0("contrast=", meta.cmds[1])
    meta.cmds <- sprintf("%s%s", strrep(" ", 8), meta.cmds)
    meta.cmds[length(meta.cmds)] <- paste0(meta.cmds[length(meta.cmds)], ",")

    resultify <- function(copy) {
        copy[["diff-metadata"]] <- meta.cmds
        if (is.null(annotation)) {
            copy[["add-annotation"]] <- NULL
        }
        if (is.null(subset.factor)) {
            copy[["subset-metadata"]] <- NULL
        }
        if (!merge.metadata) {
            copy[["merge-metadata"]] <- NULL
        } else {
            copy[["no-merge-metadata"]] <- NULL
        }
        copy
    }

    # At last, we can get to finally running each of the methods. 
    methods <- match.arg(methods, c("mroast", "fry", "camera", "romer"), several.ok=TRUE) 
    save.chunk.names <- character(0)
    rseed <- function() sample(.Machine$integer.max, 1)

    if ("mroast" %in% methods || "fry" %in% methods) {
        # No need to include weights as they are picked up by getEAWP in limma::.lmEffects.
        extra.args <- c(eb.args, lm.args)
        extra.arg.string <- .stringify_extra_args(extra.args)

        if ("mroast" %in% methods) {
            save.name <- "save-mroast"
            save.chunk.names <- c(save.chunk.names, save.name)
            parsed$mroast <- replacePlaceholders(
                resultify(parsed$mroast),
                list(
                    SAVING_CHUNK_NAME=save.name,
                    MROAST_ARGS=extra.arg.string,
                    SEED=deparseToString(rseed())
                ),
                error=FALSE
            )
        } else {
            parsed$mroast <- NULL
        }

        if ("fry" %in% methods) {
            save.name <- "save-fry"
            save.chunk.names <- c(save.chunk.names, save.name)
            parsed$fry <- replacePlaceholders(
                resultify(parsed$fry),
                list(
                    SAVING_CHUNK_NAME=save.name,
                    FRY_ARGS=extra.arg.string
                ),
                error=FALSE
            )
        } else {
            parsed$fry <- NULL
        }
    } else {
        parsed$mroast <- NULL
        parsed$fry <- NULL
    }

    if ("camera" %in% methods) {
        # No need to include weights as they are picked up by getEAWP in limma::.lmEffects.
        extra.args <- character(0)
        if (trend) {
            extra.args <- c(extra.args, "trend.var=TRUE")
        }

        save.name <- "save-camera"
        save.chunk.names <- c(save.chunk.names, save.name)
        parsed$camera <- replacePlaceholders(
            resultify(parsed$camera),
            list(
                SAVING_CHUNK_NAME=save.name,
                CAMERA_ARGS=.stringify_extra_args(extra.args)
            ),
            error=FALSE
        )
    } else {
        parsed$camera <- NULL
    }

    if ("romer" %in% methods) {
        # romer only accepts array weights, so we'll just have to make do.
        extra.args <- character(0)
        if (quality) {
            extra.args <- c(extra.args, "array.weights=fit$targets$sample.weights")
        }
        if (!is.null(dc.block)) {
            extra.args <- c(extra.args, "block=dc.block", "correlation=fit$correlation")
        }

        save.name <- "save-romer"
        save.chunk.names <- c(save.chunk.names, save.name)
        parsed$romer <- replacePlaceholders(
            resultify(parsed$romer),
            list(
                SAVING_CHUNK_NAME=save.name,
                ROMER_ARGS=.stringify_extra_args(extra.args),
                SEED=deparseToString(rseed())
            ),
            error=FALSE
        )
    } else {
        parsed$romer <- NULL
    }

    contents <- replacePlaceholders(parsed, replacements)
    writeRmd(contents, file=fname)

    if (dry.run) {
        return(NULL)
    }

    if (save.results) {
        skip.chunks <- NULL
    } else {
        skip.chunks <- c("save-directory", save.chunk.names)
    }

    if (suppress.plots) {
        skip.chunks <- c(
            skip.chunks,
            "plot-norm-factors",
            "plot-sa"
        )
    }

    env <- new.env()
    compileReport(fname, skip.chunks=skip.chunks, env=env, suppress.plots=suppress.plots)
    env$all.results
}

.stringify_extra_args <- function(extra.args) {
    if (length(extra.args)) {
        paste(sprintf(",\n    %s", extra.args), collapse="")
    } else {
        ""
    }
}

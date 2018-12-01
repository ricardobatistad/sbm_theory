scriptFilenameTemplate <- "runscript_%s_%s_%d_%d_%.2f.sh"

scriptTemplate <-
"#!/bin/bash
#SBATCH -e output_%s_%s/%d/%d/%.2f/errs/slurm-%%A-%%a.err
#SBATCH -o output_%s_%s/%d/%d/%.2f/outs/slurm-%%A-%%a.out
#SBATCH --array=%s
#SBATCH -c 2
#SBATCH --mem=%dG
#SBATCH -p volfovskylab,statdept-low,herringlab-low,common
module load R/3.4.4
export IPYTHONDIR=/tmp
Rscript --no-save --no-restore prelims.R %s %s %d %d %d"


#' Generate array for slurm arra job
#'
#' Relies on triage job to count from 0 onwards. Meant for R. If python, need to deal
#' with fact that indices start at 0.
#' @param group_mode standard or overlap
#' @param n_inst single n
#' @param d_inst single d
#' @param member_param single p or s
#' @param samps how many samples will we process
#' @param group how many indices to process per task ID
#' @param step step to get around slurm's 50k constraint
#' @param skip how many graphs (presumably done) to skip
#' @return Set of scripts
#' @export
arrayGen <- function(group_mode, n_inst, d_inst, member_param, samps, group, step, skip){

  arrayStr <- ""


  if (group_mode == "std") {

    prameters_ndp <- dplyr::filter(graph_stats_100, n == n_inst, d == d_inst, p == member_param)

    for (lambda_inst in Lambda) {

      prameters_ndpl <- dplyr::filter(prameters_ndp, lambda == lambda_inst)

      for (sign_inst in Sign) {

        prameters_ndpls <- dplyr::filter(prameters_ndpl, sign == sign_inst)

        if (all(prameters_ndpls$validIdx)) {

          beg  <- prameters_ndpls$graphNum[1] + skip
          end  <- (beg - 1) + samps

          low_bound <- if (group == 1)  beg - step  else (beg - 1)/group
          up_bound  <- if (group == 1)  end - step  else end/group - 1

          if (low_bound == up_bound) {
            arrayStr <- sprintf("%s%d,", arrayStr, low_bound)
          } else {
            arrayStr <- sprintf("%s%d-", arrayStr, low_bound)
            arrayStr <- sprintf("%s%d,", arrayStr, up_bound)
          }

        }

      }

    }


  } else {

    prameters_nds <- dplyr::filter(parameters_overlap_100, n == n_inst, d == d_inst, s == member_param)

    C <- unique(prameters_nds$c)

    for (c_inst in C) {

      prameters_ndsc <- dplyr::filter(prameters_nds, c == c_inst)

      beg  <- prameters_ndsc$graphNum[1] + skip
      end  <- (beg - 1) + samps

      low_bound <- if (group == 1)  beg - step  else (beg - 1)/group
      up_bound  <- if (group == 1)  end - step  else end/group - 1

      if (low_bound == up_bound) {
        arrayStr <- sprintf("%s%d,", arrayStr, low_bound)
      } else {
        arrayStr <- sprintf("%s%d-", arrayStr, low_bound)
        arrayStr <- sprintf("%s%d,", arrayStr, up_bound)
      }

    }

  }

  return(substr(arrayStr, 1, nchar(arrayStr) - 1))  # removes last superfluous comma

}



#' Generate slurm scripts
#'
#' @param param_mode semi or non
#' @param group_mode std or overlap
#' @param n single value
#' @param d single value
#' @param p single value
#' @param memReq single value
#' @param arrayStr string of task indices
#' @param group size of agglomeration
#' @param step size of step
#' @export
script_generator <- function(param_mode, group_mode, HARVEST_NITER, n, d, member_param, memReq, arrayStr, group, step) {

  scriptName <- sprintf(scriptFilenameTemplate, param_mode, group_mode, n, d, member_param)

  script <- sprintf(scriptTemplate,
                    param_mode, group_mode, n, d, member_param,
                    param_mode, group_mode, n, d, member_param,
                    arrayStr,
                    memReq,
                    param_mode, group_mode, HARVEST_NITER, group, step)

  scriptPath <- paste("scripts", scriptName, sep = "/")
  write(script, file = scriptPath) #scripts folder

}






#' Batch generate all scripts
#'
#' @param param_mode semi or non
#' @param group_mode std or overlap
#' @param HARVEST_NITER number of harvest iterations
#' @param samps number of samps to process
#' @param group how many indices to process per task
#' @param step step to get around slurm's 50k constraint
#' @param skip number of graphs to skip
#' @return Set of scripts
#' @export
scripts_all <- function(param_mode, group_mode, HARVEST_NITER, samps = 100, group = 1, step = 0, skip = 0) {


  if (group_mode == "std") {

    i = 0
    for (n in N) {
      for (d in D) {
        for (p in P) {
          i = i + 1
          arrayStr <- arrayGen(group_mode, n, d, p, samps, group, step, skip)
          script_generator(param_mode, group_mode, HARVEST_NITER,
                           n, d, p, memReq[i], arrayStr, group, step)
        }
      }
    }


  } else {

    S_sans_1 <- S[!S %in% 1.00]

    i = 0
    for (d in D_1k) {
      for (s in S_sans_1) {
        i <- i + 1
        arrayStr <- arrayGen(group_mode, 1000, d, s, samps, group, step, skip)
        script_generator(param_mode, group_mode, HARVEST_NITER,
                         1000, d, s, memReq_overlap$mem[i], arrayStr, group, step)
      }
    }

    for (d in D_10k) {
      for (s in S_sans_1) {
        i <- i + 1
        arrayStr <- arrayGen(group_mode, 10000, d, s, samps, group, step, skip)
        script_generator(param_mode, group_mode, HARVEST_NITER,
                         10000, d, s, memReq_overlap$mem[i], arrayStr, group, step)
      }
    }

  }

}


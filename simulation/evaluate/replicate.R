# Perform one step of many replications
replicate_once <- function(methods_list, params, knowledge, truth_global) {
	r <- generate_data(params)
	dat <- r$data
	truth_combined <- list(.global = truth_global, .local = r$truth_local)
	# Also get oracle knowledge from data generation method
	knowledge[[".local"]] <- r$knowledge_local
	# iterate through intermediate helper functions and update knowledge
	for(i in methods_list[[".intermediate"]]) {
		knowledge[[".intermediate"]][[i]] <- get(i)(dat, knowledge = knowledge)
	}
	# iterate through all methods
	estimation <- list()
	for(m in names(methods_list)) {
		if(m == ".intermediate") next
		estimation[[m]] <- storage(
			output = get(methods_list[[m]])(dat, knowledge = knowledge),
			truth_combined = truth_combined
		)
	}
	estimation
}

# function for replicating the simulation process
# and store outputs
replication <- function(methods_list, params, truth_global, B, cl = NULL) {
	knowledge <- list()
	source("generate/knowledge_global.R", chdir = TRUE, local = TRUE)
	# Primary oracle knowledge based on fixed parameters
	knowledge[[".global"]] <- knowledge_global(params)
	if(is.null(cl)) {
		# no progress bar
		results_temp <- lapply(as.list(1:B), function(i) replicate_once(methods_list, params,  knowledge, truth_global))
	} else {
		library(parallel)
		clusterExport(cl, c("params", "knowledge", "truth_global"), envir = environment())
		library(pbapply)
		results_temp <- pblapply(as.list(1:B), function(i) replicate_once(methods_list, params, knowledge, truth_global), cl = cl)
		# results_temp <- parLapply(cl, as.list(1:B), function(i) replicate_once(methods_list, params, knowledge, truth_global))
	}
	results <- list()
	for(m in names(methods_list)) {
		results[[m]] <- Map(f = function(x) x[[m]], results_temp)
	}
	results
}

# function for evaluating a list of methods,
# evaluated by a list of evaluators
evaluate <- function(methods_list, evaluators_list, params, B, cl = NULL) {
	# find truth based on parameters. This is global in the sense that
	# it does not depend on replication.
	source("generate/truth_global.R", chdir = TRUE, local = TRUE)
	truth_global <- truth_global(params)
	# replicate the process and store
	replication <- replication(methods_list, params, truth_global, B, cl)
	# read all evaluators
	source("evaluate/evaluators.R", chdir = TRUE, local = TRUE)
	# evaluate on replicated output
	performace <- list()
	for(m in names(methods_list)) {
		if(m == ".intermediate") next
		performace[[m]] <- list()
		for(e in names(evaluators_list)) {
			performace[[m]][[e]] <- get(evaluators_list[[e]])(replication[[m]], truth_global)
		}
	}
	data.frame(performace)
}

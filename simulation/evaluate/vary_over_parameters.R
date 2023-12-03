vary_over_parameters <- function(vary_params, methods_list, evaluators_list, params, B = 5000, ncores = 1) {
	source("evaluate/replicate.R", chdir = TRUE, local = TRUE)
	# read all methods and set up cluster
	library(parallel)
	if(ncores == 1) {
		cl <- NULL
	} else {
		cl <- makeCluster(min(detectCores() - 1, ncores))
		clusterSetRNGStream(cl, 2023)
	}
	if(is.null(cl)) {
		source("utils/storage.R", chdir = TRUE, local = TRUE)
		source("generate/generate_data.R", chdir = TRUE, local = TRUE)
		for(m in names(methods_list)) {
			if(m == ".intermediate") {
				for(i in methods_list[[m]]) {
					source(paste0("methods/intermediate/", i, ".R"), chdir = TRUE, local = TRUE)
				}
			} else {
				source(paste0("methods/", methods_list[[m]], ".R"), chdir = TRUE, local = TRUE)
			}
		}
	} else {
		library(parallel)
		clusterEvalQ(cl, source("utils/storage.R", chdir = TRUE, local = TRUE))
		clusterEvalQ(cl, source("generate/generate_data.R", chdir = TRUE, local = TRUE))
		clusterExport(cl, "methods_list", envir = environment())
		clusterExport(cl, c("replicate_once"), envir = environment())
		clusterEvalQ(cl,
			for(m in names(methods_list)) {
				if(m == ".intermediate") {
					for(i in methods_list[[m]]) {
						source(paste0("methods/intermediate/", i, ".R"), chdir = TRUE, local = TRUE)
					}
				} else {
					source(paste0("methods/", methods_list[[m]], ".R"), chdir = TRUE, local = TRUE)
				}
			}
		)
	}
	# this will store the results for each variation of the parameters
	vary_length <- length(vary_params[[1]])
	if(vary_length == 0) {
		print(111)
		t <- evaluate(methods_list, evaluators_list, params, B, cl)
		if(!is.null(cl)) {
			stopCluster(cl)
		}
		return(t)
	}
	results_list <- list()
	library(progress)
	pb <- progress_bar$new(format = "  Progress [:bar] :percent in :elapsed. ETA: :eta", total = vary_length, clear = FALSE)
	pb$tick(0)
	for(iter_num in 1:vary_length) {
		for(v in names(vary_params)) {
			if(v %in% names(params)) {
				params[[v]] <- vary_params[[v]][[iter_num]]
			} else {
				print("Warning: variable not present...")
			}
		}
		t <- evaluate(methods_list, evaluators_list, params, B, cl)
		for(v in names(vary_params)) {
			t[[v]] <- vary_params[[v]][[iter_num]]
		}
		results_list[[iter_num]] <- t
		pb$tick()
	}
	if(!is.null(cl)) {
		stopCluster(cl)
	}
	results_list
}

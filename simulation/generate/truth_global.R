# generate truth from params that will be used in evaluation
# it is passed to evaluators
# in storage function, it is accesible by truth[[".global"]]
truth_global <- function(params) {
	params[["P"]]
}


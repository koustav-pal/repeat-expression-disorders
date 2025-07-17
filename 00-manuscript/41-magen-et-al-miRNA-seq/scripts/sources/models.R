

return_ridge_composite_score <- function(discovery_df, rcis, concat_col_name = "ridge", seed){
	set.seed(seed)
	model_fit <- cv.glmnet(x = as.matrix(discovery_df[,rcis]), y = discovery_df[,"status"], alpha = 0, family = "binomial")
	lambda <- model_fit$lambda.min
	coefs <- coef(model_fit, s = "lambda.min")
	intercept <- as.numeric(coefs[1])
	discovery_df[,paste(concat_col_name, "log_odds", sep = "_")] <- as.vector(as.matrix(discovery_df[,rcis]) %*% as.vector(coefs[-1])) + intercept
	discovery_df[,paste(concat_col_name, "prob", sep = "_")] <- plogis(discovery_df[,paste(concat_col_name, "log_odds", sep = "_")])
	return(discovery_df)
}


predict_ridge_composite_score_using_model <- function(predict_df, rcis, model, seed, concat_col_name){
	set.seed(seed)
  a_matrix <- as.matrix(predict_df[, rcis, drop = FALSE])
  # Get log-odds and probability from glmnet model
  predict_df[,paste(concat_col_name, "log_odds", sep = "_")] <- as.vector(
  	predict(model, newx = a_matrix, s = "lambda.min", type = "link"))
  predict_df[,paste(concat_col_name, "prob", sep = "_")] <- as.vector(
  	predict(model, newx = a_matrix, s = "lambda.min", type = "response"))
	return(predict_df)
}


return_iterated_ridge_composite_score <- function(discovery_df, rcis, iterate_n, concat_col_name = "ridge", seed_list = NA){
	if(length(seed_list) != iterate_n & all(is.na(seed_list))){
		seed_list <- sample.int(.Machine$integer.max, iterate_n)		
	}
	score_list <- lapply(seq(1, iterate_n), function(x){
		message(x)
		set.seed(seed_list[x])
		model_fit <- cv.glmnet(x = as.matrix(discovery_df[,rcis]), y = discovery_df[,"status"], alpha = 0, family = "binomial")
		lambda <- model_fit$lambda.min
		coefs <- coef(model_fit, s = "lambda.min")
		intercept <- as.numeric(coefs[1])
		ridge_log_odds <- as.vector(as.matrix(discovery_df[,rcis]) %*% as.vector(coefs[-1])) + intercept
		ridge_prob <- plogis(ridge_log_odds)
		return(list(log_odds = ridge_log_odds, prob = ridge_prob, seed = seed_list[x], model = model_fit))
	})
	discovery_df[,paste(concat_col_name, "log_odds", sep = "_")] <- rowMedians(do.call(cbind, lapply(score_list, function(x){
		x$log_odds
	})))
	discovery_df[,paste(concat_col_name, "prob", sep = "_")] <- rowMedians(do.call(cbind, lapply(score_list, function(x){
		x$prob
	})))
	seed_list <- do.call(c, lapply(score_list, function(x){
		x$seed
	}))
	models <- lapply(score_list, function(x){
		x$model
	})
	return(list(prediction = discovery_df, seed_list = seed_list, model_list = models))
}

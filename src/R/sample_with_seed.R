sample_with_seed <- function(seed_id, x) {
    set.seed(seed_id)
    res <- sample(x)
    return(res)
}
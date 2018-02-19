#' Add prior information to a data wrapper
#'
#' Note that the given data wrapper requires a trajectory and expression values
#' to have been added already.
#'
#' @param data_wrapper A data wrapper to extend upon.
#'
#' @export
#'
#' @importFrom testthat expect_true
#' @importFrom dynutils is_wrapper_with_trajectory is_wrapper_with_expression
add_prior_information_to_wrapper <- function(
  data_wrapper
) {
  # check data wrapper
  testthat::expect_true(dynutils::is_wrapper_with_trajectory(data_wrapper))
  testthat::expect_true(dynutils::is_wrapper_with_expression(data_wrapper))

  # compute prior information and add it to the wrapper
  data_wrapper$prior_information <-
    with(data_wrapper, generate_prior_information(
      milestone_ids,
      milestone_network,
      progressions,
      milestone_percentages,
      counts,
      feature_info,
      cell_info
    ))

  # return the updated data wrapper
  data_wrapper
}

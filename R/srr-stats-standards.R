#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G1.2} I attempted to add lifecycle statements with the
#'   lifecycle package but was getting ROxygen errors when doing
#'   `usethis::use_lifecycle()`. It would be nice to get some guidance on the
#'   best way to do this.
#' @srrstatsNA {G1.6} This software has no performance claims nowadays.
#' @srrstatsNA {G2.3b} It has nothing that parameters are strictly
#'   case-sensitive.
#' @srrstatsNA {G2.4c,G2.4d,G2.4e} It has nothing that should be character and
#'   factor in my functions.
#' @srrstatsNA {G2.5} There are no inputs of type 'factor' in this software.
#' @srrstatsNA {G2.9} I do not belive I am doing any type conversions that
#'   could result in lost data.
#' @srrstatsNA {G2.10} The software has not tabular inputs.
#' @srrstatsNA {G2.14c} In order to ensure the accuracy of the
#'   calculation, it does not replace missing ('NA') data.
#' @srrstatsNA {G3.1,G3.1a} The covariance calculation is not involved in the
#'   software.
#' @srrstatsNA {G4.0} This Software which does not export outputs to be written
#'   to local files.
#' @srrstatsNA {G5.2b} There is no single or unique trigger condition for a
#'   function.
#' @srrstatsNA {G5.3} As the probability of this happening is low and only one
#'   outcome exists, a test is not considered necessary
#' @srrstatsNA {G5.4,G5.4a,G5.4b,G5.4c} We double the skeleton of comparisons
#'   using binding frameworks.
#' @srrstatsNA {G5.8d} Input a Zero-length data which is meaningless.
#' @srrstatsNA {G5.9,G5.9a,G5.9b} The software was intended to demonstrate the
#'   computability of the algorithm and the entity framework, so we did not
#'   conduct interference tests.
#' @srrstatsNA {G5.10,G5.11,G5.11a,G5.12} I have not implemented any extended
#'   tests
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @noRd
NULL

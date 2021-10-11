#' @title 
#' Japanese Breast Cancer Data (JBC)
#' @description 
#' Data on breast cancer incidence in prefects in Japan across 5 time periods. 
#' 
#' @docType data
#' 
#' 
#' 
#' @format An object of class "data.frame" with 1040 rows and 9 variables:
#' \itemize{
#'     \item id: identifier for each centroid location.
#'     \item period: time period at which each centroid location was measured (there are 5 total time periods, coded).
#'     \item death: count of the number of observed incident breast cancer deaths in each centroid-time period.
#'     \item expdeath: expected number of incident breast cancer deaths in each centroid-time period (age-standardized).
#'     \item covar1: simulated covariate 1 (will be unpenalized in model).
#'     \item covar2: simulated covariate 2 (will be unpenalized in model).
#'     \item covar3: simulated covariate 3 (will be unpenalized in model).
#'     \item covar4: simulated covariate 4 (will be unpenalized in model).
#'     \item covar5: simulated covariate 5 (will be unpenalized in model).
#' }
#' 
#' @keywords datasets
#' 
#' @examples 
#' data(jbc)
"jbc"


#' @title 
#' Polygon borders
#' @description 
#' Polygon borders for Japanese Breast Cancer dataset (JBC)
#' 
#' @docType data
#' 
#' 
#' 
#' @format An object of class "data.frame" with 1460 rows and 2 variables:
#' \itemize{
#'     \item V1: x-coordinate of the vector containing the coordinates of the vertices of the polygon.
#'     \item V2: y-coordinate of the vector containing the coordinates of the vertices of the polygon.
#' }
#' 
#' @keywords datasets
#' 
#' @examples 
#' data(japan.poly2)
"japan.poly2"


#' @title 
#' Prefect borders
#' @description 
#' Japan prefect borders for Japanese Breast Cancer dataset (JBC)
#' 
#' @docType data
#' 
#' 
#' 
#' @format An object of class "data.frame" with 108 rows and 4 variables:
#' \itemize{
#'     \item x1: x-coordinate of points from which to draw.
#'     \item y1: y-coordinate of points from which to draw.
#'     \item x2: x-coordinate of points to which to draw.
#'     \item y2: y-coordinate of points to which to draw.
#' }
#' 
#' @keywords datasets
#' 
#' @examples 
#' data(japan.prefect2)
"japan.prefect2"




#' @title 
#' UTM coordinates
#' @description 
#' UTM coordinates for Japanese Breast Cancer dataset (JBC)
#' 
#' @docType data
#' 
#' 
#' 
#' @format An object of class "data.frame" with 208 rows and 3 variables:
#' \itemize{
#'     \item id: identifier for each centroid location.
#'     \item utmx: easting (x-coordinate).
#'     \item utmy: northing (y-coordinate).
#' }
#' 
#' @keywords datasets
#' 
#' @examples 
#' data(utmJapan)
"utmJapan"
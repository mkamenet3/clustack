#' @title 
#' Japanese Breast Cancer Data 2 (JBC2)
#' @description 
#' Data on breast cancer incidence in prefects in Japan across 5 time periods. A fake cluster has been placed in center 150, with a radius of < 15km, in periods 3 through 5. A relative risk of 2 was used to create this cluster.
#' 
#' @docType data
#' 
#' 
#' 
#' @format An object of class "data.frame" with 1040 rows and 9 variables:
#' \itemize{
#'     \item Observed: count of the number of observed incident breast cancer deaths in each centroid-time period.
#'     \item Expected: expected number of incident breast cancer deaths in each centroid-time period (age-standardized).
#'     \item Time: time period at which each centroid location was measured (there are 5 total time periods, factor).
#' }
#' 
#' @keywords datasets
#' 
#' @examples 
#' data(jbc2)
"jbc2"


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
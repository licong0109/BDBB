#' mixture
#'
#' Mixture is  simulated data, which simulated methylation level data of 200 sites under 20 conditions was generated from a mixture of four beta distributions, corresponding to three LCPs and the background.
#'The parameters of four beta distribution followed by three LCPs and the background are set to (1.5, 1.5), (2, 15), (25, 2), and (20, 15), respectively. Three patterns contain 15, 20, and 30 sites respectively, and there is no overlap. In terms of conditions, there are 10, 10, and 8
#'conditions respectively. The first pattern and the second pattern
#'have 3 conditions that overlap, and the second pattern and the
#'third pattern have 2 conditions that overlap. Then BDBB was
#'applied on the simulated data, and three patterns was
#'recognized.
#'
#'
#' @examples
#'   head(mixture)
"mixture"

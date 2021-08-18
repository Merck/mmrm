# Copyright (c) 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
#
# This file is part of the mmrm program.
#
# mmrm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Repeated CD4 counts data from AIDS clinical trial
#' 
#' A dataset containing repeated measures of CD4 counts and baseline covariates.
#' 
#' @format A data frame with 5036 rows and 6 variables
#' \describe{
#'   \item{id}{Subject ID}
#'   \item{trt}{Treatment}
#'   \item{age}{Age (years)}
#'   \item{gender}{Gender (1=M, 0=F)}
#'   \item{week}{Week}
#'   \item{logcd4}{log(CD4 count + 1)}
#' }
#' @source \url{https://content.sph.harvard.edu/fitzmaur/ala/cd4.txt}
"cd4"


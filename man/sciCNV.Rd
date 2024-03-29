% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sciCNV.R
\name{sciCNV}
\alias{sciCNV}
\title{sciCNV function}
\usage{
sciCNV(
  norm.mat,
  ave.ctrl,
  gen.Loc,
  No.test,
  sharpness,
  baseline_adj = FALSE,
  baseline = 0
)
}
\arguments{
\item{norm.mat:}{Matrix of normalized single-cell gene expression (e.g. normalized scRNA-seq). Formatted as follows:
1st column = gene name list.
1st Row = cell identifier.
Cells are arranged with test cells in leftward columns, followed by control cells (normal diploid cells of matching lineage) in rightward columns.}

\item{ave.ctrl:}{It is the average expression of cells per gene in control population.}

\item{gen.Loc:}{It is the matrix of assigning chromosome number to each gene sorted based on chromosme number, starts and ends.}

\item{No.test:}{number of test cells}

\item{sharpness:}{a variable that adjusts the resolution used for the sciCNV-curve calculation (by defining a moving window size over which gene expression values are averaged).
Default =1.0 (range approximately 0.6-1.4).
A lower sharpness can be used to offset data sparsity and provide more reliable detection of large CNVs.
A higher sharpness provides more resolution for the detection of smaller CNVs but is more susceptible to noise from data sparsity and thus requires greater data density.}

\item{baseline_adj:}{The baseline adjustment is only applied to test cells if TRUE is specified. Default is FALSE.}

\item{baseline:}{An optional correction to adjust the CNV zero setpoint (copy number gain =0) and improve CNV detection when CN gains and losses are substantially unbalanced.
Consider using for markedly hyperdiploid or hyodiploid cells with multiple trisomies or monosomies where the chromosome number deviates substantially from 46.
Default = 0. Ideal setting: a fraction representing the net fractional genomic change from diploidy (using a positive fraction for gain and a negative fraction for net genomic loss).
If the CNV profile is unknown, run the CNV analysis with default zero setting, review the preliminary CNV profile, and consider re-running with baseline ccorrection.}
}
\value{
Provides the iCNV profile of all single-cells at higher sinsitivity, eficiency and accuracy
}
\description{
Inferring copy number variations of scRNA-seq data at single-cell level
}
\examples{
sciCNV_mat <- sciCNV(norm.mat=normlized_data, ave.ctrl=mean_control, gen.Loc, No.test=100)

}
\author{
Ali Mahdipour-Shirayeh, Princess Margaret Cancer Centre, University of Toronto
}

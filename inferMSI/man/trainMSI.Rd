\name{trainMSI}
\alias{trainMSI}

\title{
MSI model training
}

\description{
Create a model for MSI inference using normal data and cancer samples with known MSI status.
}

\usage{
trainMSI(germline, somatic, msi, markers = NULL, train.prop = 0.5, validate.prop = 0.2, seed = 2308947)
}

\arguments{
  \item{germline}{Character vector containing paths to .repeatseq files
    for normal samples}
  \item{somatic}{Character vector containing paths to .repeatseq files for
    cancer samples}
  \item{msi}{Logical vector containing known MSI status for tumor
    samples -- 1:1 correspondence with files in somatic}
  \item{markers}{Named character vector of marker names where
    names(markers) is the repeatseq identifier
    (e.g. "" is the marker name and
    "2:51288512-51288554 2_21.5_2_51288512_51288554_86_51_48_0_0_1.0_AC"
    is the repeatseq identifier). See ... for an example.}
  \item{train.prop}{Proportion of data to be used in training. The
    remainder (i.e. 1 - train.prop) will be used to test the final model.}
  \item{validate.prop}{Proportion of training data to use in cross
    indivdiual validation steps during training (i.e. candidate model will
    be built using 1 - validate.prop of the training set and validated on
    validate.prop of the training set).}
  \item{seed}{Ingeger value to be pass on to set.seed()}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (normal, tumor)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

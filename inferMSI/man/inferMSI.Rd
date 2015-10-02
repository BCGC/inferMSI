\name{inferMSI}
\alias{inferMSI}

\title{Infer MSI status of a tumor sample}

\description{Infer MSI status of a tumor sample given a repeatseq file
  and a model of what normal samples look like.}

\usage{
inferMSI(file, model, cut.off)
}

\arguments{
  \item{file}{A repeatseq output file ending with ".repeatseq"}
  \item{model}{A model object returned from trainMSI()}
  \item{cut.off}{A p-value cutoff for classifying a sample as MSI}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{https://github.com/adaptivegenome/repeatseq}

\author{Randall Johnson}

\note{
  Example Snakemake pipelines for aligning reads and generating
  repeatseq files can be found at https://github.com/johnsonra/inferMSI/tree/master/Scripts
}

\examples{
# coming soon
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

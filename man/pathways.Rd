\name{pathways}
\alias{pathways}
\docType{data}
\title{Cancer related pathways }
\description{A list of nine cancer related pathways corresponding to the van t'Veer data.
Each element contains
a vector gene names corresponding to those in the data set.
}
\usage{data(pathways)}
\format{
The format is:  \cr
\code{List of 9}       \cr
 \code{$ androgen_receptor_signaling: chr [1:72] "AW025529" "NM_001648" "NM_001753" "NM_003298" ...}\cr
 \code{$ apoptosis                  : chr [1:187] "AB033060" "NM_002341" "NM_002342" "AI769763" ...}   \cr
 \code{$ cell_cycle_control         : chr [1:31] "NM_001759" "NM_001760" "NM_001786" "NM_001789" ...} \cr
 \code{$ notch_delta_signalling     : chr [1:34] "NM_002405" "AL133036" "NM_003260" "NM_004316" ...} \cr
 \code{$ p53_signalling             : chr [1:33] "NM_002307" "NM_002392" "NM_003352" "NM_002745" ...}\cr
 \code{$ ras_signalling             : chr [1:266] "D25274" "AI033397" "NM_003029" "NM_001626" ...}\cr
 \code{$ tgf_beta_signaling         : chr [1:82] "NM_003036" "AI090812" "AI697699" "AI760298" ...}   \cr
 \code{$ tight_junction_signaling   : chr [1:326] "D25274" "AA604213" "AF018081" "NM_003005" ...}  \cr
 \code{$ wnt_signaling              : chr [1:176] "AB033058" "AB033087" "NM_003012" "NM_003014" ...}
  }
\examples{
data(pathways)
#str(pathways)
}
\keyword{datasets}

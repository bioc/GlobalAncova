"GlobalAncova.1gene" <-
function(xx, group, perm=10000)
{
  t.test.result      <- t.test(xx ~ group)
	F.value            <- t.test.result$statistic^2
	p.value.theo       <- 1 - pf(F.value, 1, length(group)-2)

  perm.teststat      <- mt.sample.teststat(xx, group, test="f", B=perm)
  p.value.perm       <- length(perm.teststat[perm.teststat > F.value]) / length(perm.teststat)

	test.result        <- c(F.value, p.value.perm, p.value.theo)
	names(test.result) <- c("F.value", "p.value.perm", "p.value.theo")
  return(test.result)
}


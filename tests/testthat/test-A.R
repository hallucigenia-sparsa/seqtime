library(seqtime)
context("Interaction matrix")

test_that("A has desired connectance",{
  c=0.1
  A=generateA(N=50,c=c)
  cobs=getConnectance(A)
  expect_equal(c,cobs)
})

test_that("A has desired negative edge percentage",{
  nep=70
  A=generateA(N=50,c=0.1) # random assignment of edge strength results on average in 50%
  pep=getPep(A)
  A=modifyA(A,perc=nep,strength="uniform",mode="negpercent")
  pepAm=getPep(A)
  # will not get an exact match, but expect that I have now more negative edges
  print(pep)
  print(pepAm)
  expect_true(pepAm<pep)
})

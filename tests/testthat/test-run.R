context("test-run")

test_that("tumopp runs", {
  expect_s3_class(tumopp("-N8"), "data.frame")
  .const = c("-D3", "-Chex", "-k1", "-N8")
  .alt = list(
    L = c("const", "step"),
    P = c("random", "mindrag", "roulette")
  )
  expect_message({
    argslist = make_args(alt = .alt, const = .const, each = 2L)
  }, "invalid and excluded")
  expect_length(argslist, 2L * prod(lengths(.alt)) - 2L)
  expect_s3_class(tumopp(argslist[1:2], mc.cores = 1L), "data.frame")
})

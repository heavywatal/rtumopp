test_that("tumopp runs", {
  expect_s3_class(tumopp("-N8"), "data.frame")

  .const = c(D = 3, C = "hex", k = 1, N = 8)
  .alt = list(
    L = c("const", "step", "linear"),
    P = c("random", "mindrag", "roulette")
  )
  expect_message(
    {
      argslist = make_args(alt = .alt, const = .const, each = 2L)
    },
    "invalid and excluded"
  )
  expect_equal(nrow(argslist), 2L * prod(lengths(.alt)) - 4L)
  expect_s3_class(argslist, "data.frame")
  expect_s3_class(tumopp(head(argslist, 3L), mc.cores = 1L), "data.frame")

  .prior = list(
    L = function(n = 1L) {
      sample(c("const", "step", "linear"), n, replace = TRUE)
    },
    P = function(n = 1L) {
      sample(c("random", "mindrag"), n, replace = TRUE)
    },
    d = function(n = 1L) {
      runif(n, 0.0, 0.3)
    },
    m = function(n = 1L) {
      runif(n, 0.0, 2.0)
    }
  )
  .n = 4L
  expect_silent({
    argslist = generate_args(prior = .prior, const = .const, n = .n)
  })
  expect_equal(length(argslist), .n)
})

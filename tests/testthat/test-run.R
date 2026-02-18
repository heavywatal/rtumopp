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
  expect_identical(nrow(argslist), 2L * as.integer(prod(lengths(.alt))) - 4L)
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
  expect_length(argslist, .n)
})

test_that("help and version can be displayed", {
  tumopp_version() |>
    expect_silent() |>
    expect_type("character") |>
    expect_length(1L) |>
    expect_match(
      "^v?\\d+\\.\\d+([-.]\\d+)?(-\\d+)?(-g[a-fA-F0-9]{7,40})?(-dirty)?$"
    )
  expect_message(tumopp("--help"), "Usage")
})

test_that("sanitize_cache_dir works", {
  default = tempdir()
  expect_identical(sanitize_cache_dir(""), default)
  expect_identical(sanitize_cache_dir(character(0)), default)
  expect_identical(sanitize_cache_dir(NA_character_), default)
  expect_identical(sanitize_cache_dir(NA), default)
  expect_identical(sanitize_cache_dir(NULL), default)
  expect_identical(sanitize_cache_dir(FALSE), default)
  expect_identical(sanitize_cache_dir(TRUE), "~/.cache/tumopp")
  expect_identical(sanitize_cache_dir("path"), "path")
})

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("HMM works", {

  message('getting pseudobulk ..')

  bulk = get_bulk(
      count_mat = count_mat_test,
      df = df_test,
      lambdas_ref = ref_hca,
      gtf_transcript = gtf_transcript,
      genetic_map = genetic_map_hg38,
      min_depth = 0
  )

  message('running HMM..')

  bulk = bulk %>% analyze_bulk_lnpois(t = 1e-5)
  
})



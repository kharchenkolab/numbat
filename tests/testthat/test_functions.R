
library(numbat)

test_that("Check that roman2int works as expected", {
   val.3899 <- 'MMMDCCCXCIX'
   val.3900 <- 'MMMCM'
   val.4000 <- 'MMMM'
   expect_equal(names(roman2int(val.3899)), "MMMDCCCXCIX")
   expect_equal(as.integer(roman2int(val.3899)), 3899)
   expect_equal(names(roman2int(val.3900)), "MMMCM")
   expect_equal(as.integer(roman2int(val.3900)), 3900)
   expect_equal(names(roman2int(val.4000)), "MMMM")
   expect_equal(as.integer(roman2int(val.4000)), 4000)   
})

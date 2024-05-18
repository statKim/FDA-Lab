#include <Rcpp.h>
#include "file1.h"
// #include "test.h"
// #include "trapzRcpp.h"

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
int addTwoNumbers(int a, int b) {
  return add(a, b);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
addTwoNumbers(1, 3)
# trapzRcpp(1:10, 1:10)
*/

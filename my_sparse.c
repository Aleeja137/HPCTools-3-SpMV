#include "spmv.h"
#include <stddef.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>

int my_sparse(size_t n, gsl_spmatrix *mat, gsl_vector *vec, gsl_vector *result)
{
  size_t i;
  int j, row, col;
  int n_elem_row, n_elem_total = 0;

  for(i=0; i<n; i++)
  {
    // Get how many elements are in the row i
    n_elem_row = mat->p[i+1] - mat->p[i];

    // For each element, get their column index
    for (j=0; j<n_elem_row; j++)
    {
      row = i;
      col = mat->i[n_elem_total+j];

      // Do the multiplication and cumulative sum of result on target position
      double vec_element = gsl_vector_get(vec, col);
      double mat_element = gsl_spmatrix_get(mat, row, col);
      double result_element = gsl_vector_get(result, row);
      // This works since gsl_vector_calloc initializes vector to zero, otherwise an initial set whould be done
      result_element += mat_element * vec_element;
      gsl_vector_set(result, row, result_element);
    }
    n_elem_total += n_elem_row;
  }
  
  return 0;
}

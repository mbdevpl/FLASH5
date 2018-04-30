#include "cublas_v2.h"

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
dset_pointer(
    double **output_array,
    double *input,
    int lda,
    int row, int column,
    int batch_offset,
    int batchCount, cudaStream_t stream);

void 
ddisplace_pointers(
    double **output_array,
    double **input_array, int lda,
    int row, int column, 
    int batchCount, cudaStream_t stream);
#ifdef __cplusplus
}
#endif

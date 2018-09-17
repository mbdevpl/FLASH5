!===============================================================================
! INTERFACE to DBATCHED routines
!===============================================================================

MODULE dbatchedf

    PUBLIC :: dset_pointer, ddisplace_pointers
    PUBLIC :: doff1d, doff2d, doff3d

!-----------------------------------------------------------------------
! Function Interfaces
!-----------------------------------------------------------------------

    INTERFACE

        SUBROUTINE dset_pointer(output_array, input, lda, row, column, batch_offset, batchCount, stream) &
        &   BIND(C, NAME="dset_pointer")
            USE, INTRINSIC :: ISO_C_BINDING
            TYPE(C_PTR), VALUE :: output_array
            TYPE(C_PTR), VALUE :: input
            INTEGER(C_INT), VALUE :: lda
            INTEGER(C_INT), VALUE :: row
            INTEGER(C_INT), VALUE :: column
            INTEGER(C_INT), VALUE :: batch_offset
            INTEGER(C_INT), VALUE :: batchCount
            TYPE(C_PTR), VALUE    :: stream
        END SUBROUTINE dset_pointer

        SUBROUTINE ddisplace_pointers(output_array, input, lda, row, column, batchCount, stream) &
        &   BIND(C, NAME="ddisplace_pointers")
            USE, INTRINSIC :: ISO_C_BINDING
            TYPE(C_PTR), VALUE :: output_array
            TYPE(C_PTR), VALUE :: input
            INTEGER(C_INT), VALUE :: lda
            INTEGER(C_INT), VALUE :: row
            INTEGER(C_INT), VALUE :: column
            INTEGER(C_INT), VALUE :: batchCount
            TYPE(C_PTR), VALUE    :: stream
        END SUBROUTINE ddisplace_pointers

    END INTERFACE

CONTAINS

    SUBROUTINE doff1d( ptrNew, ptrOld, inc, i)
        USE, INTRINSIC :: ISO_C_BINDING
        TYPE(C_PTR) :: ptrNew
        TYPE(C_PTR) :: ptrOld
        INTEGER :: inc, i

        INTEGER(8) :: tmp_ptrNew

        tmp_ptrNew = TRANSFER( ptrOld, tmp_ptrNew )
        tmp_ptrNew = tmp_ptrNew + (i-1) * inc * sizeof(0.0d0)
        ptrNew = TRANSFER( tmp_ptrNew, ptrNew )

    END SUBROUTINE doff1d

    SUBROUTINE doff2d( ptrNew, ptrOld, lda, i, j)
        USE, INTRINSIC :: ISO_C_BINDING
        TYPE(C_PTR) :: ptrNew
        TYPE(C_PTR) :: ptrOld
        INTEGER :: lda, i, j

        INTEGER(8) :: tmp_ptrNew

        tmp_ptrNew = TRANSFER( ptrOld, tmp_ptrNew )
        tmp_ptrNew = tmp_ptrNew + ((j-1) * lda + (i-1)) * sizeof(0.0d0)
        ptrNew = TRANSFER( tmp_ptrNew, ptrNew )

    END SUBROUTINE doff2d

    SUBROUTINE doff3d( ptrNew, ptrOld, block, lda, i, j, k)
        USE, INTRINSIC :: ISO_C_BINDING
        TYPE(C_PTR) :: ptrNew
        TYPE(C_PTR) :: ptrOld
        INTEGER :: block, lda, i, j, k

        INTEGER(8) :: tmp_ptrNew

        tmp_ptrNew = TRANSFER( ptrOld, tmp_ptrNew )
        tmp_ptrNew = tmp_ptrNew + ((k-1) * block + (j-1) * lda + (i-1)) * sizeof(0.0d0)
        ptrNew = TRANSFER( tmp_ptrNew, ptrNew )

    END SUBROUTINE doff3d

END MODULE dbatchedf

///  @file parallel_matrix_op.c
///  @brief the standard implementations for functions in parallel_matrix_op.h
///
///  the standard implementations for functions in parallel_matrix_op.h
///

#include "parallel_matrix_op.h"
#include "blas.h"
#include "blas_comm.h"

///
/// @brief  Calculate the corresponding index in an array for an upper-triangle(UT) matrix.
///
/// @param[in]  i_row     - the i-th row in an upper-triangle matrix.
/// @param[in]  j_col     - the j-th column in an upper-triangle matrix.
/// @param[in]  dim       - the dimension of the upper-triangle matrix, i.e., an dim x dim matrix.
/// @return  the corresponding index in an array storage.
///
unsigned PQCLEAN_RAINBOWIACYCLIC_CLEAN_idx_of_trimat(unsigned i_row, unsigned j_col, unsigned dim) {
    return (dim + dim - i_row + 1) * i_row / 2 + j_col - i_row;
}

///
/// @brief  Calculate the corresponding index in an array for an upper-triangle or lower-triangle matrix.
///
/// @param[in]  i_row     - the i-th row in a triangle matrix.
/// @param[in]  j_col     - the j-th column in a triangle matrix.
/// @param[in]  dim       - the dimension of the triangle matrix, i.e., an dim x dim matrix.
/// @return    the corresponding index in an array storage.
///
static inline unsigned idx_of_2trimat(unsigned i_row, unsigned j_col, unsigned n_var) {
    if (i_row > j_col) {
        return PQCLEAN_RAINBOWIACYCLIC_CLEAN_idx_of_trimat(j_col, i_row, n_var);
    }
    return PQCLEAN_RAINBOWIACYCLIC_CLEAN_idx_of_trimat(i_row, j_col, n_var);
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_UpperTrianglize(unsigned char *btriC, const unsigned char *bA, unsigned Awidth, unsigned size_batch) {
    unsigned char *runningC = btriC;
    unsigned Aheight = Awidth;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < i; j++) {
            unsigned idx = PQCLEAN_RAINBOWIACYCLIC_CLEAN_idx_of_trimat(j, i, Aheight);
            gf256v_add(btriC + idx * size_batch, bA + size_batch * (i * Awidth + j), size_batch);
        }
        gf256v_add(runningC, bA + size_batch * (i * Awidth + i), size_batch * (Aheight - i));
        runningC += size_batch * (Aheight - i);
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_trimat_madd_gf16(unsigned char *bC, const unsigned char *btriA,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch) {
    unsigned Awidth = Bheight;
    unsigned Aheight = Awidth;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (k < i) {
                    continue;
                }
                gf16v_madd(bC, &btriA[(k - i) * size_batch], PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(&B[j * size_Bcolvec], k), size_batch);
            }
            bC += size_batch;
        }
        btriA += (Aheight - i) * size_batch;
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_trimatTr_madd_gf16(unsigned char *bC, const unsigned char *btriA,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch) {
    unsigned Aheight = Bheight;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (i < k) {
                    continue;
                }
                gf16v_madd(bC, &btriA[size_batch * (PQCLEAN_RAINBOWIACYCLIC_CLEAN_idx_of_trimat(k, i, Aheight))], PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(&B[j * size_Bcolvec], k), size_batch);
            }
            bC += size_batch;
        }
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_2trimat_madd_gf16(unsigned char *bC, const unsigned char *btriA,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch) {
    unsigned Aheight = Bheight;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                if (i == k) {
                    continue;
                }
                gf16v_madd(bC, &btriA[size_batch * (idx_of_2trimat(i, k, Aheight))], PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(&B[j * size_Bcolvec], k), size_batch);
            }
            bC += size_batch;
        }
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_matTr_madd_gf16(unsigned char *bC, const unsigned char *A_to_tr, unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
        const unsigned char *bB, unsigned Bwidth, unsigned size_batch) {
    unsigned Atr_height = Awidth;
    unsigned Atr_width = Aheight;
    for (unsigned i = 0; i < Atr_height; i++) {
        for (unsigned j = 0; j < Atr_width; j++) {
            gf16v_madd(bC, &bB[j * Bwidth * size_batch], PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(&A_to_tr[size_Acolvec * i], j), size_batch * Bwidth);
        }
        bC += size_batch * Bwidth;
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_bmatTr_madd_gf16(unsigned char *bC, const unsigned char *bA_to_tr, unsigned Awidth_before_tr,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch) {
    const unsigned char *bA = bA_to_tr;
    unsigned Aheight = Awidth_before_tr;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                gf16v_madd(bC, &bA[size_batch * (i + k * Aheight)], PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(&B[j * size_Bcolvec], k), size_batch);
            }
            bC += size_batch;
        }
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_mat_madd_gf16(unsigned char *bC, const unsigned char *bA, unsigned Aheight,
        const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch) {
    unsigned Awidth = Bheight;
    for (unsigned i = 0; i < Aheight; i++) {
        for (unsigned j = 0; j < Bwidth; j++) {
            for (unsigned k = 0; k < Bheight; k++) {
                gf16v_madd(bC, &bA[k * size_batch], PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(&B[j * size_Bcolvec], k), size_batch);
            }
            bC += size_batch;
        }
        bA += (Awidth) * size_batch;
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_quad_recmat_eval_gf16(unsigned char *z, const unsigned char *y, unsigned dim_y, const unsigned char *mat,
        const unsigned char *x, unsigned dim_x, unsigned size_batch) {
    unsigned char tmp[128];

    unsigned char _x[128];
    for (unsigned i = 0; i < dim_x; i++) {
        _x[i] = PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(x, i);
    }
    unsigned char _y[128];
    for (unsigned i = 0; i < dim_y; i++) {
        _y[i] = PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(y, i);
    }

    PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf256v_set_zero(z, size_batch);
    for (unsigned i = 0; i < dim_y; i++) {
        PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf256v_set_zero(tmp, size_batch);
        for (unsigned j = 0; j < dim_x; j++) {
            gf16v_madd(tmp, mat, _x[j], size_batch);
            mat += size_batch;
        }
        gf16v_madd(z, tmp, _y[i], size_batch);
    }
}

void PQCLEAN_RAINBOWIACYCLIC_CLEAN_batch_quad_trimat_eval_gf16(unsigned char *y, const unsigned char *trimat, const unsigned char *x, unsigned dim, unsigned size_batch) {
    unsigned char tmp[256];

    unsigned char _x[256];
    for (unsigned i = 0; i < dim; i++) {
        _x[i] = PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf16v_get_ele(x, i);
    }

    PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf256v_set_zero(y, size_batch);
    for (unsigned i = 0; i < dim; i++) {
        PQCLEAN_RAINBOWIACYCLIC_CLEAN_gf256v_set_zero(tmp, size_batch);
        for (unsigned j = i; j < dim; j++) {
            gf16v_madd(tmp, trimat, _x[j], size_batch);
            trimat += size_batch;
        }
        gf16v_madd(y, tmp, _x[i], size_batch);
    }
}

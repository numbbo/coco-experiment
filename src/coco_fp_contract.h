/**
 * @file coco_fp_contract.h
 * @brief Disable floating-point contraction for the whole amalgamation.
 *
 * Fusing a * b + c into a single FMA instruction changes results in the last bits. Compilers do this by
 * default on some platforms (e.g. clang on arm64), which would make benchmark function values depend on
 * compiler and platform. This file must be the first one included in the amalgamation so that the pragma
 * covers all COCO code.
 *
 * GCC ignores the STDC pragma, therefore builds using GCC must additionally pass -ffp-contract=off.
 */

#ifndef COCOEX_COCO_FP_CONTRACT_H
#define COCOEX_COCO_FP_CONTRACT_H

#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(_MSC_VER)
#pragma fp_contract(off)
#endif

#endif

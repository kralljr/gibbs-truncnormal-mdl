#ifndef __CSV_H__
#define __CSV_H__

#include <stddef.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

size_t csv_count_fields(const char *s);
const char *csv_read_double_fields(double *y, size_t n, const char *s);

__END_DECLS

#endif

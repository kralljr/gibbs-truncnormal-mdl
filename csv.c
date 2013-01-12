#include <stdlib.h>

/* You should never use these CSV parsing functions for anything else. */

size_t
csv_count_fields(const char *s)
{
    size_t i = 1;
    const char *c;

    for (c = s; *c != '\0'; c++) {
        if (*c == ',') {
            i++;
        }
    }

    return i;
}

const char *
csv_read_double_fields(double *y, size_t n, const char *s)
{
    const char *curptr = s;
    char *endptr = NULL;
    size_t i;

    /* Rely on strtod's behavior that it will stop reading when it gets to a
     * character that it doesn't recognize as part of a double, i.e., it should
     * stop when it gets to a comma.  As long as the string only contains parts
     * of doubles, whitespace, and comma's, we should be OK parsing it.
     */
    for (i = 0, curptr = s, endptr = NULL; i < n; i++) {
        y[i] = strtod(curptr, &endptr);

        if (curptr == endptr) {
            /* We got an error.  Return the current position so the caller can 
             * print it or whatever.
             */
            return curptr;
        } else if (*endptr == '\0' && i + 1 < n) {
            /* Looks like we're out of characters.  That's bad if we haven't 
             * read enough fields yet.
             */
            return curptr;
        } else {
            /* Skip the comma */
            curptr = endptr + 1;
        }
    }

    return NULL;
}


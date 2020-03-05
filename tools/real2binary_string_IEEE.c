#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#if defined(__LP64__) || defined(_LP64)
# define BUILD_64   1
#endif

/* constants for word and double-word size */
#define WDSZ 64
#define DWSZ 128

char *fmt_binstr (unsigned long n, unsigned char sz, unsigned char szs, char sep);
void show_fltmem (float f);
float xstrtof (char *str);
int main (int argc, char** argv) {

    if (argc < 2) {
        fprintf (stderr, "error: insufficient input. Usage: %s float\n", argv[0]);
        return 1;
    }
    float fvalue = xstrtof (argv[1]);
    show_fltmem (fvalue);
    return 0;
}

/** single-precision float in memory
 *  output the float, equivalent unsigned int, and
 *  binary representation of the number in memory
 */
void show_fltmem (float f)
{
    unsigned int i = *(unsigned int *)&f;
    printf ("%s\n", fmt_binstr (i, 32, 32, '-'));
}

/** returns pointer to formatted binary representation of 'n' zero padded to 'sz'.
 *  returns pointer to string contianing formatted binary representation of
 *  unsigned 64-bit (or less ) value zero padded to 'sz' digits with char
 *  'sep' placed every 'szs' digits. (e.g. 10001010 -> 1000-1010).
 */
char *fmt_binstr (unsigned long n, unsigned char sz, unsigned char szs, char sep) {
    static char s[DWSZ + 1] = {0};
    char *p = s + DWSZ;
    unsigned char i;
    for (i = 0; i < sz; i++) {
        p--;
 //       if (i > 0 && szs > 0 && i % szs == 0)
 //           *p-- = sep;
        *p = (n >> i & 1) ? '1' : '0';
    }
    return p;
}
/** string to float with error checking. */
float xstrtof (char *str)
{
    char *endptr = NULL;
    errno = 0;

    float val = strtof (str, &endptr);

    /* Check for various possible errors */
    if ((errno == ERANGE && (val == HUGE_VALF || val == HUGE_VALL)) ||
        (errno != 0 && val == 0)) {
        perror ("strtof");
        exit (EXIT_FAILURE);
    }

    if (endptr == str) {
        fprintf (stderr, "No digits were found\n");
        exit (EXIT_FAILURE);
    }
    return val;
}


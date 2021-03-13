/* The MIT License (MIT)

   Copyright (c) 2012-2015 Yecheng Fu <cofyc.jackson@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef ARGPARSE_H
#define ARGPARSE_H

/* For c++ compatibility */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

struct argparse;
struct argparse_option;

typedef int argparse_callback (struct argparse *self,
                               const struct argparse_option *option);

enum argparse_flag {
    ARGPARSE_STOP_AT_NON_OPTION = 1,
};

enum argparse_option_type {
    /* special */
    ARGPARSE_OPT_END,
    ARGPARSE_OPT_GROUP,
    /* options with no arguments */
    ARGPARSE_OPT_BOOLEAN,
    ARGPARSE_OPT_BIT,
    /* options with arguments (optional or required) */
    ARGPARSE_OPT_INTEGER,
    ARGPARSE_OPT_FLOAT,
    ARGPARSE_OPT_STRING,
};

enum argparse_option_flags {
    OPT_NONEG = 1,              /* disable negation */
};

/**
 *  argparse option
 *
 *  `type`:
 *    holds the type of the option, you must have an ARGPARSE_OPT_END last in your
 *    array.
 *
 *  `short_name`:
 *    the character to use as a short option name, '\0' if none.
 *
 *  `long_name`:
 *    the long option name, without the leading dash, NULL if none.
 *
 *  `value`:
 *    stores pointer to the value to be filled.
 *
 *  `help`:
 *    the short help message associated to what the option does.
 *    Must never be NULL (except for ARGPARSE_OPT_END).
 *
 *  `callback`:
 *    function is called when corresponding argument is parsed.
 *
 *  `data`:
 *    associated data. Callbacks can use it like they want.
 *
 *  `flags`:
 *    option flags.
 */
struct argparse_option {
    enum argparse_option_type type;
    const char short_name;
    const char *long_name;
    void *value;
    const char *help;
    argparse_callback *callback;
    intptr_t data;
    int flags;
};

struct argparse {
    // user supplied
    const struct argparse_option *options;
    const char *const *usages;
    int flags;
    const char *description;    // a description after usage
    const char *epilog;         // a description at the end
    // internal context
    int argc;
    const char **argv;
    const char **out;
    int cpidx;
    const char *optvalue;       // current option value
};

// built-in callbacks
int argparse_help_cb(struct argparse *self,
                     const struct argparse_option *option);

// built-in option macros
#define OPT_END()        { ARGPARSE_OPT_END, 0, NULL, NULL, 0, NULL, 0, 0 }
#define OPT_BOOLEAN(...) { ARGPARSE_OPT_BOOLEAN, __VA_ARGS__ }
#define OPT_BIT(...)     { ARGPARSE_OPT_BIT, __VA_ARGS__ }
#define OPT_INTEGER(...) { ARGPARSE_OPT_INTEGER, __VA_ARGS__ }
#define OPT_FLOAT(...)   { ARGPARSE_OPT_FLOAT, __VA_ARGS__ }
#define OPT_STRING(...)  { ARGPARSE_OPT_STRING, __VA_ARGS__ }
#define OPT_GROUP(h)     { ARGPARSE_OPT_GROUP, 0, NULL, NULL, h, NULL, 0, 0 }
#define OPT_HELP(show)   OPT_BOOLEAN('h', "help", NULL,                 \
                                     "show this help message and exit", \
                                     argparse_help_cb, show, OPT_NONEG)

int argparse_init(struct argparse *self, struct argparse_option *options,
                  const char *const *usages, int flags);
void argparse_describe(struct argparse *self, const char *description,
                       const char *epilog);
int argparse_parse(struct argparse *self, int argc, const char **argv);
void argparse_usage(struct argparse *self, bool print_description);

#ifdef __cplusplus
}
#endif

#endif

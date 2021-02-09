#ifndef _STDIO_H_
#include <stdio.h>
#endif

typedef struct
{
    size_t  count;
    char    **filename;
    FILE    **fp;
}   file_list_t;

void    usage(char *argv[]);
void    open_files(char *list_filename, file_list_t *file_list, char *mode);
void    build_matrix(file_list_t *file_list, char *matrix_file);

